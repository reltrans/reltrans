import re
import dataclasses
import textwrap
import pathlib
import sys
import argparse


@dataclasses.dataclass()
class SourceLine:
    line: str
    indent: int

    comment: bool = False
    has_comment: bool = False
    source: bool = False

    # Needs a continuation character
    cont: bool = False

    # Is a subroutine call
    call: Bool = False
    # Is a write call
    write: Bool = False
    # Is an assignment operation
    assignment: Bool = False
    is_open: Bool = False

    # Starts with `subroutine` or `function`
    function_declaration: Bool = False
    # Declares a variable
    declaration: Bool = False

    @staticmethod
    def from_string(line: str) -> "SourceLine":
        stripped = line.lstrip()
        indent = len(line) - len(stripped)
        stripped = stripped.rstrip()

        is_comment = re.match(r"^!.*", stripped) is not None

        has_comment = re.match(r".*!.*", stripped) is not None

        is_cont = re.match(r".*&$", stripped) is not None
        if not is_comment and is_cont:
            stripped = stripped.rstrip("&").rstrip()

        is_call = re.match(r"^call .*", stripped) is not None
        is_write = re.match(r"^write .*", stripped) is not None
        is_assignment = re.match(r"^.*=.*", stripped) is not None
        is_open = re.match(r"^open\s.*", stripped) is not None
        is_function_declaration = (
            re.match(r"^(subroutine|function)\s.*", stripped) is not None
        )
        is_var_declaration = re.match(r"^.*::.*", stripped) is not None

        return SourceLine(
            line=stripped,
            comment=is_comment,
            has_comment=is_comment or has_comment,
            source=not is_comment,
            indent=indent,
            cont=is_cont,
            call=is_call,
            write=is_write,
            assignment=is_assignment and not is_open,
            is_open=is_open,
            function_declaration=is_function_declaration,
            declaration=is_var_declaration,
        )

    def __len__(self) -> int:
        return len(self.line) + self.indent

    def startswith(self, s: str) -> bool:
        """
        Same as `str.startswith` but strips any leading spaces for the
        comparison.
        """
        return self.line.strip().startswith(s)

    def source_wrappable(self) -> bool:
        """
        Is this a line of code that can be line wrapped?
        """
        return self.call or self.function_declaration or self.declaration

    def word_wrap(self, width) -> list["SourceLine"]:
        return [
            dataclasses.replace(self, line=wrapped)
            for wrapped in textwrap.wrap(
                self.line,
                width=width - self.indent - 2,
            )
        ]

    def maths_wrap(
        self, width, maths_tokens=["\\*\\*", "\\=", "/", ",", "\\+", "\\-", "\\*"], indent_width=4
    ) -> "SourceLine":
        kwargs = dict(cont=True, indent=self.indent)

        all_tokens = [re.split(re.compile("(" + "|".join(maths_tokens) + ")"), self.line)]
        all_tokens = [i for i in iter(*all_tokens)]

        lines = []
        current_line = all_tokens[0].strip()
        spacer = True

        i = 1
        while i < len(all_tokens):
            token = all_tokens[i].strip()
            if len(token) == 0:
                i += 1
                continue

            # Lookaheads
            if i + 1 < len(all_tokens):
                if all_tokens[i + 1].strip() == ",":
                    token += ","
                    i += 1


            if len(current_line) + len(token) + self.indent + indent_width + 1 >= width:
                lines.append(dataclasses.replace(self, line=current_line.strip(), **kwargs))
                if len(lines) == 1:
                    kwargs["indent"] += indent_width
                current_line = ""

            if token == "**":
                current_line += token
                spacer = False
            elif token == "-":
                if current_line.endswith("(") or any(current_line.endswith(op) for op in maths_tokens):
                    current_line += token
                    spacer = False
                else:
                    current_line += " " + token
            elif spacer == False:
                current_line += token
                spacer = True
            else:
                current_line += " " + token

            i += 1

        lines.append(dataclasses.replace(self, line=current_line.strip(), **kwargs))
        return lines

    def source_wrap(self, width, sep=",", indent_width=4) -> "SourceLine":
        tokens = self.line.split(sep)

        def _add_sep(new: str, i: int) -> str:
            if (i > 0) and (i != len(tokens)) and len(new) > 0:
                return new + sep
            return new

        kwargs = dict(cont=True, indent=self.indent)
        lines = []
        new = ""
        for i, t in enumerate(tokens):
            if len(new) + len(t) + kwargs["indent"] + len(sep) + indent_width >= width:
                new = _add_sep(new, i)
                if new:
                    lines.append(dataclasses.replace(self, line=new, **kwargs))
                new = ""
                if len(lines) == 1:
                    kwargs["indent"] += indent_width
            new = _add_sep(new, i)
            new += t

        if new:
            lines.append(dataclasses.replace(self, line=new, **kwargs))

        lines[-1].cont = False

        return lines


def source_lines(lines: list[str]) -> list[SourceLine]:
    return [SourceLine.from_string(line) for line in lines]


class Formatter:
    def __init__(self, text_width=80):
        self.text_width = text_width

    def format_string(self, s: str) -> str:
        self.lines = source_lines(s.splitlines())

        self._merge_continuation()
        self._adjust_spacing()
        self._wrap_text()

        return self._reassemble()

    def _adjust_spacing(self):
        """
        Fix things like `hello(  world )` -> `hello(world)`, or `if(...)` -> `if (...)`
        """
        lines = []
        for line in self.lines:
            if line.comment:
                lines.append(line)
            else:
                new_text = line.line
                if line.declaration:
                    new_text = re.sub(r"\s*::\s*", " :: ", new_text)
                new_text = new_text.replace("if(", "if (")
                new_text = re.sub(r"\(\s+", "(", new_text)
                new_text = re.sub(r"\s+\)", ")", new_text)
                new_text = re.sub(r"\s\s+", " ", new_text)
                new_text = re.sub(r"\s*=(?!>)\s*", " = ", new_text)
                new_text = re.sub(r"\s*,\s*", ", ", new_text)
                new_text = re.sub(r"\s*(?!\*)\*\s*", " * ", new_text)
                lines.append(dataclasses.replace(line, line=new_text))
        self.lines = lines

    def _merge_continuation(self):
        lines = []
        for line in self.lines:
            if len(lines) > 0 and lines[-1].cont:
                lines[-1].line += line.line
                lines[-1].cont = line.cont
            else:
                lines.append(line)
        self.lines = lines

    def _reassemble(self) -> str:
        lines = []
        for line in self.lines:
            text = line.line
            if text == "":
                lines.append(text)
                continue
            text = (" " * line.indent) + text
            if line.cont:
                text = text.ljust(self.text_width - 1) + "&"
            lines.append(text)
        return "\n".join(lines)

    def _wrap_text(self):
        lines = []

        for l in self.lines:
            if len(l) > self.text_width:
                if l.comment:
                    # lines += l.word_wrap(self.text_width)
                    lines.append(l)
                elif not l.has_comment and l.source_wrappable():
                    lines += l.source_wrap(self.text_width)
                elif not l.has_comment and l.assignment:
                    lines += l.maths_wrap(self.text_width)
                    lines[-1].cont = False
                else:
                    lines.append(l)
            else:
                lines.append(l)

        self.lines = lines


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--inplace", help="Overwrite the file", action="store_true"
    )
    parser.add_argument("filename", help="The file to apply the formatter to.")
    args = parser.parse_args()

    path = pathlib.Path(args.filename)

    fmt = Formatter()
    text = path.read_text()
    formatted_text = fmt.format_string(text)

    if args.inplace:
        path.write_text(formatted_text)
        print(f"Formatted {path}")
    else:
        print(formatted_text)
