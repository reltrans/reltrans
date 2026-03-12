# Contributing to reltrans

Welcome, and thanks for wanting to contribute! This guide explains the coding
style and conventions we use across the reltrans codebase.

The code has been written by several people over time, so there are places
where the style is a bit inconsistent. The goal of this document is to give us
a shared baseline going forward — both for new code and when cleaning up old
code.

If something is unclear or you have a suggestion, feel free to open an issue!

---

## Table of Contents

- [Golden rules](#golden-rules)
- [Python](#python)
  - [Formatting](#formatting)
  - [Naming things](#naming-things)
  - [Type hints](#type-hints)
  - [Docstrings](#docstrings)
  - [Imports](#imports)
  - [Calling Fortran from Python (ctypes)](#calling-fortran-from-python-ctypes)
- [Fortran](#fortran)
  - [Formatting](#fortran-formatting)
  - [Naming things](#fortran-naming)
  - [Comments](#fortran-comments)
- [Writing tests](#writing-tests)
- [Git workflow](#git-workflow)

---

## Golden rules

Before getting into the specifics, here are three principles that apply to
every part of the codebase:

1. **Be descriptive, not clever.** A name like `energy_min` is almost always
   better than `emin` or `x`. The next person reading the code (including
   future you) will thank you.

2. **Consistency beats personal preference.** Even if you prefer a different
   style, follow the conventions here so the codebase stays easy to navigate
   for everyone.

3. **Explain the *why*, not the *what*.** Comments and docstrings are most
   valuable when they explain decisions, not when they restate what the code
   already makes obvious.

---

## Python

### Formatting

We follow [PEP 8](https://peps.python.org/pep-0008/) and use `black` for
auto-formatting. The line length limit is **88 characters**.

- 4 spaces for indentation — no tabs.
- One blank line between methods inside a class; two blank lines between
  top-level definitions (functions, classes).
- Always end files with a newline.

Before committing, run:

```bash
black .
flake8 .
```

### Naming things

| What | Style | Example |
|------|-------|---------|
| Files and modules | `lower_snake_case` | `f2py_interface.py` |
| Functions | `lower_snake_case` | `reltrans_dcp`, `gen_wrap` |
| Classes | `UpperCamelCase` | `Reltrans`, `DcpParameters` |
| Module-level constants | `UPPER_SNAKE_CASE` | `MAX_PARAMS` |
| Local variables | `lower_snake_case` | `photon_flux`, `n_bins` |
| Private helpers | `_lower_snake_case` | `_wrap_call`, `_get_snapshot` |

**A note on ctypes pointer aliases** — there are currently two naming styles
in the codebase (`type_float_p` in `f2py_interface.py` vs `f_float` in
`test/wrapper.py`). Going forward, use the shorter `f_<type>` style:

```python
f_float  = ct.POINTER(ct.c_float)
f_int    = ct.POINTER(ct.c_int)
f_double = ct.POINTER(ct.c_double)
```

**Wrapper function names** should reflect the astrophysics model they wrap,
written in plain `lower_snake_case`. Avoid cryptic abbreviations:

```python
# use the full model name in lower_snake_case
reltrans_pl(ear, params)
reltrans_dcp(ear, params)
reltrans_dbl(ear, params)

# not like this — cryptic and inconsistently capitalised
wPL(...)
wDCp(...)
```

### Type hints

All public functions must have type hints on both arguments and the return
value. Use `np.ndarray` for NumPy arrays, and `typing` for anything more
complex:

```python
from typing import Callable
import numpy as np

def gen_wrap(
    ear: np.ndarray,
    params: np.ndarray,
    func: Callable[..., None],
) -> np.ndarray:
    ...
```

### Docstrings

Every public function, method, and class needs a docstring. We use
**Google style**, which is already present in `f2py_interface.py`:

```python
def reltrans_dcp(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    """
    Calculate the transfer function for Disc-Corona (DCp) geometry.

    Args:
        ear (np.ndarray): Energy bin edges in keV.
        params (np.ndarray): Model parameters as float32.

    Returns:
        np.ndarray: Calculated photon flux (float32).
    """
    return gen_wrap(ear, params, _w_dcp)
```

Private helpers (prefixed with `_`) don't need a full docstring, but a brief
one-liner helps a lot when the purpose isn't immediately obvious.

### Imports

Group imports in this order, with a blank line between each group:

```python
# 1. Standard library
import os
import dataclasses
from typing import Callable

# 2. Third-party
import ctypes as ct
import numpy as np

# 3. Local/project
from . import wrapper
```

Never use wildcard imports (`from module import *`) — they make it very hard
to see where things come from.

### Calling Fortran from Python (ctypes)

A few conventions to keep the Fortran interface consistent:

- **Define pointer aliases once**, at the top of the interface module, and
  reuse them everywhere:

  ```python
  f_float  = ct.POINTER(ct.c_float)
  f_int    = ct.POINTER(ct.c_int)
  f_double = ct.POINTER(ct.c_double)
  ```

- **Set `.argtypes` and `.restype` explicitly** for every Fortran symbol,
  right after loading the library. This prevents hard-to-debug type errors:

  ```python
  lib.tdreltransdcp_.argtypes = [f_float, f_int, f_float, f_int, f_float]
  lib.tdreltransdcp_.restype  = None
  ```

- **Don't repeat the ctypes boilerplate** in every wrapper. Use the shared
  `_wrap_call` helper and write a thin wrapper on top of it.

- **Parameter dataclasses** (like `DcpParameters`) should live in `wrapper.py`
  and provide a `.to_numpy_array() -> np.ndarray` method that returns `float32`.

---

## Fortran

### Formatting <a name="fortran-formatting"></a>

- Always use free-form Fortran (`.f90`). Never use fixed-form (`.f`).
- 4-space indentation, no tabs.
- Keep lines under **132 characters** (the Fortran standard maximum).
- Always use `implicit none` in every subroutine, function, and module.
- Declare all variables at the top of each program unit, grouped logically
  with a brief comment explaining what they are.

### Naming things <a name="fortran-naming"></a>

| What | Style | Example |
|------|-------|---------|
| Subroutines / functions | `lower_snake_case` | `td_reltrans_dcp` |
| Modules | `lower_snake_case` | `reltrans_utils` |
| Local variables | `lower_snake_case` | `energy_min`, `n_bins` |
| Loop counters | Descriptive when in long routines | `i_energy` instead of bare `i` |
| Constants (`parameter`) | `UPPER_SNAKE_CASE` | `PI`, `MAX_PARAMS` |
| Arrays | Descriptive, often plural | `photon_flux`, `energy_grid` |

### Comments <a name="fortran-comments"></a>

- Use `!` for all comments, with **a space after the `!`**:

  ```fortran
  ! This is easy to read.
  !This works but is harder to scan.
  ```

- Inline comments on variable declarations are encouraged — they serve as
  quick documentation without needing a separate block:

  ```fortran
  real    :: energy_min   ! lower energy bound of the grid [keV]
  real    :: energy_max   ! upper energy bound of the grid [keV]
  integer :: n_bins       ! number of energy bins
  ```

- **Don't leave commented-out dead code in the file.** Delete it and use
  `git` if you ever need to recover it.

---

## Writing tests

Tests live in `test/` and are run with `pytest`. Here are the conventions:

- Name test files `test_<feature>.py` (e.g. `test_basics.py`).
- Every test function must start with `test_` and have a name that clearly
  describes what it checks:

  ```python
  # the name should tell you what the test is checking
  def test_basic_invocation(reltrans, assert_snapshot): ...
  def test_re_im_returns_time_lag(reltrans, assert_snapshot): ...

  # these tell you nothing about what is being tested
  def test1(): ...
  def myTest(): ...
  ```

- Add a one-line docstring to each test. This shows up nicely in pytest output
  and makes reading test results much easier:

  ```python
  def test_basic_invocation(reltrans, assert_snapshot):
      """Smoke test: check that the default DCP parameters produce output."""
      ...
  ```

- Use the fixtures defined in `conftest.py` (`reltrans`, `assert_snapshot`,
  `envars`, `telescope`) instead of loading the library directly in your test.

- Try to test one concept per function. If you need to check multiple scenarios,
  use the `name=` argument to `assert_snapshot` to keep them separated and
  clearly labelled.

---

## Git workflow

### Branches

Pick a prefix that matches what you're doing:

| Work type | Branch name |
|-----------|-------------|
| New feature | `feature/<short-description>` |
| Bug fix | `fix/<short-description>` |
| Documentation | `docs/<short-description>` |
| Refactoring | `refactor/<short-description>` |

### Commit messages

We follow [Conventional Commits](https://www.conventionalcommits.org/). The
format is:

```
<type>(<scope>): <short summary, written as a command>

[Optional: a sentence or two explaining why, not what]

[Optional: Closes #issue-number]
```

Common types: `feat`, `fix`, `docs`, `test`, `refactor`, `chore`.

Some examples of good commit messages:

```
feat(wrapper): add reltrans_x model wrapper
fix(f2py_interface): correct argtypes for simrelt wrapper
docs: add CONTRIBUTING with code conventions (closes #37)
refactor(wrapper): standardise ctypes pointer alias names to f_<type>
```

### Pull requests

- Keep each PR focused on one logical change — it's much easier to review.
- Make sure all existing tests pass before opening the PR.
- Mention the relevant issue in the PR description (e.g. `Closes #37`) so
  GitHub links them automatically.
