import os
import platform
import pathlib


def main():
    print("Reltrans Build Diagnostics")
    print("------------------------------")

    system = platform.system()
    headas = os.environ.get("HEADAS")

    project_root = pathlib.Path(__file__).resolve().parent.parent
    lib_name = "libreltrans.so" if system == "Linux" else "libreltrans.dylib"
    expected_path = project_root / "build" / "lib" / lib_name

    print(f"Operating System : {system}")
    print(f"HEADAS           : {headas}")
    print(f"Expected Library : {expected_path}")
    print(f"Library Exists   : {expected_path.exists()}")

    if not expected_path.exists():
        print("\nLibrary not found.")
        print("To build:")
        print("  1. Ensure HEASOFT is installed and sourced.")
        print("  2. Run `make` in the reltrans root directory.")


if __name__ == "__main__":
    main()

