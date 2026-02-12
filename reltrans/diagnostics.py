import os
import platform
import pathlib
import shutil


def main():
    print("Reltrans Build Diagnostics")
    print("---------------------------\n")

    # System information
    system = platform.system()
    arch = platform.machine()
    headas = os.environ.get("HEADAS")

    project_root = pathlib.Path(__file__).resolve().parent.parent
    lib_name = "libreltrans.so" if system == "Linux" else "libreltrans.dylib"
    expected_path = project_root / "build" / "lib" / lib_name

    print(f"Operating System : {system}")
    print(f"Architecture     : {arch}")
    print(f"HEADAS           : {headas if headas else 'NOT SET'}")
    print(f"Expected Library : {expected_path}")
    print(f"Library Exists   : {expected_path.exists()}")

    print("\nToolchain Checks:")

    # gfortran
    if shutil.which("gfortran"):
        print("gfortran         : Found")
    else:
        print("gfortran         : NOT FOUND")

    # xspec
    if shutil.which("xspec"):
        print("xspec            : Found")
    else:
        print("xspec            : NOT FOUND")

    # FFTW header check
    fftw_paths = [
        "/opt/homebrew/include/fftw3.f03",
        "/usr/local/include/fftw3.f03",
    ]

    fftw_found = any(pathlib.Path(p).exists() for p in fftw_paths)
    print("FFTW header      :", "Found" if fftw_found else "NOT FOUND")

    # Helpful guidance if library missing
    if not expected_path.exists():
        print("\nLibrary not found.")
        print("To build reltrans:")
        print("  1. Ensure HEASoft is installed and sourced (HEADAS set).")
        print("  2. Ensure gfortran and FFTW are installed.")
        print("  3. Run `make` in the reltrans root directory.")


if __name__ == "__main__":
    main()
