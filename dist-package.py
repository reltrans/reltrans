import os
import sys
import subprocess
import shutil

import xspectrampoline_helpers

if sys.platform.startswith("linux"):
    TARGET = "linux"
    SHARED_EXT = "so"
elif sys.platform == "darwin":
    TARGET = "macos"
    SHARED_EXT = "dylib"
else:
    raise Exception(f"Unsupported platform: '{sys.platform}'")

BUILD_DIR = "build"
PACKAGE_NAME = "pyreltrans"
MODEL_NAME = "reltrans"

HEADAS = xspectrampoline_helpers.get_HEADAS()
HEADAS_INCLUDE = os.path.join(HEADAS, "include")


def get_ldflags():
    if TARGET == "linux":
        return ["-lm", "-lpthread"]
    elif TARGET == "macos":
        import site
        site_packages_path = site.getsitepackages()[0]
        return ["-lgfortran", f"-Wl,-rpath,{site_packages_path}/xspectrampoline/LibXSPEC_v6_35_1/lib/"]
    else:
        raise "unreachable"


def get_fflags() -> list[str]:
    if TARGET == "linux":
        return ["-shared", "-export-dynamic"]
    elif TARGET == "macos":
        return ["-dynamiclib"]
    else:
        raise "unreachable"


def cmd_fortran() -> list[str]:
    cmd = [
        os.environ.get("FC", "gfortran"),
        f"-I{BUILD_DIR}",
        f"-I{HEADAS_INCLUDE}",
        f"-I{BUILD_DIR}/cache",
        f"-J{BUILD_DIR}/cache",
        "-DHAVE_INLINE",
        "-fPIC",
        "-fno-automatic",
        "-fno-second-underscore",
        "-fno-omit-frame-pointer",
        # TODO: is this needed?
        "-fopenmp",
    ]
    cmd += get_fflags()
    return cmd


def compile(*source_files: str) -> str:
    if not os.path.isdir(BUILD_DIR):
        os.mkdir(BUILD_DIR)

    if not os.path.isdir(os.path.join(BUILD_DIR, "lib")):
        os.mkdir(os.path.join(BUILD_DIR, "lib"))

    dir_lib_model = os.path.join(BUILD_DIR, "lib", PACKAGE_NAME)

    if not os.path.isdir(dir_lib_model):
        os.mkdir(dir_lib_model)

    if os.path.isdir(os.path.join(BUILD_DIR, "cache")):
        shutil.rmtree(os.path.join(BUILD_DIR, "cache"))

    os.mkdir(os.path.join(BUILD_DIR, "cache"))

    output_path = os.path.join(dir_lib_model, f"lib{MODEL_NAME}.{SHARED_EXT}")

    cmd = cmd_fortran()
    cmd += source_files
    cmd += ["-o", output_path]
    cmd += get_ldflags()
    cmd += xspectrampoline_helpers.get_linker_flags(
        ["XSFunctions", "XSModel", "fftw3", "cfitsio"],
        # Use relative rpaths so that when it is installed into the Python
        # envrionment it can find the xspectrampoline libraries.
        rpath_relative=True,
    )

    subprocess.run(cmd, check=True)
    shutil.copyfile(
        "./xspec/lmodel_reltrans.dat", os.path.join(dir_lib_model, "lmodel.dat")
    )
    return output_path


if __name__ == "__main__":
    lib_path = compile("wrappers.f90")
    print(f"All targets compiled to path '{BUILD_DIR}'")

    artifact_dir = xspectrampoline_helpers.get_artifact_dir(PACKAGE_NAME)

    print(f"""
To install into your Python environment use:

    python3 -m pip install .

Reltrans can then be used through Python.
    """)
