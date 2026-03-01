import subprocess
import os
import shutil
import xspectrampoline_helpers
from setuptools import Extension, setup

BUILD_DIR = "pybuild"
HEADAS = xspectrampoline_helpers.get_HEADAS()
HEADAS_INCLUDE = os.path.join(HEADAS, "include")


def linux_fflags():
    return ["-shared", "-export-dynamic"]


def linux_ldflags():
    return ["-lm", "-lpthread"]


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
    cmd += linux_fflags()
    return cmd


def compile(*source_files: str):
    if os.path.isdir(BUILD_DIR):
        shutil.rmtree(BUILD_DIR)
    os.mkdir(BUILD_DIR)
    os.mkdir(os.path.join(BUILD_DIR, "cache"))

    wrapper_o = os.path.join(BUILD_DIR, "wrappers.o")
    cmd = cmd_fortran()
    cmd += source_files
    cmd += ["-c", "-o", wrapper_o]

    print(cmd)
    subprocess.run(cmd, check=True)

    cmd = cmd_fortran()
    cmd += [wrapper_o]
    cmd += ["-o", os.path.join(BUILD_DIR, "pyreltrans")]
    cmd += linux_ldflags()
    cmd += xspectrampoline_helpers.get_linker_flags(
        ["XSFunctions", "XSModel", "fftw3", "cfitsio"],
        rpath_relative=True,
    )

    print(cmd)
    subprocess.run(cmd, check=True)


compile("wrappers.f90")
