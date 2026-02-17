import os
import platform
import dataclasses
import ctypes as ct
import pathlib
import warnings

import numpy as np

f_double = ct.POINTER(ct.c_double)
f_float = ct.POINTER(ct.c_float)
f_int = ct.POINTER(ct.c_int)


def get_reltrans_library_path(lib_name="libreltrans") -> str:
    """
    Get the reltrans library path as a string. Checks common locations from the
    reltrans root directory. This can be overwritten using the `RELTRANS_PATH`
    environment variable.

    If no path is set, this function will return the name of the shared library
    with the hope that `LoadLibrary` will be able to resolve it in the linker
    path.
    """
    lib_path = os.environ.get("RELTRANS_PATH", None)
    if lib_path:
        warnings.warn(f"Using RELTRANS_PATH variable: {lib_path}")
        return lib_path

    build_dir = pathlib.Path(pathlib.Path(__file__).parent.parent) / "build" / "lib"

    system = platform.system()
    if system == "Linux":
        lib_name = lib_name + ".so"
    elif system == "Darwin":
        lib_name = lib_name + ".dylib"
    else:
        raise Exception("Unsupported OS " + system)

    lib_path = build_dir / lib_name
    print(str(lib_path.absolute()))
    if lib_path.is_file():
        return str(lib_path.absolute())

    return lib_name


def _wrap_call(f, energy: np.ndarray, params: np.ndarray) -> np.ndarray:
    ne = len(energy) - 1
    output = np.zeros(ne, dtype=np.float32)
    f(
        energy.ctypes.data_as(f_float),
        ct.byref(ct.c_int(ne)),
        params.ctypes.data_as(f_float),
        ct.byref(ct.c_int(1)),
        output.ctypes.data_as(f_float),
    )
    return output


@dataclasses.dataclass
class DCP_Parameters:
    # Lamp post height
    h: float = 6.0
    # Spin
    a: float = 0.998
    # Inclination (degrees)
    inc: float = 30.0
    # Inner radius
    rin: float = -1.0
    # Outer radius
    rout: float = 1e3
    # Cosmological redshift
    zcos: float = 0.0
    # Photon index
    gamma: float = 2.0
    # logξ ionisation parameter
    logxi: float = 3.0
    # Iron abundance
    afe: float = 1.0
    # Electron abundance
    lognep: float = 15.0
    # Electron temperature in observer frame
    kte: float = 60.0
    # Hydrogen column density
    nh: float = 0.0
    # Boosting factor (ad-hoc normalisation)
    boost: float = 1.0
    # Black hole mass in solar units
    mass: float = 4.6e7
    # Lowest frequency in band
    flo_hz: float = 0.0
    # Highest frequency in band
    fhi_hz: float = 0.0
    # 1 -> Re, 2 -> Im, 3 -> modulus, 4 -> time lag, 5 -> folded modulus, 6 -> folded time lag
    re_im: float = 1.0
    del_a: float = 0.0
    del_ab: float = 0.0
    g: float = 0.0
    telescope_response: float = 1.0

    def to_numpy_array(self) -> np.ndarray:
        return np.array(dataclasses.astuple(self), dtype=np.float32)


class Reltrans:
    def __init__(self, path=None):
        self.lib_reltrans = ct.cdll.LoadLibrary(path or get_reltrans_library_path())
        self.lib_reltrans.tdreltransdcp_.argtypes = [
            f_float,
            f_int,
            f_float,
            f_int,
            f_float,
        ]
        self.lib_reltrans.tdreltransdcp_.restype = None

    def dcp(self, energy: np.ndarray, parameters: DCP_Parameters) -> np.ndarray:
        """A wrapper around the XSPEC interface of reltransDcp"""
        return _wrap_call(
            self.lib_reltrans.tdreltransdcp_,
            energy.astype(np.float32),
            parameters.to_numpy_array(),
        )


__all__ = [DCP_Parameters, Reltrans, get_reltrans_library_path]
