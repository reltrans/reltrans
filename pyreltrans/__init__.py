import os
import platform
import dataclasses
import ctypes as ct
import pathlib
import warnings
import importlib.resources
import xspectrampoline_helpers as helpers

import numpy as np

f_double = ct.POINTER(ct.c_double)
f_float = ct.POINTER(ct.c_float)
f_int = ct.POINTER(ct.c_int)

_pyreltrans_dir = helpers.get_artifact_dir("pyreltrans")


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

    lib_path = _pyreltrans_dir / f"{lib_name}.{helpers.SHARED_LIB_EXT}"
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


def _wrap_getrgrid(f, rnmin, rnmax, mueff, nro, nphi):
    # convert scalars
    rnmin_c = ct.c_double(rnmin)
    rnmax_c = ct.c_double(rnmax)
    mueff_c = ct.c_double(mueff)
    nro_c = ct.c_int(nro)
    nphi_c = ct.c_int(nphi)
    # allocate outputs
    rn = np.zeros(nro, dtype=np.float64)
    domega = np.zeros(nro, dtype=np.float64)
    # call Fortran
    f(
        ct.byref(rnmin_c),
        ct.byref(rnmax_c),
        ct.byref(mueff_c),
        ct.byref(nro_c),
        ct.byref(nphi_c),
        rn.ctypes.data_as(f_double),
        domega.ctypes.data_as(f_double),
    )
    return rn, domega


def _wrap_grtrace(f, nphi, rn, mueff, mu0, spin, rmin, rout, mudisk, d):
    nro = len(rn)
    # ctypes scalars
    nro_c = ct.c_int(nro)
    nphi_c = ct.c_int(nphi)

    mueff_c = ct.c_double(mueff)
    mu0_c = ct.c_double(mu0)
    spin_c = ct.c_double(spin)
    rmin_c = ct.c_double(rmin)
    rout_c = ct.c_double(rout)
    mudisk_c = ct.c_double(mudisk)
    d_c = ct.c_double(d)
    rn = np.asarray(rn, dtype=np.float64)
    f(
        ct.byref(nro_c),
        ct.byref(nphi_c),
        rn.ctypes.data_as(f_double),
        ct.byref(mueff_c),
        ct.byref(mu0_c),
        ct.byref(spin_c),
        ct.byref(rmin_c),
        ct.byref(rout_c),
        ct.byref(mudisk_c),
        ct.byref(d_c),
    )


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


@dataclasses.dataclass
class Dbl_Parameters:
    # First Lamp post height
    h1: float = 6.0
    # Second Lamp post height
    h2: float = 50.0
    # Spin
    A: float = 0.997
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
    # Time-averaged normalization ratio C1 / C2 between the two lampposts and sets continuum cutoff and disk disk ionization
    eta_0: float = 1.0
    # Fourier frequency dependent normalization ratio C1(νc)/ C2(νc)
    eta: float = 1.0
    # propagation speed delay between the two lampposts if they are coherent (0 if incoherent)
    beta_p: float = 0.0
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
    del_ab1: float = 0.0
    g1: float = 0.0
    del_ab2: float = 0.0
    g2: float = 0.0
    telescope_response: float = 1.0

    def to_numpy_array(self) -> np.ndarray:
        return np.array(dataclasses.astuple(self), dtype=np.float32)


@dataclasses.dataclass
class rtdist_Parameters:
    # Lamp post height
    h: float = 6.0
    # Spin
    a: float = 0.996
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
    # distance of the soource
    Dist: float = 1e5
    # Iron abundance
    afe: float = 1.0
    # Electron abundance
    lognep: float = 15.0
    # Electron temperature in observer frame
    kte: float = 60.0
    # Hydrogen column density
    nh: float = 0.0
    # Boosting factor
    boost: float = 1.0
    # Black hole mass in solar units
    mass: float = 4.6e7
    # emissivity parameter 1
    honr: float = 0.02
    # emissivity parameter 1
    b1: float = 0.0
    # emissivity parameter 2
    b2: float = 0.0
    # Highest frequency in band
    flo_hz: float = 0.0
    # Highest frequency in band
    fhi_hz: float = 0.0
    # 1 -> Re, 2 -> Im, 3 -> modulus, 4 -> time lag, 5 -> folded modulus, 6 -> folded time lag
    re_im: float = -1.0
    del_a: float = 0.0
    del_ab: float = 0.0
    g: float = 0.0
    # Anorm
    Anorm: float = 7e-05
    telescope_response: float = 1.0

    def to_numpy_array(self) -> np.ndarray:
        return np.array(dataclasses.astuple(self), dtype=np.float32)


class Reltrans:
    def __init__(self, path=None, **kwargs):
        self.lib_reltrans = ct.cdll.LoadLibrary(
            path or get_reltrans_library_path(**kwargs)
        )
        self.lib_reltrans.tdreltransdcp_.argtypes = [
            f_float,
            f_int,
            f_float,
            f_int,
            f_float,
        ]
        self.lib_reltrans.tdreltransdcp_.restype = None

        self.lib_reltrans.getrgrid_.argtypes = [
            ct.POINTER(ct.c_double),  # rnmin
            ct.POINTER(ct.c_double),  # rnmax
            ct.POINTER(ct.c_double),  # mueff
            ct.POINTER(ct.c_int),  # nro
            ct.POINTER(ct.c_int),  # nphi
            f_double,  # rn
            f_double,  # domega
        ]
        self.lib_reltrans.getrgrid_.restype = None

        self.lib_reltrans.grtrace_.argtypes = [
            ct.POINTER(ct.c_int),  # nro
            ct.POINTER(ct.c_int),  # nphi
            f_double,  # rn
            ct.POINTER(ct.c_double),  # mueff
            ct.POINTER(ct.c_double),  # mu0
            ct.POINTER(ct.c_double),  # spin
            ct.POINTER(ct.c_double),  # rmin
            ct.POINTER(ct.c_double),  # rout
            ct.POINTER(ct.c_double),  # mudisk
            ct.POINTER(ct.c_double),  # d
        ]
        self.lib_reltrans.grtrace_.restype = None

    def dcp(self, energy: np.ndarray, parameters: DCP_Parameters) -> np.ndarray:
        """A wrapper around the XSPEC interface of reltransDcp"""
        return _wrap_call(
            self.lib_reltrans.tdreltransdcp_,
            energy.astype(np.float32),
            parameters.to_numpy_array(),
        )

    def dbl_lamp(self, energy: np.ndarray, parameters: Dbl_Parameters) -> np.ndarray:
        """A wrapper around the XSPEC interface of reltransDbl"""
        return _wrap_call(
            self.lib_reltrans.tdreltransdbl_,
            energy.astype(np.float32),
            parameters.to_numpy_array(),
        )

    def rtdist(self, energy: np.ndarray, parameters: rtdist_Parameters) -> np.ndarray:
        """A wrapper around the XSPEC interface of rtdist"""
        return _wrap_call(
            self.lib_reltrans.tdrtdist_,
            energy.astype(np.float32),
            parameters.to_numpy_array(),
        )

    def reset(self):
        # print(f"NAME OF THE LIBRARY {self.lib_reltrans._name}")
        self.lib_reltrans.reset_reltrans()

    def getrgrid(self, rnmin, rnmax, mueff, nro, nphi):
        return _wrap_getrgrid(
            self.lib_reltrans.getrgrid_, rnmin, rnmax, mueff, nro, nphi
        )

    def grtrace(self, nphi, rn, mueff, mu0, spin, rmin, rout, mudisk, d):
        _wrap_grtrace(
            self.lib_reltrans.grtrace_,
            nphi,
            rn,
            mueff,
            mu0,
            spin,
            rmin,
            rout,
            mudisk,
            d,
        )

    def get_re(self, nro, nphi):
        out = np.zeros(nro * nphi, dtype=np.float64)
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.get_re(
            out.ctypes.data_as(f_double),
            ct.byref(nro_c),
            ct.byref(nphi_c),
        )
        # reshape back to 2D
        # return out.reshape((nro, nphi), order="F")
        return out

    def set_re(self, nro, nphi):
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.allocate_re(ct.byref(nro_c), ct.byref(nphi_c))

    def get_taudo(self, nro, nphi):
        out = np.zeros(nro * nphi, dtype=np.float64)
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.get_taudo(
            out.ctypes.data_as(f_double),
            ct.byref(nro_c),
            ct.byref(nphi_c),
        )
        return out

    def set_taudo(self, nro, nphi):
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.allocate_taudo(ct.byref(nro_c), ct.byref(nphi_c))

    def get_pem(self, nro, nphi):
        out = np.zeros(nro * nphi, dtype=np.float64)
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.get_pem(
            out.ctypes.data_as(f_double),
            ct.byref(nro_c),
            ct.byref(nphi_c),
        )
        return out

    def set_pem(self, nro, nphi):
        nro_c = ct.c_int(nro)
        nphi_c = ct.c_int(nphi)
        self.lib_reltrans.allocate_pem(ct.byref(nro_c), ct.byref(nphi_c))


__all__ = [
    DCP_Parameters,
    Reltrans,
    get_reltrans_library_path,
    Dbl_Parameters,
    rtdist_Parameters,
]
