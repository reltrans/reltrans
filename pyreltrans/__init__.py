import os
import platform
import dataclasses
import ctypes as ct
import pathlib
import warnings
import importlib.resources

import numpy as np

f_double = ct.POINTER(ct.c_double)
f_float = ct.POINTER(ct.c_float)
f_int = ct.POINTER(ct.c_int)


def get_reltrans_library_path(lib_name="libreltrans") -> str:
    """
    Get the reltrans library path as a string. Checks common locations from the
    (py)reltrans root directory. This can be overwritten using the
    `RELTRANS_PATH` environment variable.

    If no path is set, this function will return the name of the shared library
    with the hope that `LoadLibrary` will be able to resolve it in the linker
    path.
    """
    lib_path = os.environ.get("RELTRANS_PATH", None)
    if lib_path:
        warnings.warn(f"Using RELTRANS_PATH variable: {lib_path}")
        return lib_path

    # Is pyreltrans installed?
    try:
        import xspectrampoline_helpers as helpers
    except ModuleNotFoundError:
        warnings.warn(
            "Could not import `xspectrampoline_helpers`. Is `xspectrampoline` pip installed?"
        )
    else:
        _pyreltrans_dir = helpers.get_artifact_dir("pyreltrans")
        lib_path = _pyreltrans_dir / f"{lib_name}.{helpers.SHARED_LIB_EXT}"
        if lib_path.is_file():
            return str(lib_path.absolute())

    # Try to use what is in the `build` directory
    build_dir = pathlib.Path(pathlib.Path(__file__).parent.parent) / "build" / "lib"

    system = platform.system()
    if system == "Linux":
        lib_name = lib_name + ".so"
    elif system == "Darwin":
        lib_name = lib_name + ".dylib"
    else:
        raise Exception("Unsupported OS " + system)

    lib_path = build_dir / lib_name
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


def _wrap_trace_disk_observer(f, nphi, rn, mueff, mu0, spin, rmin, rout, mudisk, d):
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

def is_c_double(a):
    b = a if isinstance(a, ct.c_double) else ct.c_double(a)
    return b

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
    a: float = 0.997
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

@dataclasses.dataclass
class Simrelt_Parameters:
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
    flo_hz: float = 1e-5
    # Highest frequency in band
    fhi_hz: float = 1e-3
    # Squared coherence
    coherence_squared: float = 1.0
    del_A: float = 0.0
    del_AB: float = 0.0
    g: float = 0.0
    Anorm: float = 0.002
    exposure_time: float = 150000
    power: float = 10.0
    use_telescope_response: float = 1

    def to_numpy_array(self) -> np.ndarray:
        return np.array(dataclasses.astuple(self), dtype=np.float32)


class Reltrans:
    def __init__(self, path=None, **kwargs):
        self.lib_reltrans = ct.cdll.LoadLibrary(
            path or get_reltrans_library_path(**kwargs)
        )
        self.lib_reltrans.FNINIT()
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

        self.lib_reltrans.get_raytrace_coords.argtypes = [
            ct.POINTER(ct.c_double),  # p
            f_double,  # f1234(4)
            ct.POINTER(ct.c_double),  # lambda
            ct.POINTER(ct.c_double),  # q
            ct.POINTER(ct.c_double),  # sinobs
            ct.POINTER(ct.c_double),  # muobs
            ct.POINTER(ct.c_double),  # a_spin
            ct.POINTER(ct.c_double),  # robs
            ct.POINTER(ct.c_double),  # scal
            ct.POINTER(ct.c_double),  # radi
            ct.POINTER(ct.c_double),  # mu
            ct.POINTER(ct.c_double),  # phi
            ct.POINTER(ct.c_double),  # time
            ct.POINTER(ct.c_double),  # sigma
        ]
        self.lib_reltrans.get_raytrace_coords.restype = None

        self.lib_reltrans.trace_disk_observer.argtypes = [
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
        self.lib_reltrans.trace_disk_observer.restype = None

        self.lib_reltrans.get_needresp2.argtypes = [ct.POINTER(ct.c_int)]
        self.lib_reltrans.get_needresp2.restype = None

        self.lib_reltrans.test_kerrz_trace.argtypes = [
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
        ]

        self.lib_reltrans.test_kerrz_trace_lamppost.argtypes = [
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
        ]


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

    def simrelt(self, energy: np.ndarray, parameters: Dbl_Parameters) -> np.ndarray:
        return _wrap_call(
            self.lib_reltrans.simrelt_,
            energy.astype(np.float32),
            parameters.to_numpy_array(),
        )

    def trace_disk_observer(self, nphi, rn, mueff, mu0, spin, rmin, rout, mudisk, d):
        _wrap_trace_disk_observer(
            self.lib_reltrans.trace_disk_observer,
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

    def get_raytrace_coords(self,p,f1234,lambda_,q,sinobs,muobs,a_spin,
             robs, scal):
    # ctypes scalars
        f1234_c  = np.asarray(f1234, dtype=np.float64)
        p_c      = is_c_double(p      )
        lambda_c = is_c_double(lambda_)
        q_c      = is_c_double(q      )
        sinobs_c = is_c_double(sinobs )
        muobs_c  = is_c_double(muobs  )
        aspin_c  = is_c_double(a_spin )
        robs_c   = is_c_double(robs   )
        scal_c   = is_c_double(scal   )
        radi_c   = ct.c_double(0      ) #out
        mu_c     = ct.c_double(0      ) #out
        phi_c    = ct.c_double(0      ) #out
        time_c   = ct.c_double(0      ) #out
        sigma_c  = ct.c_double(0      ) #out
        self.lib_reltrans.get_raytrace_coords(
            ct.byref(p_c     ),
            f1234_c.ctypes.data_as(f_double),
            ct.byref(lambda_c),
            ct.byref(q_c     ),
            ct.byref(sinobs_c),
            ct.byref(muobs_c ),
            ct.byref(aspin_c ),
            ct.byref(robs_c  ),
            ct.byref(scal_c  ),
            ct.byref(radi_c  ),
            ct.byref(mu_c    ),
            ct.byref(phi_c   ),
            ct.byref(time_c  ),
            ct.byref(sigma_c )
            )
        return radi_c, mu_c, phi_c, time_c, sigma_c

    def p_disk_crossing(self,f1234,lambda_,q,sins,mus,a_spin,h,
          scal,mudisk,r_max,r_min):
    # ctypes scalars
        f1234_c  = np.asarray(f1234, dtype=np.float64)

        lambda_c = is_c_double(lambda_)
        q_c      = is_c_double(q      )
        sins_c   = is_c_double(sins   )
        mus_c    = is_c_double(mus    )
        aspin_c  = is_c_double(a_spin )
        h_c      = is_c_double(h      )
        scal_c   = is_c_double(scal   )
        mudisk_c = is_c_double(mudisk )
        rmax_c   = is_c_double(r_max  )
        rmin_c   = is_c_double(r_min  )
        out      = ct.c_double(0      ) #output
        self.lib_reltrans.p_disk_crossing(
            f1234_c.ctypes.data_as(f_double),
            ct.byref(lambda_c),
            ct.byref(q_c     ),
            ct.byref(sins_c  ),
            ct.byref(mus_c   ),
            ct.byref(aspin_c ),
            ct.byref(h_c     ),
            ct.byref(scal_c  ),
            ct.byref(mudisk_c),
            ct.byref(rmax_c  ),
            ct.byref(rmin_c  ),
            ct.byref(out     )
            )
        return out

    def p_coord_at_infinity(self,f1234,lambda_,q,sins,mus,a_spin,h,scal):
    # ctypes scalars
        f1234_c  = np.asarray(f1234, dtype=np.float64)
        lambda_c = is_c_double(lambda_)
        q_c      = is_c_double(q      )
        sins_c   = is_c_double(sins   )
        mus_c    = is_c_double(mus    )
        aspin_c  = is_c_double(a_spin )
        h_c      = is_c_double(h      )
        scal_c   = is_c_double(scal   )
    # vectors 
        out = np.zeros(4, dtype=np.float64)
        self.lib_reltrans.p_coord_at_infinity(
            f1234_c.ctypes.data_as(f_double),
            ct.byref(lambda_c),
            ct.byref(q_c     ),
            ct.byref(sins_c  ),
            ct.byref(mus_c   ),
            ct.byref(aspin_c ),
            ct.byref(h_c     ),
            ct.byref(scal_c  ),
            out.ctypes.data_as(f_double)
            )
        return out

    def constants_of_motion(self,alpha,beta,robs,sinobs,muobs,a_spin,scal):
    # ctypes scalars
        alpha_c  = is_c_double(alpha  )
        beta_c   = is_c_double(beta   )
        robs_c   = is_c_double(robs   )
        sinobs_c = is_c_double(sinobs )
        muobs_c  = is_c_double(muobs  )
        aspin_c  = is_c_double(a_spin )
        scal_c   = is_c_double(scal   )
        lambda_c = ct.c_double(0.0    ) #output
        q_c      = ct.c_double(0.0    ) #output
    # vectors
        # velocity_c = np.asarray(velocity, dtype=np.float64)
        velocity_c = np.zeros(3, dtype=np.float64)
        f1234_c  = np.zeros(4, dtype=np.float64) #output
        self.lib_reltrans.constants_of_motion(
            ct.byref(alpha_c ),
            ct.byref(beta_c  ),
            ct.byref(robs_c  ),
            ct.byref(sinobs_c),
            ct.byref(muobs_c ),
            ct.byref(aspin_c ),
            ct.byref(scal_c  ),
            velocity_c.ctypes.data_as(f_double),
            f1234_c.ctypes.data_as(f_double),
            ct.byref(lambda_c),
            ct.byref(q_c     )            
            )
        return f1234_c,lambda_c, q_c

    def initial_photon(self,pr,pt,pp,sins,mus,a_spin,h):
    # ctypes scalars
        alpha_c  = is_c_double(pr     )
        beta_c   = is_c_double(pt     )
        robs_c   = is_c_double(pp     )                                  
        sins_c   = is_c_double(sins   )
        mus_c    = is_c_double(mus    )
        aspin_c  = is_c_double(a_spin )
        h_c      = is_c_double(h      )
        lambda_c = ct.c_double(0      ) #output
        q_c      = ct.c_double(0      ) #output
    # vectors 
        velocity_c = np.zeros(3, dtype=np.float64)
        # velocity_c = np.asarray(velocity, dtype=np.float64)
        f1234_c  = np.zeros(4, dtype=np.float64) #output
        self.lib_reltrans.initial_photon(
            ct.byref(alpha_c ),
            ct.byref(beta_c  ),
            ct.byref(robs_c  ),
            ct.byref(sins_c  ),
            ct.byref(mus_c   ),
            ct.byref(aspin_c ),
            ct.byref(h_c     ),
            velocity_c.ctypes.data_as(f_double),
            ct.byref(lambda_c),
            ct.byref(q_c     ),
            f1234_c.ctypes.data_as(f_double)
            )
        return f1234_c,lambda_c,q_c

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

    def get_needresp2(self):
        flag = ct.c_int()
        self.lib_reltrans.get_needresp2(ct.byref(flag))
        return bool(flag.value)

    def kerrz_trace(self, spin, cos0, alpha, beta) -> tuple[float,float,float,float]:
        spin = ct.c_double(spin)
        cos0 = ct.c_double(cos0)
        alpha = ct.c_double(alpha)
        beta = ct.c_double(beta)
        t = ct.c_double(0.0)
        r = ct.c_double(0.0)
        th = ct.c_double(0.0)
        ph = ct.c_double(0.0)
        self.lib_reltrans.test_kerrz_trace(
            ct.byref(spin),
            ct.byref(cos0),
            ct.byref(alpha),
            ct.byref(beta),
            # Outputs:
            ct.byref(t),
            ct.byref(r),
            ct.byref(th),
            ct.byref(ph),
        )
        return (t.value, r.value, th.value, ph.value)

    def kerrz_trace_lamppost(self, spin, h, delta_s) -> tuple[float,float,float,float]:
        spin = ct.c_double(spin)
        h = ct.c_double(h)
        delta_s = ct.c_double(delta_s)
        t = ct.c_double(0.0)
        r = ct.c_double(0.0)
        th = ct.c_double(0.0)
        ph = ct.c_double(0.0)
        self.lib_reltrans.test_kerrz_trace_lamppost(
            ct.byref(spin),
            ct.byref(h),
            ct.byref(delta_s),
            # Outputs:
            ct.byref(t),
            ct.byref(r),
            ct.byref(th),
            ct.byref(ph),
        )
        return (t.value, r.value, th.value, ph.value)


__all__ = [
    DCP_Parameters,
    Reltrans,
    get_reltrans_library_path,
    Dbl_Parameters,
    rtdist_Parameters,
]
