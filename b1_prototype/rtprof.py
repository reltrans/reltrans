"""Profiling helpers for reltrans (timers added in rt_timing.f90)."""
import os, sys, time, ctypes as ct
import numpy as np
RT = os.environ.get("RELTRANS_ROOT", os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
BUILD = os.environ.get("RELTRANS_BUILD", "build")
LIBEXT = "dylib" if sys.platform == "darwin" else "so"
sys.path.insert(0, RT)
os.environ.setdefault("RELTRANS_TABLES", RT + "/cache/tables")
import pyreltrans
from pyreltrans import DCP_Parameters

NAMES = ["total", "rtrans", "trace", "getdcos", "sum_gr", "sum_flat", "getlens",
         "conv_total", "restframe", "fftconv", "initcont", "post", "grid"]

def load(**env):
    for k, v in env.items():
        os.environ[k] = str(v)
    r = pyreltrans.Reltrans(path=f"{RT}/{BUILD}/lib/libreltrans.{LIBEXT}")
    return r

def timers(r):
    out = (ct.c_double * 16)(); calls = (ct.c_int * 16)()
    r.lib_reltrans.rt_timing_get(out, calls)
    return {n: (out[i], calls[i]) for i, n in enumerate(NAMES)}

def reset(r):
    r.lib_reltrans.rt_timing_reset()

def telescope_env():
    ci = RT + "/cache/instrument-files/"
    os.environ["RMF_SET"] = ci + "nicer-rmf6s-teamonly-array50.rmf"
    os.environ["ARF_SET"] = ci + "nicer-consim135p-teamonly-array50.arf"
    os.environ["EMIN_REF"] = "0.3"; os.environ["EMAX_REF"] = "10.0"
