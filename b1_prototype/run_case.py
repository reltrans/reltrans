"""Evaluate a fixed set of reltrans outputs at the grid given by env vars; save npz.
usage: python run_case.py out.npz  (env RT_NRO etc. set by caller)"""
import sys, time, copy
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from rtprof import *
telescope_env()
r = load()
E = np.logspace(-1, 2, 501)          # as in test_basics
geoms = {
  "default": dict(),                                   # a=.998 i=30 h=6
  "i70":     dict(inc=70.0),
  "a05h3":   dict(a=0.5, h=3.0),
  "i10rin3": dict(inc=10.0, rin=-3.0, h=10.0),
}
flav = {
  "dc":      dict(),
  "re":      dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=1.0),
  "im":      dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=2.0),
  "lag":     dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0),
  "re_hf":   dict(mass=10.0, flo_hz=5.0, fhi_hz=20.0, re_im=5.0),
  "lag_hf":  dict(mass=10.0, flo_hz=5.0, fhi_hz=20.0, re_im=6.0),
  "refl":    dict(boost=-1.0),
}
out = {}
t0 = time.perf_counter()
for g, gp in geoms.items():
    for f, fp in flav.items():
        p = DCP_Parameters(**{**gp, **fp})
        out[f"{g}/{f}"] = r.dcp(E, p).astype(np.float64)
out["_time"] = np.array(time.perf_counter() - t0)
np.savez(sys.argv[1], **out)
