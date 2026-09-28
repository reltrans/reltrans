import sys, time, copy, json
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from rtprof import *

telescope_env()
r = load()
E = np.logspace(-1, 2, 501)

flavours = {
  "time_avg":        dict(),
  "cross_lagE_nf3":  dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=1.0),
  "cross_lagE_nf12": dict(mass=10.0, flo_hz=0.5, fhi_hz=5.0, re_im=1.0),
  "lag_freq":        dict(mass=10.0, flo_hz=0.1, fhi_hz=10.0, re_im=7.0),
}
# Call types: new (a,i) -> trace; new h -> rtrans without trace; new logxi -> conv only; repeat -> cached
def call_seq(base):
    p0 = DCP_Parameters(**base)
    seq = []
    for k, (a, inc) in enumerate([(0.9, 35.0), (0.7, 50.0), (0.95, 20.0)]):
        p = copy.copy(p0); p.a = a; p.inc = inc; seq.append(("new_a_inc", p))
        q = copy.copy(p); q.h = 8.0 + k; seq.append(("new_h", q))
        q2 = copy.copy(q); q2.rin = -1.5; seq.append(("new_rin", q2))
        s = copy.copy(q2); s.logxi = 2.5 + 0.1*k; seq.append(("new_logxi", s))
        s2 = copy.copy(s); s2.gamma = 1.9 + 0.05*k; seq.append(("new_gamma", s2))
    return seq

# warm-up (table load, fftw plans)
t = time.perf_counter(); r.dcp(E, DCP_Parameters()); print("first call incl. table load: %.2f s" % (time.perf_counter()-t))
res = {}
for fl, base in flavours.items():
    # prime flavour (frequency grid change)
    r.dcp(E, DCP_Parameters(**base))
    agg = {}
    for kind, p in call_seq(base):
        reset(r)
        t = time.perf_counter(); r.dcp(E, p); wall = time.perf_counter() - t
        tm = timers(r)
        d = {n: tm[n][0] for n in NAMES}; d["wall"] = wall
        agg.setdefault(kind, []).append(d)
    res[fl] = {k: {n: float(np.median([d[n] for d in v])) for n in v[0]} for k, v in agg.items()}
json.dump(res, open("r0_profile.json", "w"), indent=1)
cols = ["wall", "total", "trace", "getdcos", "sum_gr", "sum_flat", "getlens", "restframe", "fftconv", "initcont", "post"]
for fl, kk in res.items():
    print("\n==", fl)
    print("%-10s" % "call" + "".join("%10s" % c for c in cols))
    for k, d in kk.items():
        print("%-10s" % k + "".join("%10.4f" % d[c] for c in cols))
