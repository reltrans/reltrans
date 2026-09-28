import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1up import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
d = dict(np.load("results/default_hi.npz"))
ps = [(fl, DCP_Parameters(mass=4.6e7, flo_hz=5e-5, fhi_hz=1.5e-4, re_im=r)) for fl, r in (("re",1.0),("im",2.0),("lag",4.0))]
H.call(E, ps[0][1]); s = H.state(); print(">> fi", s["fi"])
def run(K):
    out = {}
    k = K.arrays()
    for fl, p in ps: out[fl] = H.call(E, p, kernel=k)
    return out
def errs(a, b):
    return "  ".join("%s %.1e" % (fl, (np.max(np.abs(a[fl]-b[fl]))/np.max(np.abs(b[fl])))) for fl in a)
refs = {}
for (nto, nro) in ((128, 32), (256, 128), (512, 256)):
    import b1up as B; 
    K, nr, nd = b1up_kernel(H, s, 128, 6, 8, 8, nth_out=nto, nrho_out=nro); refs[(nto,nro)] = run(K)
base = {fl: d["b1up_256x8_U8x8/" + fl] for fl in ("re","im","lag")}
for k, v in refs.items(): print(">> b1up128x6 outer", k, "vs study ref:", errs(v, base), " vs 512x256 outer:", errs(v, refs[(512,256)]))
