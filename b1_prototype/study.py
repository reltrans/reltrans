"""Error-vs-cost study. Usage: python study.py GEOM FSET  -> results/GEOM_FSET.npz"""
import sys, os, time, json
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from methods import *; from b1 import *; from b1up import *; from ctf import *; from b1c import *

GEOMS = {
  "default":  dict(),
  "i70":      dict(inc=70.0),
  "a05h3":    dict(a=0.5, h=3.0),
  "i10h10r3": dict(inc=10.0, h=10.0, rin=-3.0),
  "a09i50h20":dict(a=0.9, inc=50.0, h=20.0),
}
FSETS = {
  "dc": [("dc", dict()), ("refl", dict(boost=-1.0))],
  "lo": [("re", dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=1.0)),
         ("im", dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=2.0)),
         ("lag", dict(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0))],
  "hi": [("re", dict(mass=4.6e7, flo_hz=5e-5, fhi_hz=1.5e-4, re_im=1.0)),
         ("im", dict(mass=4.6e7, flo_hz=5e-5, fhi_hz=1.5e-4, re_im=2.0)),
         ("lag", dict(mass=4.6e7, flo_hz=5e-5, fhi_hz=1.5e-4, re_im=4.0))],
}
def pix(n):
    def f(H, s):
        smp, nr = pixel_samples(H, s, n, n, max(50, n // 2), max(50, n // 2))
        return kernel_from_samples(H, s, smp, deposit_nearest), nr
    return f
def b1gl(nth, ngl, dep):
    def f(H, s):
        smp, nr, bad = b1_gl_samples(H, s, nth, ngl, nth, 8)
        return kernel_from_samples(H, s, smp, dep), nr
    return f
def b1tri(nth, nsub):
    def f(H, s):
        K, nr, nd = b1_tri_kernel(H, s, nth, nsub); return K, nr
    return f
def b1u(nth, ns, U, Us):
    def f(H, s):
        K, nr, nd = b1up_kernel(H, s, nth, ns, U, Us); return K, nr
    return f
def ctfm(nth, nr_, nchi, Ur):
    def f(H, s):
        K, nr, nd = ctf_kernel(H, s, nth, nr_, nchi, Ur, 1); return K, nr
    return f

METHODS = [("ref", "b1up_256x8_U8x8", b1u(256, 8, 8, 8))]
for n in (100, 141, 200, 283, 400, 566, 800, 1131, 1600):
    METHODS.append(("pix", f"pix_{n}", pix(n)))
for c in ((32, 2), (64, 3), (128, 4)):
    METHODS.append(("b1gl_cic", "b1gl_cic_%dx%d" % c, b1gl(*c, deposit_cic)))
for c in ((32, 2), (64, 4), (128, 8), (256, 16)):
    METHODS.append(("b1tri", "b1tri_%dx%d" % c, b1tri(*c)))
for c in ((24, 3, 4, 4), (32, 3, 4, 4), (32, 3, 8, 8), (48, 3, 8, 8), (64, 4, 4, 4), (64, 4, 8, 8), (96, 5, 8, 8), (128, 6, 8, 8)):
    METHODS.append(("b1up", "b1up_%dx%d_U%dx%d" % c, b1u(*c)))
for c in ((32, 3, 32, 4), (64, 3, 64, 4), (64, 4, 128, 8), (96, 5, 192, 8)):
    METHODS.append(("ctf", "ctf_%dx%dx%d_U%d" % c, ctfm(*c)))

def b1cm(nth, nc, ns, U, Us):
    def f(H, s):
        C = LineCache(H, s, nth, nc); K, nd = b1c_kernel(H, s, C, ns, U, Us); return K, C.ntrace
    return f
for c in ((32, 16, 3, 4, 4), (48, 24, 4, 4, 4), (64, 24, 4, 4, 4), (64, 32, 4, 8, 8), (96, 48, 4, 8, 8), (128, 48, 5, 8, 8)):
    METHODS.append(("b1c", "b1c_%dx%d_n%d_U%dx%d" % c, b1cm(*c)))

if __name__ == "__main__":
    geom, fset = sys.argv[1], sys.argv[2]
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
    os.makedirs("results", exist_ok=True)
    path = f"results/{geom}_{fset}.npz"
    res = dict(np.load(path, allow_pickle=True)) if os.path.exists(path) else {}
    params = [(fl, DCP_Parameters(**{**GEOMS[geom], **fp})) for fl, fp in FSETS[fset]]
    base = {fl: H.call(E, p) for fl, p in params}          # reltrans itself (pixel 200)
    for fl in base: res[f"reltrans/{fl}"] = base[fl]
    H.call(E, params[0][1]); s = H.state()
    for fam, name, fn in METHODS:
        if only and fam not in only and name not in only: continue
        if f"{name}/cost" in res: continue
        Harness.nq = 0; t0 = time.time()
        K, nrays = fn(H, s)
        tpy = time.time() - t0
        k0, k1 = K.arrays()
        for fl, p in params:
            res[f"{name}/{fl}"] = H.call(E, p, kernel=(k0, k1))
        res[f"{name}/cost"] = np.array([nrays, Harness.nq, K.ntri, K.ndeposits, s["nf"], tpy])
        if fam == "ref":
            res["ref/kernel_absW0"] = np.abs(k0).sum(axis=1)   # (NEX, xe)
        np.savez(path, **res)
        print(f"{geom} {fset} {name}: rays {nrays} nq {Harness.nq} ntri {K.ntri} ndep {K.ndeposits} ({tpy:.1f}s)", flush=True)
