import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1 import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
p = DCP_Parameters(boost=-1.0)
H.call(E, p); s = H.state()
# flat mapping check
al = np.array([0., 300., 0., 500.]); be = np.array([-300*s["mu0"], 0., 400*s["mu0"], 0.])
print(">> flat re", H.flat(al, be))
res = {}
def run(name, K, ntr):
    k0, k1 = K.arrays(); out = H.call(E, p, kernel=(k0, k1)); res[name] = (out, k0, ntr, K.ndeposits)
t=time.time(); smp, nr = pixel_samples(H, s); run("pix200", kernel_from_samples(H, s, smp), nr)
smp, nr = pixel_samples(H, s, 800, 800); run("pix800", kernel_from_samples(H, s, smp), nr)
for dep in (deposit_nearest, deposit_cic):
    smp, nr, bad = b1_gl_samples(H, s, 64, 3); print(">> gl bad", bad)
    run("b1gl64x3_" + dep.__name__[8:], kernel_from_samples(H, s, smp, dep), nr)
for nth, ns in ((64, 4), (128, 8), (256, 16)):
    K, nr, nd = b1_tri_kernel(H, s, nth, ns); run(f"b1tri{nth}x{ns}", K, nr)
ref = res["b1tri256x16"]
for k, (out, k0, nr, nd) in res.items():
    e = np.max(np.abs(out / ref[0] - 1))
    kz = np.abs(k0).sum(0).sum(0); kzr = np.abs(ref[1]).sum(0).sum(0)
    print(">> %-22s rays %8d deposits %9d  max|dS/S| %.2e  zone-flux max rel diff %.2e" % (k, nr, nd, e, np.max(np.abs(kz / kzr - 1))))
np.savez("exp1.npz", **{k: v[0] for k, v in res.items()})
