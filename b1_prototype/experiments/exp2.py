import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1up import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
p = DCP_Parameters(boost=-1.0)
H.call(E, p); s = H.state()
# fourier_up self-test
th = 2*np.pi*(np.arange(16)+0.5)/16; thf = 2*np.pi*(np.arange(64)+0.5)/64
print(">> fourier_up err", np.max(np.abs(fourier_up(np.cos(3*th)+np.sin(th), 4) - (np.cos(3*thf)+np.sin(thf)))))
res = {}
def run(name, K, ntr):
    k0, k1 = K.arrays(); out = H.call(E, p, kernel=(k0, k1)); res[name] = (out, k0, ntr, K.ndeposits)
for cfg in [(32,3,4,4),(64,4,4,4),(64,4,8,8),(128,6,4,4),(128,6,8,8),(256,8,4,4),(256,8,8,8)]:
    t=time.time(); K, nr, nd = b1up_kernel(H, s, *cfg); run("b1up%dx%d_U%dx%d"%cfg, K, nr); print(">>", cfg, "%.1fs"%(time.time()-t))
ref = res["b1up256x8_U8x8"]
old = np.load("exp1.npz")
for k in old.files: res[k] = (old[k], None, 0, 0)
for k, (out, k0, nr, nd) in res.items():
    print(">> %-22s rays %8d deposits %9d  max|dS/S| %.2e  rms %.2e" % (k, nr, nd, np.max(np.abs(out/ref[0]-1)), np.sqrt(np.mean((out/ref[0]-1)**2))))
np.savez("exp2.npz", **{k: v[0] for k, v in res.items()})
