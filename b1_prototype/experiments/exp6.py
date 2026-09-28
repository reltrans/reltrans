import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1c import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
for geom in ("default", "i70"):
    d = dict(np.load(f"results/{geom}_dc.npz")); ref = d["b1up_256x8_U8x8/refl"]
    p = DCP_Parameters(boost=-1.0, **({} if geom == "default" else dict(inc=70.0))); H.call(E, p); s = H.state()
    C = LineCache(H, s, 64, 32)
    for (ns, U, Us) in ((4,4,4),(4,4,8),(4,2,8),(4,2,16),(4,1,16),(6,2,8),(6,4,8),(4,8,8),(4,4,16)):
        K, nd = b1c_kernel(H, s, C, ns, U, Us); out = H.call(E, p, kernel=K.arrays())
        print(">> %-8s ns %d Uth %2d Us %d: ntri %7d ndep %7d refl max %.2e rms %.2e" % (geom, ns, U, Us, K.ntri, K.ndeposits, np.max(np.abs(out/ref-1)), np.sqrt(np.mean((out/ref-1)**2))))
