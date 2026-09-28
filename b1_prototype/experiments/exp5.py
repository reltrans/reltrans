import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1c import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
for geom in ("default", "i70"):
    d = dict(np.load(f"results/{geom}_dc.npz")); ref = d["b1up_256x8_U8x8/refl"]
    p = DCP_Parameters(boost=-1.0, **({} if geom == "default" else dict(inc=70.0))); H.call(E, p); s = H.state()
    for nth, nc in ((32, 16), (48, 24), (64, 24), (64, 32), (96, 48)):
        C = LineCache(H, s, nth, nc)
        for (ns, U, Us) in ((3, 4, 4), (4, 8, 8)):
            Harness.nq = 0; K, nd = b1c_kernel(H, s, C, ns, U, Us); out = H.call(E, p, kernel=K.arrays())
            print(">> %-8s cache %dx%d (rays %d, bad %d) ns %d U %d: ntri %d nq %d  refl max %.2e" % (geom, nth, nc, C.ntrace, C.bad, ns, U, K.ntri, Harness.nq, np.max(np.abs(out/ref-1))))
