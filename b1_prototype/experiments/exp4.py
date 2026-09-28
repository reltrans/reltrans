import sys; sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), ".."))
from methods import *; from b1up import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
for gp in (dict(), dict(inc=70.0)):
    p = DCP_Parameters(boost=-1.0, **gp); H.call(E, p); s = H.state()
    K,_,_ = b1up_kernel(H, s, 256, 8, 8, 8); ref = H.call(E, p, kernel=K.arrays())
    for cfg in [(32,3,4,4),(32,3,8,8),(64,4,8,8),(128,6,8,8)]:
        for iq in (False, True):
            Harness.nq=0; K,nr,nd = b1up_kernel(H, s, *cfg, interp_q=iq); out = H.call(E, p, kernel=K.arrays())
            print(">>", gp, cfg, "interp_q", iq, "nq", Harness.nq, "max %.2e rms %.2e" % (np.max(np.abs(out/ref-1)), np.sqrt(np.mean((out/ref-1)**2))))
