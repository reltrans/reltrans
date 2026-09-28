import sys; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from methods import *
telescope_env()
H = Harness()
E = np.logspace(-1, 2, 501)
for name, p in [("dc", DCP_Parameters()), ("lag", DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)),
                ("i70", DCP_Parameters(inc=70.0, rin=-2.0))]:
    ref = H.call(E, p)
    s = H.state(); w = H.kernel()
    smp, nr = pixel_samples(H, s)
    K = kernel_from_samples(H, s, smp)
    k0, k1 = K.arrays()
    d0 = np.abs(k0 - w[0]).max() / np.abs(w[0]).max(); d1 = np.abs(k1 - w[1]).max() / np.abs(w[1]).max()
    out = H.call(E, p, kernel=(k0, k1))
    print(">>", name, "nf", s["nf"], "kernel diff W0 %.2e W1 %.2e" % (d0, d1), "spectrum diff %.2e" % np.max(np.abs(out - ref) / np.abs(ref).max()),
          "nsamples", sum(x["re"].size for x in smp), "rays", nr)
