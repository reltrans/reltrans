import sys; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from methods import *
telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
d = dict(np.load("results/default_lo.npz"))
p = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0); H.call(E, p); s = H.state()
def pixel_nogap(H, s, n):
    a, mu0, mueff = s["a"], s["mu0"], s["mueff"]
    rn, dom = getrgrid(rfunc(a, mu0), 330.0, mueff, int(n*1.02), n)
    phin = (np.arange(1, n + 1) - 0.5) * 2 * np.pi / n
    al = rn[:, None] * np.sin(phin); be = -rn[:, None] * np.cos(phin) * mueff
    re, td, st = H.trace(a, mu0, al.ravel(), be.ravel()); re = re.reshape(al.shape); td = td.reshape(al.shape)
    hit = (st.reshape(al.shape) == 0) & (re > s["risco"]) & (re >= s["rin"]) & (re <= RNMAX)
    w = np.broadcast_to(dom[:, None], al.shape)
    gr = dict(alpha=al[hit], beta=be[hit], re=re[hit], taudo=td[hit], w=w[hit], nonrel=False)
    smp, _ = pixel_samples(H, s, 10, 10, n // 2, n // 2)
    return [gr, smp[1]]
for n in (400, 800):
    K = kernel_from_samples(H, s, pixel_nogap(H, s, n)); out = H.call(E, p, kernel=K.arrays())
    r = d["b1up_256x8_U8x8/lag"]
    print(">> nogap pix", n, "lag err %.2e" % (np.max(np.abs(out - r)) / np.max(np.abs(r))), " (with gap: %.2e)" % (np.max(np.abs(d[f'pix_{n}/lag'] - r)) / np.max(np.abs(r))))
