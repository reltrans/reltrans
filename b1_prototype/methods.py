import numpy as np
from rtk import *

def rfunc(a, mu0):
    if a > 0.8:
        r = 1.5 + 0.5 * mu0**5.5; r = min(r, -0.1 + 5.6 * mu0)
    else:
        r = 3.0 + 0.5 * mu0**5.5; r = min(r, -0.2 + 10.0 * mu0)
    return max(0.1, r)

def getrgrid(rnmin, rnmax, mueff, nro, nphi):
    i = np.arange(1, nro + 1)
    dlogr = np.log10(rnmax / rnmin) / nro
    rar = np.concatenate([[rnmin], 10.0 ** (np.log10(rnmin) + i * dlogr)])
    rn = 0.5 * (rar[1:] + rar[:-1])
    dom = rn * (rar[1:] - rar[:-1]) * mueff * 2 * np.pi / nphi
    return rn, dom

def _early_exit_mask(hit):
    """hit: (nro, nphi) bool of in-disc samples. reltrans loops rings from
    outermost inwards and stops at the first empty ring after the disc was seen."""
    keep = np.zeros(hit.shape[0], bool); seen = False
    for i in range(hit.shape[0] - 1, -1, -1):
        any_hit = hit[i].any()
        keep[i] = True
        if seen and not any_hit:
            keep[i] = False; keep[:i] = False; break
        if any_hit: seen = True
    return keep

def pixel_samples(H, s, nro=200, nphi=200, nron=100, nphin=100):
    """The current reltrans camera. Returns list of sample dicts."""
    a, mu0, mueff = s["a"], s["mu0"], s["mueff"]
    out = []; nrays = 0
    for nonrel, (r0, r1, n1, n2) in ((False, (rfunc(a, mu0), RNMAX, nro, nphi)),
                                     (True, (RNMAX, s["rout"], nron, nphin))):
        rn, dom = getrgrid(r0, r1, mueff, n1, n2)
        phin = (np.arange(1, n2 + 1) - 0.5) * 2 * np.pi / n2
        al = rn[:, None] * np.sin(phin)[None, :]; be = -rn[:, None] * np.cos(phin)[None, :] * mueff
        if nonrel:
            re = H.flat(al.ravel(), be.ravel()).reshape(al.shape); td = np.zeros_like(re)
            ok = np.ones(al.shape, bool)
        else:
            re, td, st = H.trace(a, mu0, al.ravel(), be.ravel()); nrays += al.size
            re = re.reshape(al.shape); td = td.reshape(al.shape)
            ok = (st.reshape(al.shape) == 0) & (re > s["risco"]) & (re < s["rout"])
        hit = ok & (re >= s["rin"]) & (re <= s["rout"])
        hit &= _early_exit_mask(hit)[:, None]
        w = np.broadcast_to(dom[:, None], al.shape)
        out.append(dict(alpha=al[hit], beta=be[hit], re=re[hit], taudo=td[hit], w=w[hit], nonrel=nonrel))
    return out, nrays

def kernel_from_samples(H, s, samples, deposit=deposit_nearest):
    K = Kernel(s)
    for sm in samples:
        if sm["re"].size == 0: continue
        q = H.quantities(sm["alpha"], sm["beta"], sm["re"], sm["taudo"], sm["nonrel"])
        rb = sm.get("rbin")
        if rb is None: rb = rbin_of(sm["re"], s["rin"], s["xe"])
        deposit(K, q, rb, sm["w"])
    return K
