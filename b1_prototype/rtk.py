"""Harness: compute reltrans reflection kernels in Python from arbitrary
image-plane samples, and inject them back into reltrans (C hooks in subroutines/rt_hooks.f90)."""
import os, sys, time, ctypes as ct
import numpy as np
RT = os.environ.get("RELTRANS_ROOT", os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
BUILD = os.environ.get("RELTRANS_BUILD", "build")
LIBEXT = "dylib" if sys.platform == "darwin" else "so"
sys.path.insert(0, RT)
os.environ.setdefault("RELTRANS_TABLES", RT + "/cache/tables")
import pyreltrans
from pyreltrans import DCP_Parameters

dp = np.ctypeslib.ndpointer(np.float64, flags="C")
ip = np.ctypeslib.ndpointer(np.int32, flags="C")
cp = np.ctypeslib.ndpointer(np.complex64, flags="F")
I = ct.POINTER(ct.c_int); D = ct.POINTER(ct.c_double)

NEX = 4096
DLOGE = np.float32(np.log10(np.float32(3e3) / np.float32(1e-2)) / np.float32(NEX))
RNMAX = 300.0

def telescope_env():
    ci = RT + "/cache/instrument-files/"
    os.environ["RMF_SET"] = ci + "nicer-rmf6s-teamonly-array50.rmf"
    os.environ["ARF_SET"] = ci + "nicer-consim135p-teamonly-array50.arf"
    os.environ["EMIN_REF"] = "0.3"; os.environ["EMAX_REF"] = "10.0"

class Harness:
    def __init__(self, build=BUILD):
        self.r = pyreltrans.Reltrans(path=f"{RT}/{build}/lib/libreltrans.{LIBEXT}")
        L = self.lib = self.r.lib_reltrans
        L.rt_trace_batch.argtypes = [I, D, D, dp, dp, dp, dp, ip]
        L.rt_sample_quantities.argtypes = [I, I, dp, dp, dp, dp, dp, dp, dp, dp, dp]
        L.rt_flat_batch.argtypes = [I, dp, dp, dp]
        L.rt_state_info.argtypes = [I, I, I] + [D] * 10
        L.rt_get_freqs.argtypes = [I, dp]
        L.rt_get_kernel.argtypes = [I, I, I, cp, cp, cp, cp]
        L.rt_set_kernel.argtypes = [I, I, I, cp, cp]
        L.rt_timing_get.argtypes = [ct.c_void_p, ct.c_void_p]

    # ---- model calls
    def call(self, E, p, kernel=None):
        if kernel is None:
            self.lib.rt_clear_kernel()
        else:
            w0, w1 = kernel
            ne, nf, xe = w0.shape
            self.lib.rt_set_kernel(ct.byref(ct.c_int(ne)), ct.byref(ct.c_int(nf)), ct.byref(ct.c_int(xe)),
                                   np.asfortranarray(w0, np.complex64), np.asfortranarray(w1, np.complex64))
        out = self.r.dcp(E, p).astype(np.float64)
        self.lib.rt_clear_kernel()
        return out

    def state(self):
        ints = [ct.c_int() for _ in range(3)]; ds = [ct.c_double() for _ in range(10)]
        self.lib.rt_state_info(*[ct.byref(x) for x in ints + ds])
        s = dict(zip(["ne", "nf", "xe"], [x.value for x in ints]))
        s.update(zip(["a", "mu0", "rin", "rout", "h", "mudisk", "risco", "mueff", "zcos", "gamma"], [x.value for x in ds]))
        f = np.zeros(s["nf"]); self.lib.rt_get_freqs(ct.byref(ct.c_int(s["nf"])), f); s["fi"] = f
        return s

    def kernel(self):
        s = self.state(); sh = (s["ne"], s["nf"], s["xe"])
        ws = [np.zeros(sh, np.complex64, order="F") for _ in range(4)]
        self.lib.rt_get_kernel(*[ct.byref(ct.c_int(x)) for x in sh], *ws)
        return ws

    # ---- geometry
    def trace(self, spin, mu0, alpha, beta):
        alpha = np.ascontiguousarray(alpha, np.float64); beta = np.ascontiguousarray(beta, np.float64)
        n = alpha.size; r = np.empty(n); t = np.empty(n); st = np.empty(n, np.int32)
        self.lib.rt_trace_batch(ct.byref(ct.c_int(n)), ct.byref(ct.c_double(spin)), ct.byref(ct.c_double(mu0)),
                                alpha, beta, r, t, st)
        return r, t, st

    def flat(self, alpha, beta):
        alpha = np.ascontiguousarray(alpha, np.float64); beta = np.ascontiguousarray(beta, np.float64)
        n = alpha.size; r = np.empty(n)
        self.lib.rt_flat_batch(ct.byref(ct.c_int(n)), alpha, beta, r)
        return r

    nq = 0
    def quantities(self, alpha, beta, re, taudo, nonrel=False):
        Harness.nq += np.size(alpha)
        a = [np.ascontiguousarray(x, np.float64) for x in (alpha, beta, re, taudo)]
        n = a[0].size; outs = [np.empty(n) for _ in range(5)]
        self.lib.rt_sample_quantities(ct.byref(ct.c_int(n)), ct.byref(ct.c_int(int(nonrel))), *a, *outs)
        return dict(zip(["g", "dfe", "tau", "lgsd", "mue"], outs))

# ---------------------------------------------------------------- deposition
def gpos(g, zcos):
    """continuous bin coordinate: bin index k (1-based) covers x in (k-1, k]"""
    return (np.log10(g / (1.0 + zcos)) / DLOGE).astype(np.float64) + NEX // 2

def rbin_of(re, rin, xe):
    dlogr = np.log10(RNMAX / rin) / (xe - 1)
    return np.clip(np.ceil(np.log10(re / rin) / dlogr), 1, xe).astype(np.int64)

class Kernel:
    """Accumulates W0 and W1 (single lamppost, me = 1)."""
    def __init__(self, s):
        self.s = s; self.nf = s["nf"]; self.xe = s["xe"]
        self.w0 = np.zeros((NEX * self.xe, self.nf), np.complex128)
        self.w1 = np.zeros_like(self.w0)
        self.ndeposits = 0
        self.ntri = 0

    def add(self, gb, rb, w, tau, lgsd):
        """gb, rb: 1-based bins (int arrays, same length as w). w already includes
        dfe * domega * fraction."""
        ok = (gb >= 1) & (gb <= NEX)
        idx = (rb - 1) * NEX + (gb - 1)
        f = self.s["fi"]
        for j in range(self.nf):
            ph = np.exp(1j * (2 * np.pi * tau * f[j]).astype(np.float32).astype(np.float64))
            v = w * ph
            self.w0[:, j] += np.bincount(idx, weights=v.real, minlength=NEX * self.xe) + \
                1j * np.bincount(idx, weights=v.imag, minlength=NEX * self.xe)
            v1 = v * lgsd
            self.w1[:, j] += np.bincount(idx, weights=v1.real, minlength=NEX * self.xe) + \
                1j * np.bincount(idx, weights=v1.imag, minlength=NEX * self.xe)
        self.ndeposits += w.size

    def arrays(self):
        sh = (self.xe, NEX, self.nf)
        w0 = self.w0.reshape(sh).transpose(1, 2, 0); w1 = self.w1.reshape(sh).transpose(1, 2, 0)
        return np.asfortranarray(w0.astype(np.complex64)), np.asfortranarray(w1.astype(np.complex64))

def deposit_nearest(K, q, rb, domega):
    K.ntri += q["g"].size
    gb = np.clip(np.ceil(gpos(q["g"], K.s["zcos"])), 1, NEX).astype(np.int64)
    K.add(gb, rb, q["dfe"] * domega, q["tau"], q["lgsd"])

def deposit_cic(K, q, rb, domega):
    K.ntri += q["g"].size
    x = gpos(q["g"], K.s["zcos"]) - 0.5      # bin k centre at x = k - 0.5 -> coordinate k-1
    k0 = np.floor(x).astype(np.int64); fr = x - k0
    w = q["dfe"] * domega
    for kk, ww in ((k0 + 1, w * (1 - fr)), (k0 + 2, w * fr)):
        K.add(np.clip(kk, 1, NEX), rb, ww, q["tau"], q["lgsd"])

def deposit_tent(K, x3, W, rb, tau, lgsd):
    """Exact deposition of triangles with linear bin coordinate x (3 vertices,
    array (n,3)) and total weight W (already includes dfe*domega) into bins
    k covering (k-1, k].  tau, lgsd: triangle means."""
    x = np.sort(x3, axis=1)
    x1, x2, x3_ = x[:, 0], x[:, 1], x[:, 2]
    eps = 1e-9
    d31 = np.maximum(x3_ - x1, eps)
    d21 = np.maximum(x2 - x1, eps * 0.5); d32 = np.maximum(x3_ - x2, eps * 0.5)
    k_lo = np.floor(x1).astype(np.int64) + 1          # first bin touched
    k_hi = np.ceil(x3_).astype(np.int64)              # last bin touched
    k_hi = np.maximum(k_hi, k_lo)
    nb = k_hi - k_lo + 1
    tri = np.repeat(np.arange(x.shape[0]), nb)
    kk = np.repeat(k_lo, nb) + (np.arange(nb.sum()) - np.repeat(np.cumsum(nb) - nb, nb))
    def cdf(y, i):
        a1, a2, a3 = x1[i], x2[i], x3_[i]
        y = np.clip(y, a1, a3)
        lower = (y - a1) ** 2 / (d31[i] * d21[i])
        upper = 1.0 - (a3 - y) ** 2 / (d31[i] * d32[i])
        c = np.where(y <= a2, lower, upper)
        return np.clip(c, 0.0, 1.0)
    K.ntri += x.shape[0]
    frac = cdf(kk.astype(float), tri) - cdf(kk - 1.0, tri)
    # degenerate triangles (all x equal within eps): put everything in one bin
    frac = np.where(np.repeat(d31 <= 2 * eps, nb), 1.0, frac)
    m = frac > 0
    tri, kk, frac = tri[m], kk[m], frac[m]
    K.add(np.clip(kk, 1, NEX), rb[tri], W[tri] * frac, tau[tri], lgsd[tri])
    return tri.size
