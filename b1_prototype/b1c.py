"""B1 with a per-(a,i) line cache: each image angle carries a Chebyshev-Lobatto
set of traced rays between the ISCO contour and the r = 300 contour; any rin,
zone edge and node is then found by interpolation, without new traces."""
import numpy as np
from rtk import *
from contour import contours
from b1 import zone_edges
from b1up import lobatto, bary_matrix, fourier_up, tri_grid_deposit, outer_tri

def cheb_lobatto(n):
    return 0.5 * (1 - np.cos(np.pi * np.arange(n) / (n - 1)))      # [0,1]

class LineCache:
    def __init__(self, H, s, nth=64, nc=24):
        a, mu0 = s["a"], s["mu0"]; self.mus = mus = max(mu0, 0.05)
        self.th = th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
        self.Rlo, self.Rhi = s["risco"] * (1 + 1e-9), RNMAX
        rho, T = contours(H, a, mu0, mus, th, np.array([self.Rlo, self.Rhi]))
        self.u0, self.u1 = np.log(rho[0]), np.log(rho[1])
        self.x = cheb_lobatto(nc)
        u = self.u0[None, :] + self.x[:, None] * (self.u1 - self.u0)[None, :]
        rr = np.exp(u)
        r, t, st = H.trace(a, mu0, (rr * np.sin(th)).ravel(), (-rr * np.cos(th) * mus).ravel())
        T.ntrace += rr.size
        self.lr = np.log(np.maximum(r.reshape(rr.shape), 1e-300)); self.td = t.reshape(rr.shape)
        self.lr[0] = np.log(self.Rlo); self.lr[-1] = np.log(self.Rhi)
        self.ntrace = T.ntrace; self.bad = int((st != 0).sum())
        # barycentric weights for Chebyshev-Lobatto
        w = (-1.0) ** np.arange(nc); w[0] *= 0.5; w[-1] *= 0.5; self.w = w
    def interp(self, xq, Y):
        """xq: (m, nth) positions in [0,1]; Y: (nc, nth) -> (m, nth)"""
        d = xq[:, None, :] - self.x[None, :, None]            # (m, nc, nth)
        exact = np.abs(d) < 1e-15
        d = np.where(exact, 1.0, d)
        t = self.w[None, :, None] / d
        out = (t * Y[None]).sum(1) / t.sum(1)
        hit = exact.any(1)
        if hit.any():
            k = np.argmax(exact, axis=1)
            out = np.where(hit, np.take_along_axis(Y[None].repeat(xq.shape[0], 0), k[:, None, :], 1)[:, 0, :], out)
        return out
    def x_of_r(self, R):
        """Solve lr(x) = log R on every line (Newton + bisection safeguard)."""
        target = np.log(R)
        lo = np.zeros_like(self.u0); hi = np.ones_like(self.u0)
        x = np.full_like(self.u0, 0.5)
        # initial guess by linear interpolation in the node values
        for j in range(self.th.size):
            x[j] = np.interp(target, np.maximum.accumulate(self.lr[:, j]), self.x)
        for _ in range(30):
            f = self.interp(x[None], self.lr)[0] - target
            h = 1e-7
            df = (self.interp((x + h)[None], self.lr)[0] - self.interp((x - h)[None], self.lr)[0]) / (2 * h)
            lo = np.where(f < 0, x, lo); hi = np.where(f >= 0, x, hi)
            xn = x - f / df
            xn = np.where((xn <= lo) | (xn >= hi) | ~np.isfinite(xn), 0.5 * (lo + hi), xn)
            if np.max(np.abs(xn - x)) < 1e-14: x = xn; break
            x = xn
        return x

def b1c_kernel(H, s, cache, ns=4, Uth=4, Us=4, nth_out=128, nrho_out=32, interp_q=True):
    mus = cache.mus; th = cache.th; nth = th.size
    Rs = zone_edges(s); xe = s["xe"]
    xs = np.stack([cache.x_of_r(R) if R < cache.Rhi * (1 - 1e-12) else np.ones(nth) for R in Rs], 0)
    xs[0] = cache.x_of_r(Rs[0]) if Rs[0] > cache.Rlo else 0.0
    ucont = cache.u0[None, :] + xs * (cache.u1 - cache.u0)[None, :]     # contours in log rho
    lrho_f = fourier_up(ucont, Uth)
    thf = 2 * np.pi * (np.arange(nth * Uth) + 0.5) / (nth * Uth)
    sn = lobatto(ns); sf = np.linspace(0, 1, (ns - 1) * Us + 1); M = bary_matrix(sn, sf)
    K = Kernel(s); nd = 0
    for k in range(1, xe):
        u0, u1 = ucont[k - 1], ucont[k]
        uu = u0[None, :] + sn[:, None] * (u1 - u0)[None, :]
        xq = (uu - cache.u0[None, :]) / (cache.u1 - cache.u0)[None, :]
        re = np.exp(cache.interp(xq, cache.lr)); td = cache.interp(xq, cache.td)
        re[0], re[-1] = Rs[k - 1], Rs[k]
        rrn = np.exp(uu)
        qn = H.quantities((rrn * np.sin(th)).ravel(), (-rrn * np.cos(th) * mus).ravel(), re.ravel(), td.ravel(), False)
        up = lambda y: M @ fourier_up(y.reshape(uu.shape), Uth)
        Xf = up(gpos(qn["g"], s["zcos"])); lFf = up(np.log(qn["dfe"])); tauf = up(qn["tau"]); lgf = up(qn["lgsd"])
        uf0, uf1 = lrho_f[k - 1], lrho_f[k]
        rhof = np.exp(uf0[None, :] + sf[:, None] * (uf1 - uf0)[None, :])
        FJ = np.exp(lFf) * mus * rhof**2 * (uf1 - uf0)[None, :]
        nd += tri_grid_deposit(K, Xf, FJ, tauf, lgf, k, 2 * np.pi / (nth * Uth), 1.0 / ((ns - 1) * Us))
    nd += outer_tri(H, K, s, nth_out, nrho_out)
    return K, nd
