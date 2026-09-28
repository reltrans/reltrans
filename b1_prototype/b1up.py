"""B1 with Lobatto nodes per zone, spectral upsampling (Fourier in theta,
polynomial in s) and exact tent deposition on triangles."""
import numpy as np
from rtk import *
from contour import contours, Tracer
from b1 import zone_edges

def lobatto(n):
    """n Gauss-Lobatto-Legendre nodes on [0,1] (n >= 2)."""
    if n == 2: return np.array([0.0, 1.0])
    from numpy.polynomial import legendre as L
    c = np.zeros(n); c[-1] = 1            # P_{n-1}
    x = np.sort(L.legroots(L.legder(c)))
    return 0.5 * (np.concatenate([[-1.0], x, [1.0]]) + 1)

def bary_matrix(xn, xe):
    """Barycentric Lagrange interpolation matrix from nodes xn to points xe."""
    w = np.array([1.0 / np.prod(xn[j] - np.delete(xn, j)) for j in range(len(xn))])
    M = np.zeros((len(xe), len(xn)))
    for i, x in enumerate(xe):
        d = x - xn
        k = np.nonzero(np.abs(d) < 1e-14)[0]
        if k.size: M[i, k[0]] = 1.0; continue
        t = w / d; M[i] = t / t.sum()
    return M

_FM = {}
def fourier_matrix(n, U):
    """Trigonometric interpolation matrix from n midpoint nodes to n*U midpoint nodes."""
    key = (n, U)
    if key in _FM: return _FM[key]
    m = n * U
    th = 2 * np.pi * (np.arange(n) + 0.5) / n
    thf = 2 * np.pi * (np.arange(m) + 0.5) / m
    k = np.fft.fftfreq(n, 1.0 / n)                 # integer wavenumbers
    wk = np.ones(n)
    if n % 2 == 0:
        k[n // 2] = n // 2; wk[n // 2] = 0.5       # Nyquist split into +-n/2 (cos part)
    # value at thf = (1/n) sum_j y_j sum_k w_k cos(k (thf - th_j))  (real, symmetric k)
    d = thf[:, None] - th[None, :]
    Mx = np.zeros((m, n))
    for kk, w in zip(k, wk):
        if kk < 0: continue
        fac = 1.0 if kk == 0 else 2.0
        if n % 2 == 0 and kk == n // 2: fac = 1.0
        Mx += fac * np.cos(kk * d)
    Mx /= n
    _FM[key] = Mx
    return Mx

def fourier_up(y, U):
    if U == 1: return y
    n = y.shape[-1]
    return y @ fourier_matrix(n, U).T

def tri_grid_deposit(K, X, FJ, tau, lg, rb, dth, ds):
    """Vertex grid (ns+1, nth) periodic in theta; triangles in logical coords."""
    ns1, nth = X.shape
    j0 = np.arange(nth); j1 = (j0 + 1) % nth
    i0 = np.arange(ns1 - 1)[:, None]; i1 = i0 + 1
    A = 0.5 * dth * ds
    nd = 0
    for (pa, pb, pc) in (((i0, j0), (i1, j0), (i1, j1)), ((i0, j0), (i1, j1), (i0, j1))):
        idx = [np.broadcast_arrays(p[0], p[1][None, :]) for p in (pa, pb, pc)]
        xs = np.stack([X[ii, jj].ravel() for ii, jj in idx], 1)
        W = A * np.mean([FJ[ii, jj].ravel() for ii, jj in idx], 0)
        tm = np.mean([tau[ii, jj].ravel() for ii, jj in idx], 0)
        lm = np.mean([lg[ii, jj].ravel() for ii, jj in idx], 0)
        nd += deposit_tent(K, xs, W, np.full(W.size, rb, np.int64), tm, lm)
    return nd

def outer_tri(H, K, s, nth=128, nrho=32):
    mu0 = s["mu0"]
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    u = np.linspace(np.log(RNMAX), np.log(s["rout"]), nrho + 1)
    rho = np.exp(u)[:, None] * np.ones(nth)[None, :]
    al = rho * np.sin(th); be = -rho * np.cos(th) * mu0
    re = H.flat(al.ravel(), be.ravel())
    q = H.quantities(al.ravel(), be.ravel(), re, np.zeros(re.size), True)
    sh = rho.shape
    X = gpos(q["g"], s["zcos"]).reshape(sh)
    FJ = (q["dfe"].reshape(sh)) * mu0 * rho**2 * (u[-1] - u[0])
    return tri_grid_deposit(K, X, FJ, q["tau"].reshape(sh), q["lgsd"].reshape(sh), s["xe"], 2 * np.pi / nth, 1.0 / nrho)

def b1up_kernel(H, s, nth=64, ns=4, Uth=4, Us=4, nth_out=128, nrho_out=32, retrace_edges=True, interp_q=False):
    """ns: Lobatto nodes per zone segment (incl. both contour ends)."""
    a, mu0 = s["a"], s["mu0"]; mus = max(mu0, 0.05)
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    Rs = zone_edges(s); xe = s["xe"]
    rho, T = contours(H, a, mu0, mus, th, Rs)
    # trace the contour points once more for taudo (and exact r)
    rc, tc, st = H.trace(a, mu0, (rho * np.sin(th)).ravel(), (-rho * np.cos(th) * mus).ravel())
    T.ntrace += rho.size
    tc = tc.reshape(rho.shape)
    lrho_f = fourier_up(np.log(rho), Uth)           # (xe, nth*Uth)
    thf = 2 * np.pi * (np.arange(nth * Uth) + 0.5) / (nth * Uth)
    sn = lobatto(ns); sf = np.linspace(0, 1, (ns - 1) * Us + 1)
    M = bary_matrix(sn, sf)
    K = Kernel(s); nd = 0
    for k in range(1, xe):
        u0, u1 = np.log(rho[k - 1]), np.log(rho[k])
        # nodes: rows = s nodes
        uu = u0[None, :] + sn[:, None] * (u1 - u0)[None, :]
        re = np.empty(uu.shape); td = np.empty(uu.shape)
        re[0], re[-1] = Rs[k - 1], Rs[k]; td[0], td[-1] = tc[k - 1], tc[k]
        if ns > 2:
            rr = np.exp(uu[1:-1])
            r_, t_, st_ = H.trace(a, mu0, (rr * np.sin(th)).ravel(), (-rr * np.cos(th) * mus).ravel())
            T.ntrace += rr.size
            re[1:-1] = r_.reshape(rr.shape); td[1:-1] = t_.reshape(rr.shape)
        if interp_q:
            rrn = np.exp(uu)
            qn = H.quantities((rrn * np.sin(th)).ravel(), (-rrn * np.cos(th) * mus).ravel(), re.ravel(), td.ravel(), False)
            up = lambda y: M @ fourier_up(y.reshape(uu.shape), Uth)
            Xf = up(gpos(qn["g"], s["zcos"])); lFf = up(np.log(qn["dfe"])); tauf = up(qn["tau"]); lgf = up(qn["lgsd"])
            uf0, uf1 = lrho_f[k - 1], lrho_f[k]
            ufine = uf0[None, :] + sf[:, None] * (uf1 - uf0)[None, :]
            rhof = np.exp(ufine)
            FJ = np.exp(lFf) * mus * rhof**2 * (uf1 - uf0)[None, :]
            nd += tri_grid_deposit(K, Xf, FJ, tauf, lgf, k, 2 * np.pi / (nth * Uth), 1.0 / ((ns - 1) * Us))
            continue
        # upsample: theta first (Fourier), then s (polynomial)
        lre = fourier_up(np.log(re), Uth); tdf = fourier_up(td, Uth)
        lre = M @ lre; tdf = M @ tdf
        # exact edge radii on the fine grid
        lre[0] = np.log(Rs[k - 1]); lre[-1] = np.log(Rs[k])
        uf0, uf1 = lrho_f[k - 1], lrho_f[k]
        ufine = uf0[None, :] + sf[:, None] * (uf1 - uf0)[None, :]
        rhof = np.exp(ufine)
        al = rhof * np.sin(thf); be = -rhof * np.cos(thf) * mus
        q = H.quantities(al.ravel(), be.ravel(), np.exp(lre).ravel(), tdf.ravel(), False)
        sh = rhof.shape
        X = gpos(q["g"], s["zcos"]).reshape(sh)
        FJ = q["dfe"].reshape(sh) * mus * rhof**2 * (uf1 - uf0)[None, :]
        nd += tri_grid_deposit(K, X, FJ, q["tau"].reshape(sh), q["lgsd"].reshape(sh), k,
                               2 * np.pi / (nth * Uth), 1.0 / ((ns - 1) * Us))
    nd += outer_tri(H, K, s, nth_out, nrho_out)
    return K, T.ntrace, nd
