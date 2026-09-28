"""Per-radius contour (Cunningham) method: r-contours traced by root solves,
parameterised by g* = sin^2(chi) on each branch, exact tent deposition in
(log r, chi) triangles."""
import numpy as np
from rtk import *
from contour import contours
from b1 import zone_edges
from b1up import lobatto, bary_matrix, tri_grid_deposit, outer_tri

class Fourier:
    """Trigonometric interpolant through midpoint nodes (rows = functions)."""
    def __init__(self, y):
        y = np.atleast_2d(y); self.n = n = y.shape[-1]
        th = 2 * np.pi * (np.arange(n) + 0.5) / n
        k = np.fft.fftfreq(n, 1.0 / n)
        c = np.fft.fft(y, axis=-1) / n * np.exp(-1j * k * np.pi / n)[None, :]
        if n % 2 == 0:
            c[:, n // 2] *= 0.5
            c = np.concatenate([c, c[:, n // 2:n // 2 + 1]], 1); k = np.concatenate([k, [n // 2]]); k[n // 2] = -n // 2
        self.c, self.k = c, k
    def __call__(self, row, th, der=0):
        e = np.exp(1j * np.outer(th, self.k)) * (1j * self.k) ** der
        return (e @ self.c[row]).real

def ctf_kernel(H, s, nth=64, nr=4, nchi=16, Ur=4, Uchi=2, nth_out=128, nrho_out=32):
    """nr: Lobatto radii per zone (incl. edges); nchi: chi intervals per branch."""
    a, mu0 = s["a"], s["mu0"]; mus = max(mu0, 0.05); zc = s["zcos"]
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    Rs = zone_edges(s); xe = s["xe"]
    sn = lobatto(nr)
    lR = np.log(Rs)
    radii = [np.exp(lR[k - 1] + sn * (lR[k] - lR[k - 1])) for k in range(1, xe)]
    allr = np.unique(np.concatenate(radii))
    rho, T = contours(H, a, mu0, mus, th, allr)
    rc, tc, st = H.trace(a, mu0, (rho * np.sin(th)).ravel(), (-rho * np.cos(th) * mus).ravel())
    T.ntrace += rho.size; tc = tc.reshape(rho.shape)
    FR = Fourier(np.log(rho)); FT = Fourier(tc)
    ridx = {r: i for i, r in enumerate(allr)}
    chi = 0.5 * np.pi * np.linspace(0, 1, nchi + 1)
    K = Kernel(s); nd = 0
    def gq(i, thv, rv):
        rh = np.exp(FR(i, thv)); al = rh * np.sin(thv); be = -rh * np.cos(thv) * mus
        q = H.quantities(al, be, np.full(thv.size, rv), FT(i, thv), False)
        return q, rh
    thf = 2 * np.pi * (np.arange(nth * 8) + 0.5) / (nth * 8)
    for k in range(1, xe):
        rk = radii[k - 1]; ii = [ridx[r] for r in rk]
        # per radius and branch: theta(chi), and the jacobian ingredients
        TH = np.empty((2, nr, nchi + 1)); DTH = np.empty_like(TH)
        for m, (i, rv) in enumerate(zip(ii, rk)):
            q, _ = gq(i, thf, rv); g = q["g"]
            jmin, jmax = np.argmin(g), np.argmax(g)
            def refine(j0, sgn):
                t = thf[j0]
                for _ in range(6):   # Newton on dg/dth using finite differences of the smooth interpolant
                    h = 1e-4
                    gm = gq(i, np.array([t - h, t, t + h]), rv)[0]["g"]
                    d1 = (gm[2] - gm[0]) / (2 * h); d2 = (gm[2] - 2 * gm[1] + gm[0]) / h**2
                    t = t - d1 / d2
                gm = gq(i, np.array([t - 1e-4, t, t + 1e-4]), rv)[0]["g"]
                return t, gm[1], (gm[2] - 2 * gm[1] + gm[0]) / 1e-8
            t0, g0, d20 = refine(jmin, 1); t1, g1, d21 = refine(jmax, -1)
            if t1 < t0: t1 += 2 * np.pi
            D = g1 - g0
            for b in range(2):
                # branch 0: t0 -> t1 (g rising); branch 1: t1 -> t0 + 2pi (g falling)
                ta, tb = (t0, t1) if b == 0 else (t1, t0 + 2 * np.pi)
                target = g0 + D * np.sin(chi) ** 2 if b == 0 else g1 - D * np.sin(chi) ** 2
                # initial guess by monotone interpolation on a fine sampling of the branch
                tt = np.linspace(ta, tb, 400); gg = gq(i, tt, rv)[0]["g"]
                if b == 1: gg = -gg; tgt = -target
                else: tgt = target
                gg = np.maximum.accumulate(gg)
                t = np.interp(tgt, gg, tt)
                t[0], t[-1] = ta, tb
                for _ in range(4):
                    inner = slice(1, nchi)
                    h = 1e-6
                    gp = gq(i, t[inner] + h, rv)[0]["g"]; gmn = gq(i, t[inner] - h, rv)[0]["g"]; gc = gq(i, t[inner], rv)[0]["g"]
                    d1 = (gp - gmn) / (2 * h)
                    t[inner] = np.clip(t[inner] - (gc - target[inner]) / d1, ta, tb)
                TH[b, m] = t
                # dtheta/dchi = D sin(2chi) / |g'|; endpoint limits sqrt(2D/|g''|)
                gp = gq(i, t + 1e-6, rv)[0]["g"]; gmn = gq(i, t - 1e-6, rv)[0]["g"]
                d1 = np.abs(gp - gmn) / 2e-6
                with np.errstate(divide="ignore", invalid="ignore"):
                    dt = D * np.sin(2 * chi) / d1
                dt[0] = np.sqrt(2 * D / abs(d20 if b == 0 else d21)) if True else 0
                dt[-1] = np.sqrt(2 * D / abs(d21 if b == 0 else d20))
                DTH[b, m] = dt
        # upsample in log r (polynomial, Lobatto nodes) and chi (polynomial per sub-interval is overkill: linear refine)
        sf = np.linspace(0, 1, (nr - 1) * Ur + 1); M = bary_matrix(sn, sf)
        Md = None
        lr_f = np.log(rk[0]) + sf * (np.log(rk[-1]) - np.log(rk[0]))
        chif = 0.5 * np.pi * np.linspace(0, 1, nchi * Uchi + 1)
        if k == 1: ctf_kernel.dbg = (TH.copy(), DTH.copy(), chi, rk)
        # make theta(chi) continuous across the radii of the zone (2 pi ambiguity)
        for b in range(2):
            for m in range(1, nr):
                TH[b, m] += 2 * np.pi * np.round((TH[b, m - 1, nchi // 2] - TH[b, m, nchi // 2]) / (2 * np.pi))
        for b in range(2):
            THf = M @ TH[b]; DTHf = M @ DTH[b]
            if Uchi > 1:
                THf = np.array([np.interp(chif, chi, row) for row in THf])   # (cheap; refined below via cubic?)
                DTHf = np.array([np.interp(chif, chi, row) for row in DTHf])
            # rho and d rho/d r at (theta, r_fine): interpolate log rho across the zone radii at fixed theta
            nrf, nc = THf.shape
            lrho_nodes = np.stack([FR(i, THf.ravel()) for i in ii], 0)          # (nr, nrf*nc)
            td_nodes = np.stack([FT(i, THf.ravel()) for i in ii], 0)
            # weights for value and derivative in u = log r
            u_nodes = np.log(rk); ufine = np.repeat(lr_f, nc)
            Wv = bary_matrix(u_nodes, lr_f)                                         # (nrf, nr)
            # derivative matrix by differentiating Lagrange basis numerically
            hstep = 1e-6 * (u_nodes[-1] - u_nodes[0])
            Wd = (bary_matrix(u_nodes, lr_f + hstep) - bary_matrix(u_nodes, lr_f - hstep)) / (2 * hstep)
            idx_r = np.repeat(np.arange(nrf), nc)
            lrho = np.einsum("pn,np->p", Wv[idx_r], lrho_nodes)
            dlrho = np.einsum("pn,np->p", Wd[idx_r], lrho_nodes)                   # d log rho / d log r
            td = np.einsum("pn,np->p", Wv[idx_r], td_nodes)
            rho_ = np.exp(lrho); r_ = np.exp(ufine); thv = THf.ravel()
            al = rho_ * np.sin(thv); be = -rho_ * np.cos(thv) * mus
            q = H.quantities(al, be, r_, td, False)
            sh = (nrf, nc)
            X = gpos(q["g"], zc).reshape(sh)
            if ctf_kernel.area_mode:
                nd += tri_area_deposit(K, X, q["dfe"].reshape(sh), al.reshape(sh), be.reshape(sh),
                                       q["tau"].reshape(sh), q["lgsd"].reshape(sh), k)
            else:
                J = mus * rho_**2 * np.abs(dlrho) * DTHf.ravel()
                FJ = (q["dfe"] * J).reshape(sh) * (u_nodes[-1] - u_nodes[0])
                nd += tri_open_deposit(K, X, FJ, q["tau"].reshape(sh), q["lgsd"].reshape(sh), k,
                                       1.0 / (nrf - 1), chif[1] - chif[0])
    nd += outer_tri(H, K, s, nth_out, nrho_out)
    return K, T.ntrace, nd

ctf_kernel.area_mode = False

def tri_open_deposit(K, X, FJ, tau, lg, rb, ds, dc):
    n1, n2 = X.shape
    i0 = np.arange(n1 - 1)[:, None]; i1 = i0 + 1; j0 = np.arange(n2 - 1)[None, :]; j1 = j0 + 1
    A = 0.5 * ds * dc; nd = 0
    for (pa, pb, pc) in (((i0, j0), (i1, j0), (i1, j1)), ((i0, j0), (i1, j1), (i0, j1))):
        idx = [np.broadcast_arrays(*p) for p in (pa, pb, pc)]
        xs = np.stack([X[a_, b_].ravel() for a_, b_ in idx], 1)
        W = A * np.mean([FJ[a_, b_].ravel() for a_, b_ in idx], 0)
        tm = np.mean([tau[a_, b_].ravel() for a_, b_ in idx], 0)
        lm = np.mean([lg[a_, b_].ravel() for a_, b_ in idx], 0)
        nd += deposit_tent(K, xs, W, np.full(W.size, rb, np.int64), tm, lm)
    return nd

def tri_area_deposit(K, X, F, AL, BE, tau, lg, rb):
    """Triangles on an open logical grid; weight = image-plane area x mean F."""
    n1, n2 = X.shape
    i0 = np.arange(n1 - 1)[:, None]; i1 = i0 + 1; j0 = np.arange(n2 - 1)[None, :]; j1 = j0 + 1
    nd = 0
    for (pa, pb, pc) in (((i0, j0), (i1, j0), (i1, j1)), ((i0, j0), (i1, j1), (i0, j1))):
        idx = [tuple(np.broadcast_arrays(*p)) for p in (pa, pb, pc)]
        xs = np.stack([X[t].ravel() for t in idx], 1)
        ax = [AL[t].ravel() for t in idx]; bx = [BE[t].ravel() for t in idx]
        area = 0.5 * np.abs((ax[1] - ax[0]) * (bx[2] - bx[0]) - (ax[2] - ax[0]) * (bx[1] - bx[0]))
        W = area * np.mean([F[t].ravel() for t in idx], 0)
        tm = np.mean([tau[t].ravel() for t in idx], 0)
        lm = np.mean([lg[t].ravel() for t in idx], 0)
        nd += deposit_tent(K, xs, W, np.full(W.size, rb, np.int64), tm, lm)
    return nd
