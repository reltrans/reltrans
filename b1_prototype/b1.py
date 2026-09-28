"""B1 samplers for reltrans kernels."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from rtk import *
from contour import contours, Tracer

def zone_edges(s):
    xe = s["xe"]; dl = np.log10(RNMAX / s["rin"]) / (xe - 1)
    return s["rin"] * 10 ** (dl * np.arange(xe))          # R_0 = rin ... R_{xe-1} = 300

def outer_samples(H, s, nth, nrho):
    """Zone xe: flat-space disc from 300 to rout, ellipse coords (exact r for thin disc)."""
    mu0 = s["mu0"]
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    x, w = leggauss(nrho)
    u0, u1 = np.log(RNMAX), np.log(s["rout"])
    u = 0.5 * (u1 - u0) * x + 0.5 * (u1 + u0); wu = 0.5 * (u1 - u0) * w
    rho = np.exp(u)
    al = rho[:, None] * np.sin(th)[None, :]; be = -rho[:, None] * np.cos(th)[None, :] * mu0
    wt = (mu0 * rho**2 * wu)[:, None] * np.full(nth, 2 * np.pi / nth)[None, :]
    re = H.flat(al.ravel(), be.ravel())
    return dict(alpha=al.ravel(), beta=be.ravel(), re=re, taudo=np.zeros(re.size),
                w=wt.ravel(), nonrel=True, rbin=np.full(re.size, s["xe"], np.int64))

def b1_gl_samples(H, s, nth=64, ngl=3, nth_out=64, nrho_out=8):
    a, mu0 = s["a"], s["mu0"]; mus = max(mu0, 0.05)
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    Rs = zone_edges(s)
    rho, T = contours(H, a, mu0, mus, th, Rs)
    x, w = leggauss(ngl)
    al, be, wt, rb = [], [], [], []
    for k in range(1, len(Rs)):
        u0 = np.log(rho[k - 1]); u1 = np.log(rho[k])        # (nth,)
        u = 0.5 * (u1 - u0)[None, :] * x[:, None] + 0.5 * (u1 + u0)[None, :]
        wu = 0.5 * (u1 - u0)[None, :] * w[:, None]
        rr = np.exp(u)
        al.append((rr * np.sin(th)).ravel()); be.append((-rr * np.cos(th) * mus).ravel())
        wt.append((mus * rr**2 * wu * 2 * np.pi / nth).ravel()); rb.append(np.full(rr.size, k, np.int64))
    al = np.concatenate(al); be = np.concatenate(be)
    r, t, st = H.trace(a, mu0, al, be); T.ntrace += al.size
    inner = dict(alpha=al, beta=be, re=r, taudo=t, w=np.concatenate(wt), nonrel=False, rbin=np.concatenate(rb))
    bad = (st != 0).sum()
    return [inner, outer_samples(H, s, nth_out, nrho_out)], T.ntrace, bad

def b1_tri_kernel(H, s, nth=64, nsub=4, nth_out=64, nrho_out=8, upth=1):
    """B1 vertex grid (contours at zone edges, nsub log-rho intervals per zone),
    triangles in (theta, s), exact tent deposition in log g."""
    a, mu0 = s["a"], s["mu0"]; mus = max(mu0, 0.05)
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    Rs = zone_edges(s); xe = s["xe"]
    rho, T = contours(H, a, mu0, mus, th, Rs)
    K = Kernel(s); ndep = 0
    # vertex rho for every zone: (nsub+1, nth), traced (edges reuse the contour value r = R)
    for k in range(1, xe):
        u0 = np.log(rho[k - 1]); u1 = np.log(rho[k])
        sv = np.linspace(0, 1, nsub + 1)
        u = u0[None, :] + sv[:, None] * (u1 - u0)[None, :]
        rr = np.exp(u)
        alv = rr * np.sin(th); bev = -rr * np.cos(th) * mus
        # trace interior vertices only; edge vertices have r = R exactly (taudo needs a trace)
        r, t, st = H.trace(a, mu0, alv.ravel(), bev.ravel()); T.ntrace += alv.size
        r = r.reshape(rr.shape); t = t.reshape(rr.shape)
        q = H.quantities(alv.ravel(), bev.ravel(), r.ravel(), t.ravel(), False)
        g = q["g"].reshape(rr.shape); F = q["dfe"].reshape(rr.shape)
        tau = q["tau"].reshape(rr.shape); lg = q["lgsd"].reshape(rr.shape)
        J = mus * rr**2 * (u1 - u0)[None, :]          # d Omega = J ds dth
        X = gpos(g, s["zcos"])
        FJ = F * J
        # quads (i, j) -> triangles (i,j),(i+1,j),(i+1,j+1) and (i,j),(i+1,j+1),(i,j+1)
        j0 = np.arange(nth); j1 = (j0 + 1) % nth
        i0 = np.arange(nsub)[:, None]; i1 = i0 + 1
        A = 0.5 * (1.0 / nsub) * (2 * np.pi / nth)     # logical triangle area
        for (pa, pb, pc) in (((i0, j0), (i1, j0), (i1, j1)), ((i0, j0), (i1, j1), (i0, j1))):
            idx = [np.broadcast_arrays(p[0], p[1][None, :]) for p in (pa, pb, pc)]
            xs = np.stack([X[ii, jj].ravel() for ii, jj in idx], 1)
            W = A * np.mean([FJ[ii, jj].ravel() for ii, jj in idx], 0)
            tm = np.mean([tau[ii, jj].ravel() for ii, jj in idx], 0)
            lm = np.mean([lg[ii, jj].ravel() for ii, jj in idx], 0)
            ndep += deposit_tent(K, xs, W, np.full(W.size, k, np.int64), tm, lm)
    # outer zone: GL + CIC (smooth, few bins)
    osm = outer_samples(H, s, nth_out, nrho_out)
    q = H.quantities(osm["alpha"], osm["beta"], osm["re"], osm["taudo"], True)
    deposit_cic(K, q, osm["rbin"], osm["w"])
    return K, T.ntrace, ndep
