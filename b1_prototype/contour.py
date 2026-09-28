"""Vectorised contour solver: for each image angle theta, find the image radius
rho at which the ray lands on disc radius R (primary image), in scaled polar
coordinates alpha = rho sin(th), beta = -rho cos(th) mus."""
import numpy as np

class Tracer:
    def __init__(self, H, spin, mu0, mus):
        self.H, self.spin, self.mu0, self.mus = H, spin, mu0, mus
        self.ntrace = 0
    def __call__(self, rho, th):
        al = rho * np.sin(th); be = -rho * np.cos(th) * self.mus
        r, t, st = self.H.trace(self.spin, self.mu0, al.ravel(), be.ravel())
        self.ntrace += al.size
        r = np.where(st == 0, r, 0.0).reshape(al.shape)   # misses count as r = 0
        return r, t.reshape(al.shape), st.reshape(al.shape)

def solve_contour(T, th, R, lo, hi, rlo=None, rhi=None, tol=1e-11, maxit=40):
    """Illinois on log r - log R in log rho within brackets [lo, hi] (arrays).
    r(lo) < R <= r(hi) assumed.  Returns rho, taudo at the root."""
    lR = np.log(R)
    a = np.log(lo); b = np.log(hi)
    fa = (np.log(np.maximum(rlo, 1e-300)) - lR) if rlo is not None else None
    fb = (np.log(np.maximum(rhi, 1e-300)) - lR) if rhi is not None else None
    if fa is None:
        r, _, _ = T(np.exp(a), th); fa = np.log(np.maximum(r, 1e-300)) - lR
    if fb is None:
        r, _, _ = T(np.exp(b), th); fb = np.log(np.maximum(r, 1e-300)) - lR
    fa = np.maximum(fa, -50.0)
    side = np.zeros_like(a)
    active = np.ones(a.shape, bool)
    x = 0.5 * (a + b)
    for it in range(maxit):
        # regula falsi (Illinois); bisect if the interval is huge in f
        xs = (a * fb - b * fa) / (fb - fa)
        bad = ~np.isfinite(xs) | (xs <= np.minimum(a, b)) | (xs >= np.maximum(a, b))
        xs = np.where(bad, 0.5 * (a + b), xs)
        x = np.where(active, xs, x)
        idx = np.nonzero(active)[0]
        if idx.size == 0: break
        r, _, _ = T(np.exp(x[idx]), th[idx])
        f = np.log(np.maximum(r, 1e-300)) - lR[idx]
        f = np.maximum(f, -50.0)
        fa_i, fb_i, a_i, b_i, s_i = fa[idx], fb[idx], a[idx], b[idx], side[idx]
        left = f < 0          # root is to the right of x
        # update
        a_new = np.where(left, x[idx], a_i); fa_new = np.where(left, f, fa_i)
        b_new = np.where(left, b_i, x[idx]); fb_new = np.where(left, fb_i, f)
        # Illinois modification: halve the retained end's f when the same side repeats
        fb_new = np.where(left & (s_i == -1), fb_new * 0.5, fb_new)
        fa_new = np.where(~left & (s_i == 1), fa_new * 0.5, fa_new)
        s_new = np.where(left, -1, 1)
        a[idx], fa[idx], b[idx], fb[idx], side[idx] = a_new, fa_new, b_new, fb_new, s_new
        conv = (np.abs(f) < tol) | (np.abs(b_new - a_new) < tol)
        active[idx[conv]] = False
    rho = np.exp(x)
    return rho

def contours(H, spin, mu0, mus, th, Rs, rho_min=0.3, nscan=16):
    """Solve rho_k(theta) for increasing radii Rs[k].  Returns array (K, nth)
    and the tracer (for the trace count)."""
    T = Tracer(H, spin, mu0, mus)
    nth = th.size; K = len(Rs)
    out = np.empty((K, nth))
    # inner edge: scan to bracket
    R0 = Rs[0]
    grid = np.geomspace(rho_min, (2.0 * R0 + 6.0) / mus, nscan)
    rr, _, _ = T(grid[None, :].repeat(nth, 0), th[:, None].repeat(nscan, 1))
    # last index with r < R0
    below = rr < R0
    i = nscan - 1 - np.argmax(below[:, ::-1], axis=1)       # last True
    ok = below.any(1) & (i < nscan - 1)
    if not ok.all():
        raise RuntimeError("inner-edge bracket failed for %d angles" % (~ok).sum())
    lo = grid[i]; hi = grid[i + 1]
    out[0] = solve_contour(T, th, np.full(nth, R0), lo, hi, rr[np.arange(nth), i], rr[np.arange(nth), i + 1])
    # further radii: bracket from the previous contour upward
    for k in range(1, K):
        R = Rs[k]; lo = out[k - 1].copy(); rlo = np.full(nth, Rs[k - 1])
        guess = lo * (R / Rs[k - 1])          # near-flat scaling
        hi = guess * 1.02 + 0.05
        rh, _, _ = T(hi, th)
        # expand if needed
        for _ in range(30):
            bad = rh < R
            if not bad.any(): break
            lo = np.where(bad, hi, lo); rlo = np.where(bad, rh, rlo)
            hi = np.where(bad, hi * 1.3, hi)
            rb, _, _ = T(hi[bad], th[bad]); rh[bad] = rb
        out[k] = solve_contour(T, th, np.full(nth, R), lo, hi, rlo, rh)
    return out, T
