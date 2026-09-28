import sys; sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from methods import *; from contour import contours; from b1 import zone_edges, b1_gl_samples
from b1up import lobatto, fourier_up, bary_matrix
from b1c import LineCache, cheb_lobatto
from ctf import Fourier
import textwrap
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.collections import PolyCollection, LineCollection

telescope_env(); H = Harness(); E = np.logspace(-1, 2, 501)
p = DCP_Parameters(inc=60.0); H.call(E, p); s = H.state()
a, mu0, mus = s["a"], s["mu0"], s["mu0"]
Rs = zone_edges(s)
LIM = 11.0; NZ = 7            # zoom box half-width; zone edges drawn
norm = TwoSlopeNorm(vcenter=1.0, vmin=0.55, vmax=1.3); CM = "RdBu_r"
INK, MUTED = "#0b0b0b", "#8a8984"

def ab(rho, th, m=mus): return rho * np.sin(th), -rho * np.cos(th) * m
def g_of(al, be, re, td=None):
    return H.quantities(al, be, re, np.zeros_like(re) if td is None else td, False)["g"]
# true zone-edge contours (fine, for overlays)
thF = 2 * np.pi * (np.arange(512) + 0.5) / 512
rhoF, _ = contours(H, a, mu0, mus, thF, Rs[:NZ + 1])
def draw_edges(ax, lw=1.0, col=INK, ls="-"):
    for k in range(NZ + 1):
        x, y = ab(np.append(rhoF[k], rhoF[k][0]), np.append(thF, thF[0]))
        ax.plot(x, y, color=col, lw=lw if k else 1.6, ls=ls)
def frame(ax, title, sub):
    ax.set_xlim(-LIM, LIM); ax.set_ylim(-LIM * 0.8, LIM * 0.8); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_color("#d8d7d2")
    ax.set_title(title, fontsize=11.5, loc="left", fontweight="bold")
    ax.text(0.0, -0.025, textwrap.fill(" ".join(sub.split()), 78), transform=ax.transAxes, va="top", fontsize=9, color="#52514e", linespacing=1.35)

fig = plt.figure(figsize=(19, 17))
gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.55], hspace=0.62, wspace=0.16, top=0.9, bottom=0.05, left=0.04, right=0.97)

# ---- 1 pixel grid
ax = fig.add_subplot(gs[0, 0])
smp, _ = pixel_samples(H, s, 200, 200)
gr = smp[0]; m = (np.abs(gr["alpha"]) < LIM) & (np.abs(gr["beta"]) < LIM)
zb = rbin_of(gr["re"][m], s["rin"], s["xe"])
g = g_of(gr["alpha"][m], gr["beta"][m], gr["re"][m])
ax.scatter(gr["alpha"][m], gr["beta"][m], c=g, cmap=CM, norm=norm, s=2.2, lw=0, alpha=0.35 + 0.65 * (zb % 2))
draw_edges(ax, lw=0.8)
frame(ax, "1  reltrans now: pixel grid",
      "200×200 rays on a fixed polar camera (40k rays). Each pixel's whole flux goes to the one energy bin\n"
      "containing its g, and to the zone containing its r. Alternate zones are shaded: pixel membership\n"
      "is a staircase along the true zone edges (black), and the inner edge is cut by pixel centres.")
# inset zoom on staircase
ins = ax.inset_axes([0.62, 0.62, 0.36, 0.36])
ins.scatter(gr["alpha"][m], gr["beta"][m], c=g, cmap=CM, norm=norm, s=9, lw=0, alpha=0.35 + 0.65 * (zb % 2))
draw_edges(ins, lw=1.0); ins.set_xlim(2.0, 4.4); ins.set_ylim(-3.2, -1.2); ins.set_xticks([]); ins.set_yticks([])
ins.set_aspect("equal"); ax.indicate_inset_zoom(ins, edgecolor=MUTED)

# ---- 2 B1 GL + CIC
ax = fig.add_subplot(gs[0, 1])
smp, _, _ = b1_gl_samples(H, s, 32, 3, 32, 8)
b = smp[0]; m = b["rbin"] <= NZ
g = g_of(b["alpha"][m], b["beta"][m], b["re"][m], b["taudo"][m])
draw_edges(ax, lw=0.8)
ax.scatter(b["alpha"][m], b["beta"][m], c=g, cmap=CM, norm=norm, s=22, lw=0.6, edgecolors="#fcfcfb", zorder=3)
frame(ax, "2  B1, Gauss nodes + cloud-in-cell",
      "As in kerrzbb: contours solved on the zone edges (black), 3 Gauss–Legendre nodes per zone on each\n"
      "of 32 image angles. Zones are exact and the quadrature is spectral, but each node is still a point:\n"
      "its flux is split between the two nearest energy bins. Far too few points per 1.3e-3-dex bin.")

# ---- 3 Cunningham / g*
ax = fig.add_subplot(gs[0, 2])
draw_edges(ax, lw=0.5, col=MUTED)
nth = 64; th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
radii = np.unique(np.concatenate([np.exp(np.log(Rs[k - 1]) + lobatto(3) * np.log(Rs[k] / Rs[k - 1])) for k in range(1, NZ + 1)]))
rho_r, _ = contours(H, a, mu0, mus, th, radii)
Fr = Fourier(np.log(rho_r)); thf = np.linspace(0, 2 * np.pi, 2001)
chi = 0.5 * np.pi * np.linspace(0, 1, 13)
for i, rv in enumerate(radii):
    rh = np.exp(Fr(i, thf)); al, be = ab(rh, thf); gg = g_of(al, be, np.full(thf.size, rv))
    ax.plot(al, be, color="#eda100", lw=0.9)
    j0, j1 = np.argmin(gg), np.argmax(gg); g0, g1 = gg[j0], gg[j1]
    for (ta, tb, up) in ((thf[j0], thf[j1], True), (thf[j1], thf[j0] + 2 * np.pi, False)):
        if tb < ta: tb += 2 * np.pi
        tt = np.linspace(ta, tb, 800); tw = np.mod(tt, 2 * np.pi)
        rh2 = np.exp(Fr(i, tw)); al2, be2 = ab(rh2, tw); g2 = g_of(al2, be2, np.full(tt.size, rv))
        tgt = g0 + (g1 - g0) * np.sin(chi) ** 2
        key = g2 if up else -g2; tk = tgt if up else -tgt
        tn = np.interp(tk, np.maximum.accumulate(key), tt)
        rn = np.exp(Fr(i, np.mod(tn, 2 * np.pi))); x, y = ab(rn, np.mod(tn, 2 * np.pi))
        ax.scatter(x, y, c=g0 + (g1 - g0) * np.sin(chi) ** 2, cmap=CM, norm=norm, s=16, lw=0.5, edgecolors="#fcfcfb", zorder=3)
frame(ax, "3  Per-radius contours, g* (Cunningham)",
      "One contour per disc radius (orange; 3 Lobatto radii per zone), each needing a root solve per\n"
      "point. Along each contour the nodes sit at uniform χ with g* = sin²χ, so they pile up at the red-\n"
      "and blueshift extremes. Deposited with exact triangles in (r, χ). Accurate, but costs many traces.")

# ---- 4 B1 contours + exact triangles
ax = fig.add_subplot(gs[1, 0])
nth, ns, U = 24, 3, 3
th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
rho, _ = contours(H, a, mu0, mus, th, Rs[:NZ + 1])
sn = lobatto(ns); sf = np.linspace(0, 1, (ns - 1) * U + 1); M = bary_matrix(sn, sf)
lrho_f = fourier_up(np.log(rho), U); thf = 2 * np.pi * (np.arange(nth * U) + 0.5) / (nth * U)
polys, cols, coarse = [], [], []
for k in range(1, NZ + 1):
    u0, u1 = np.log(rho[k - 1]), np.log(rho[k])
    uu = u0 + sn[:, None] * (u1 - u0); rr = np.exp(uu); al, be = ab(rr, th)
    re, td, st = H.trace(a, mu0, al.ravel(), be.ravel())
    q = H.quantities(al.ravel(), be.ravel(), re, td, False)
    X = gpos(q["g"], s["zcos"]).reshape(uu.shape)
    Xf = M @ fourier_up(X, U)
    uf = lrho_f[k - 1] + sf[:, None] * (lrho_f[k] - lrho_f[k - 1]); af, bf = ab(np.exp(uf), thf)
    gf = 10 ** ((Xf - NEX // 2) * DLOGE)
    n1, n2 = Xf.shape
    for i in range(n1 - 1):
        for j in range(n2):
            j1 = (j + 1) % n2
            for tri in (((i, j), (i + 1, j), (i + 1, j1)), ((i, j), (i + 1, j1), (i, j1))):
                polys.append([(af[t], bf[t]) for t in tri]); cols.append(np.mean([gf[t] for t in tri]))
    coarse.append((al, be))
pc = PolyCollection(polys, array=np.array(cols), cmap=CM, norm=norm, edgecolors="#fcfcfb", linewidths=0.25)
ax.add_collection(pc)
for al, be in coarse:
    ax.scatter(al, be, s=9, color=INK, zorder=3, lw=0)
draw_edges(ax, lw=0.9)
frame(ax, "4  B1 contours + exact triangle deposit",
      "Same zone-edge contours; traced nodes (black dots, 24 angles × Lobatto nodes per zone) are\n"
      "upsampled spectrally (Fourier in angle, polynomial across the zone) to a fine mesh of triangles\n"
      "at no extra ray cost. Each triangle spreads its flux exactly over the energy bins it spans.")

# ---- 5 B1 line cache
ax = fig.add_subplot(gs[1, 1])
C = LineCache(H, s, 32, 24)
draw_edges(ax, lw=0.9, col="#1baf7a")
u = C.u0[None, :] + C.x[:, None] * (C.u1 - C.u0)[None, :]
al, be = ab(np.exp(u), C.th)
segs = [np.column_stack([al[:, j], be[:, j]]) for j in range(C.th.size)]
ax.add_collection(LineCollection(segs, colors="#b9b8b2", linewidths=0.7))
gl = g_of(al.ravel(), be.ravel(), np.exp(C.lr).ravel(), C.td.ravel())
inb = (np.abs(al.ravel()) < LIM) & (np.abs(be.ravel()) < LIM)
ax.scatter(al.ravel()[inb], be.ravel()[inb], c=gl[inb], cmap=CM, norm=norm, s=18, lw=0.5, edgecolors="#fcfcfb", zorder=3)
frame(ax, "5  B1 line cache + exact triangles (recommended)",
      "Traced once per (a, i): 32 image angles × 24 Chebyshev nodes from the ISCO contour to r = 300\n"
      "(most nodes lie outside this zoom). Zone edges for any rin (green) and all mesh nodes come from\n"
      "interpolation along these lines, so rin, h and frequency changes need no new rays. Then as panel 4.")

# ---- 6 kernel for one zone: what the deposit produces
ax = fig.add_subplot(gs[1, 2])
zone = 3
Kp = kernel_from_samples(H, s, pixel_samples(H, s, 200, 200)[0])
from b1up import b1up_kernel
Kr, _, _ = b1up_kernel(H, s, 128, 6, 8, 8)
from b1c import b1c_kernel
Kc, _ = b1c_kernel(H, s, LineCache(H, s, 64, 32), 4, 8, 8)
from b1 import b1_gl_samples as _g
Kg = kernel_from_samples(H, s, _g(H, s, 32, 3, 32, 8)[0], deposit_cic)
gax = 10 ** ((np.arange(NEX) + 0.5 - NEX // 2) * DLOGE)
def prof(K): return K.w0[:, 0].real.reshape(s["xe"], NEX)[zone - 1]
ref = prof(Kr); nz = np.nonzero(ref > ref.max() * 1e-4)[0]; sl = slice(nz[0] - 8, nz[-1] + 8)
sc = 1 / ref.max()
ax.step(gax[sl], prof(Kg)[sl] * sc, where="mid", color="#e87ba4", lw=0.8, label="2  B1 Gauss + CIC")
ax.step(gax[sl], prof(Kp)[sl] * sc, where="mid", color="#2a78d6", lw=0.8, label="1  pixel grid 200²")
ax.step(gax[sl], prof(Kc)[sl] * sc, where="mid", color="#1baf7a", lw=1.4, label="5  B1 line cache")
ax.plot(gax[sl], ref[sl] * sc, color=INK, lw=1.0, ls=":", label="reference")
ax.set_ylim(0, 2.2); ax.set_xlabel("g = E_obs / E_emit"); ax.set_ylabel("kernel W0, one zone (norm.)")
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.grid(True, color="#e4e3df", lw=0.6); ax.legend(frameon=False, fontsize=8.5, loc="upper left")
ax.set_title(f"6  What each method deposits: zone {zone} (r = {Rs[zone-1]:.2f}–{Rs[zone]:.2f})", fontsize=11.5, loc="left", fontweight="bold")
ax.text(0.0, -0.13, textwrap.fill("Kernel on reltrans's internal energy grid (4096 bins, Δlog g = 1.3e-3). Point deposits give shot-noise-like spikes; the triangle deposit reproduces the smooth line profile. The spikes survive the xillver convolution as the bin-to-bin noise in the spectra.", 78), transform=ax.transAxes, va="top", fontsize=9, color="#52514e", linespacing=1.35)

# ---- 7-9 deposition schematics
from matplotlib.patches import Polygon
xv = np.array([2.6, 3.4, 6.9]); W = 1.0
def tent_frac():
    x1, x2, x3 = xv
    def cdf(y):
        y = np.clip(y, x1, x3)
        return np.where(y <= x2, (y - x1) ** 2 / ((x3 - x1) * (x2 - x1)), 1 - (x3 - y) ** 2 / ((x3 - x1) * (x3 - x2)))
    k = np.arange(1, 9); return k, cdf(k) - cdf(k - 1.0)
xc = xv.mean()
schem = [("Nearest bin (pixel grid)", "Each point: whole weight to one bin.", lambda: (np.array([int(np.ceil(xc))]), np.array([1.0]))),
         ("Cloud-in-cell (B1 Gauss nodes)", "Each point: weight split linearly between\nthe two nearest bin centres.",
          lambda: (lambda x0, f: (np.array([x0 + 1, x0 + 2]), np.array([1 - f, f])))(int(np.floor(xc - 0.5)), (xc - 0.5) - np.floor(xc - 0.5))),
         ("Exact triangle (panels 3–5)", "Log g linear over the triangle: the flux\ndensity in g is a tent between the vertex\nvalues, integrated exactly over each bin.", tent_frac)]
sub = gs[2, :].subgridspec(1, 3, wspace=0.16)
for c, (tt, dsc, fn) in enumerate(schem):
    ax = fig.add_subplot(sub[0, c]); k, fr = fn()
    ax.bar(k - 0.5, fr, width=0.92, color=["#2a78d6", "#e87ba4", "#1baf7a"][c], alpha=0.85)
    for e in range(0, 10): ax.axvline(e, color="#d8d7d2", lw=0.6, zorder=0)
    if c == 2:
        x1, x2, x3 = xv; hpk = 2 / (x3 - x1)
        ax.plot([x1, x2, x3], [0, hpk, 0], color=INK, lw=1.2)
        for x in xv: ax.plot(x, 0, "o", color=INK, ms=5, clip_on=False, zorder=4)
    else:
        ax.plot(xc, 0, "o", color=INK, ms=6, clip_on=False, zorder=4)
    ax.set_xlim(0, 9); ax.set_ylim(0, 1.05); ax.set_xticks(range(0, 10)); ax.set_xticklabels([]); ax.set_xlabel("energy bins (log g)")
    ax.set_yticks([]); [ax.spines[sp].set_visible(False) for sp in ("top", "right", "left")]
    ax.set_title(tt, fontsize=10.5, loc="left", fontweight="bold")
    ax.text(0.99, 0.97, dsc, transform=ax.transAxes, ha="right", va="top", fontsize=9, color="#52514e")
    if c == 0: row3_top = ax.get_position().y1
fig.text(0.04, row3_top + 0.035, "How flux is put into energy bins  (black dots: sample positions in log g; bars: fraction of the weight in each bin)", fontsize=12, fontweight="bold", color=INK)

cax = fig.add_axes([0.72, 0.945, 0.22, 0.008]); cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CM), cax=cax, orientation="horizontal")
cb.set_label("redshift g of sample", fontsize=9.5); cax.xaxis.set_label_position("top"); cb.outline.set_visible(False)
fig.text(0.04, 0.965, "How each method samples the image plane and bins the reflection kernel", fontsize=16, fontweight="bold")
fig.text(0.04, 0.945, "a = 0.998, i = 60°, h = 6. Image plane zoomed on the inner disc (first 7 of the 20 ionisation zones, which run from rin to r = 300).", fontsize=10, color="#52514e")
fig.text(0.04, 0.93, "Image plane (α, β): colour = redshift g of each sample. Black curves: disc-radius contours at rin (thick) and at the ionisation-zone edges.",
         fontsize=9.5, color="#52514e")
fig.savefig("methods_viz.png", dpi=100)
print("ok")
