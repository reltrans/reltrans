import sys, collections
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from analyse import load_all, T_Q
T_INTERP = 20e-9
rows = load_all()
# aggregate: worst over geometries, per (method, group)
GROUPS = {"Time-averaged spectrum (total)": ("dc", ["dc"]), "Time-averaged, reflection only": ("dc", ["refl"]),
          "Cross spectrum, worst of Re/Im/lag": (("lo", "hi"), ["re", "im", "lag"])}
FAM = [("pix", "reltrans pixel grid (nearest bin)", "#2a78d6", "o"),
       ("b1gl", "B1, Gauss nodes + cloud-in-cell", "#e87ba4", "v"),
       ("ctf", "per-radius contours, g* (Cunningham)", "#eda100", "s"),
       ("b1up", "B1 contours + exact triangle deposit", "#eb6834", "D"),
       ("b1c", "B1 line cache + exact triangle deposit", "#1baf7a", "^")]
def tcost(r, cached):
    t = r["t"]
    if r["method"].startswith("b1c"): t += T_INTERP * 0  # added below via ntri
    return t - (r["t_ray"] if cached else 0)
d = collections.defaultdict(lambda: collections.defaultdict(list))
for r in rows:
    if r["method"] in ("reltrans",) or r["method"].startswith("b1tri") or r["method"] == "b1up_256x8_U8x8": continue
    for gname, (fs, fls) in GROUPS.items():
        fs = fs if isinstance(fs, tuple) else (fs,)
        if r["fset"] in fs and r["fl"] in fls:
            d[gname][r["method"]].append(r)
fig, axes = plt.subplots(3, 2, figsize=(12, 13), sharey="row")
for gi, gname in enumerate(GROUPS):
    for ci, cached in enumerate((False, True)):
        ax = axes[gi, ci]
        for fam, lab, col, mk in FAM:
            if cached and fam in ("b1up", "ctf", "b1gl"): continue
            pts = []
            for m, rs in d[gname].items():
                if not m.startswith(fam): continue
                if len({(r["geom"]) for r in rs}) < 5: continue
                e = max(r["err"] for r in rs)
                # cost: mean over geometries and frequency sets (dc: nf=1; cross: mean of lo/hi)
                t = np.mean([r["t"] - (r["t_ray"] if cached else 0) for r in rs])
                pts.append((t * 1e3, e, m))
            if not pts: continue
            pts.sort(); t_, e_, m_ = zip(*pts)
            ax.plot(t_, e_, "-", color=col, lw=2, marker=mk, ms=8, mec="#fcfcfb", mew=1.5, label=lab)
            if fam == "pix":
                for t0, e0, m0 in pts:
                    if m0 in ("pix_200", "pix_800", "pix_1600"):
                        ax.annotate(m0.split("_")[1] + "²" + (" (default)" if m0 == "pix_200" else ""), (t0, e0), textcoords="offset points", xytext=(6, 4), fontsize=8, color="#52514e")
        ax.axhline(2e-4, color="#52514e", lw=1, ls="--")
        ax.text(0.99, 2e-4, "test rtol 2e-4", transform=ax.get_yaxis_transform(), ha="right", va="bottom", fontsize=8, color="#52514e")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.grid(True, which="major", color="#e4e3df", lw=0.6)
        ax.set_title(gname + ("\ncost with geometry cached" if cached else "\ncost incl. tracing (new a, i)"), fontsize=10)
        if ci == 0: ax.set_ylabel("max error vs reference\n(worst of 5 geometries)")
        if gi == 2: ax.set_xlabel("modelled kernel cost per call [ms]")
        for s in ("top", "right"): ax.spines[s].set_visible(False)
h, l = axes[0, 0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=3, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.965))
for ax in axes[2]: ax.annotate("pixel floor: lag error from the\nr ≈ 298–300 gap at the GR/flat seam", xy=(0.97, 0.80), xycoords="axes fraction", ha="right", fontsize=8, color="#52514e")
fig.suptitle("reltrans reflection kernels: error against cost", fontsize=13, y=0.995)
fig.text(0.01, 0.005, "Cost model (Fortran, measured): 3.8 µs/ray, 0.41 µs per per-sample quantity eval, per frequency 15 ns/triangle + 9 ns/bin deposit. "
         "Cross-spectrum costs are for nf = 3 and 6.\nReference: B1 contours, 256×8 Lobatto, 8× upsampling. Errors: spectra relative, Im and lag relative to max |value|.", fontsize=7.5, color="#52514e")
fig.tight_layout(rect=(0, 0.03, 1, 0.94)); fig.savefig("kernel_cost.png", dpi=110)
