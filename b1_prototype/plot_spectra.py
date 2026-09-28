import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
E = np.logspace(-1, 2, 501); Ec = np.sqrt(E[1:] * E[:-1]); dE = np.diff(E)
REF, B1, PIX = "b1up_256x8_U8x8", "b1c_64x32_n4_U8x8", "reltrans"
C_PIX, C_B1, C_REF, INK = "#2a78d6", "#1baf7a", "#0b0b0b", "#52514e"
GEO = [("default", "a = 0.998, i = 30°, h = 6"), ("i70", "a = 0.998, i = 70°, h = 6")]
def style(ax):
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.grid(True, color="#e4e3df", lw=0.6); ax.set_xscale("log")

def fig_timeavg(fl, title, fname):
    fig, axes = plt.subplots(3, 2, figsize=(13, 10), gridspec_kw=dict(height_ratios=[2.2, 1, 1]))
    for c, (g, lab) in enumerate(GEO):
        d = np.load(f"results/{g}_dc.npz")
        r, p, b = d[f"{REF}/{fl}"], d[f"{PIX}/{fl}"], d[f"{B1}/{fl}"]
        f = lambda x: x / dE * Ec**2
        ax = axes[0, c]
        ax.plot(Ec, f(r), color=C_REF, lw=2.5, label="converged reference")
        ax.plot(Ec, f(p), color=C_PIX, lw=1, label="reltrans now (200² pixels, 40k rays)")
        ax.plot(Ec, f(b), color=C_B1, lw=1, ls="--", label="B1 line cache 64×32 (3.7k rays)")
        ax.set_yscale("log"); ax.set_ylabel("E² F(E)  [arb.]"); ax.set_title(lab); style(ax)
        ax.set_ylim(f(r).min() * 0.7, f(r).max() * 1.5)
        for row, (lo, hi) in ((1, (0.1, 100)), (2, (3, 9))):
            ax = axes[row, c]; m = (Ec >= lo) & (Ec <= hi)
            ax.plot(Ec[m], 100 * (p / r - 1)[m], color=C_PIX, lw=1, label="reltrans now")
            ax.plot(Ec[m], 100 * (b / r - 1)[m], color=C_B1, lw=1.2, label="B1 line cache")
            ax.axhspan(-0.02, 0.02, color="#d8d7d2", alpha=0.6, lw=0)
            ax.set_ylabel("residual [%]" + ("" if row == 1 else "\nFe K zoom")); style(ax)
            if row == 2: ax.set_xlabel("energy [keV]"); ax.set_xscale("linear")
            ea, eb = np.max(np.abs(p / r - 1)[m]), np.max(np.abs(b / r - 1)[m])
            ax.text(0.01, 0.95, f"max |res|: reltrans {100*ea:.2f}%, B1 {100*eb:.3f}%", transform=ax.transAxes, va="top", fontsize=8.5, color=INK)
    axes[0, 0].legend(frameon=False, fontsize=9, loc="lower left")
    fig.text(0.5, 0.005, "Grey band: ±0.02 % (the test suite's rtol 2e-4). Same xillver convolution for all three; only the reflection kernel differs.",
             ha="center", fontsize=8.5, color=INK)
    fig.suptitle(title, fontsize=13); fig.tight_layout(rect=(0, 0.02, 1, 0.97)); fig.savefig(fname, dpi=110)

fig_timeavg("refl", "Time-averaged reflection spectrum (boost = −1): current pixel grid vs B1", "spectra_refl.png")
fig_timeavg("dc", "Time-averaged total spectrum: current pixel grid vs B1", "spectra_total.png")

# cross spectrum, AGN-like high frequency (M = 4.6e7, 5e-5 - 1.5e-4 Hz)
fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), gridspec_kw=dict(height_ratios=[2, 1]))
for c, (fl, ylab, conv) in enumerate([("re", "Re G · E²  [arb.]", lambda x: x / dE * Ec**2),
                                      ("im", "Im G · E²  [arb.]", lambda x: x / dE * Ec**2),
                                      ("lag", "lag [c/Rg units × M]  (s)", lambda x: x / dE)]):
    d = np.load("results/i70_hi.npz")
    r, p, b = (conv(d[f"{k}/{fl}"]) for k in (REF, PIX, B1))
    ax = axes[0, c]
    ax.plot(Ec, r, color=C_REF, lw=2.5, label="converged reference")
    ax.plot(Ec, p, color=C_PIX, lw=1, label="reltrans now")
    ax.plot(Ec, b, color=C_B1, lw=1, ls="--", label="B1 line cache")
    ax.set_ylabel(ylab.replace(" [c/Rg units × M]", "")); style(ax); ax.set_title({"re": "Real part", "im": "Imaginary part", "lag": "Lag-energy spectrum"}[fl])
    ax = axes[1, c]; sc = np.max(np.abs(r))
    ax.plot(Ec, 100 * (p - r) / sc, color=C_PIX, lw=1); ax.plot(Ec, 100 * (b - r) / sc, color=C_B1, lw=1.2)
    ax.axhspan(-0.02, 0.02, color="#d8d7d2", alpha=0.6, lw=0)
    ax.set_ylabel("difference [% of max]"); ax.set_xlabel("energy [keV]"); style(ax)
axes[0, 0].legend(frameon=False, fontsize=9)
fig.suptitle("Cross spectrum, a = 0.998, i = 70°, M = 4.6e7, 5e-5–1.5e-4 Hz (reference band 0.3–10 keV, NICER response)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig("spectra_cross.png", dpi=110)
print("ok")
