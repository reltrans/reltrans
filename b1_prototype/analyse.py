import sys, glob, os, re
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
T_RAY, T_Q, T_PH, T_ADD = 3.8e-6, 0.41e-6, 15e-9, 9e-9

def err(fl, x, r):
    if fl in ("im", "lag"):
        return np.max(np.abs(x - r)) / np.max(np.abs(r))
    den = np.maximum(np.abs(r), 1e-3 * np.max(np.abs(r)))
    return np.max(np.abs(x - r) / den)

def family(name):
    return name.split("_")[0] if not name.startswith("b1gl") else "b1gl"

def cost_model(name, c):
    rays, nq, ntri, ndep, nf, tpy = c
    tq = T_Q if not name.startswith("ctf") else 0.1e-6
    tint = 20e-9 * ntri / 2 if name.startswith("b1c") else 0.0   # spectral upsampling of node values
    t_q = nq * tq + tint
    return rays * T_RAY + t_q + nf * (ntri * T_PH + ndep * T_ADD), rays * T_RAY, t_q, nf * (ntri * T_PH + ndep * T_ADD)

def load_all():
    rows = []
    for path in sorted(glob.glob(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "results", "*.npz"))):
        geom, fset = os.path.basename(path)[:-4].rsplit("_", 1)
        d = dict(np.load(path))
        refname = "b1up_256x8_U8x8"
        fls = sorted({k.split("/")[1] for k in d if k.startswith(refname + "/") and not k.endswith("cost")})
        names = sorted({k.split("/")[0] for k in d if k.endswith("/cost")})
        for n in names + ["reltrans"]:
            for fl in fls:
                key = f"{n}/{fl}"
                if key not in d: continue
                e = err(fl, d[key], d[f"{refname}/{fl}"])
                c = d.get(f"{n}/cost")
                if n == "reltrans": c = d["pix_200/cost"] if "pix_200/cost" in d else None
                cm = cost_model(n, c) if c is not None else (np.nan,) * 4
                rows.append(dict(geom=geom, fset=fset, method=n, fam=family(n) if n != "reltrans" else "reltrans",
                                 fl=fl, err=e, t=cm[0], t_ray=cm[1], t_q=cm[2], t_dep=cm[3],
                                 rays=c[0] if c is not None else np.nan))
    return rows

if __name__ == "__main__":
    rows = load_all()
    import collections
    by = collections.defaultdict(list)
    for r in rows: by[(r["geom"], r["fset"], r["method"])].append(r)
    for (g, f, m), rs in sorted(by.items()):
        print("%-10s %-3s %-22s rays %8s  t_model %7.1f ms  " % (g, f, m, int(rs[0]["rays"]) if rs[0]["rays"] == rs[0]["rays"] else "-", rs[0]["t"] * 1e3) +
              "  ".join("%s %.1e" % (r["fl"], r["err"]) for r in rs))
