# B1 prototype for reltrans reflection kernels

Investigation code from 2026-09-25/26: can kerrzbb's boundary-fitted image-plane
method (B1) make reltrans's reflection transfer functions faster or more
accurate? This is a prototype, not production code. The samplers are in Python
and drive reltrans's own Fortran through a few C hooks; the downstream steps
(xillver convolution, cross spectra, folding) are unchanged reltrans.

Findings, tables and decisions are in the project doc
`claude/reltrans-b1-findings.md`. The figures in `figures/` are the ones
discussed there.

## Fortran changes (in `subroutines/`)

All are for investigation only and leave the model output unchanged (the test
suite passes, 41/41).

- `rt_timing.f90`: wall-clock section timers (`tic`/`toc`), read from C with
  `rt_timing_get` and cleared with `rt_timing_reset`. Timed sections: total,
  rtrans, trace, getdcos, GR and flat summations, getlens, convolutions
  (split into xillver `rest_frame` and FFT convolution), init_cont, post-processing.
- `rt_state.f90`: module holding a copy of the model arguments and the last
  kernels computed by `rtrans`, plus an optional externally supplied kernel.
- `rt_hooks.f90`: C-callable routines used by the Python harness:
  - `rt_trace_batch`: trace a batch of image-plane rays with kerrz;
  - `rt_sample_quantities`: per-sample g, emissivity weight, delay, log gsd
    and emission angle, reproducing `sum_multiple_lampposts` for one lamppost
    (including its flat-space delay formula);
  - `rt_flat_batch`: the straight-line mapping of the outer disc;
  - `rt_state_info`, `rt_get_freqs`, `rt_get_kernel`: read back state;
  - `rt_set_kernel` / `rt_clear_kernel`: replace rtrans's W0..W3 by an
    external kernel before the convolutions (single lamppost, `MU_ZONES=1`).
- `strans.f90`: timers, and a call to `rt_store_and_override` after the
  summations.
- `genreltrans.f90`: timers.
- `common.f90`: environment overrides `RT_NRO`, `RT_NPHI`, `RT_NRON`,
  `RT_NPHIN` for the camera grid.
- `header.h`: includes the three new files.

## Python modules

| file | contents |
|---|---|
| `rtk.py` | `Harness` (ctypes wrapper of the hooks), `Kernel` accumulator, deposition rules: nearest bin, cloud-in-cell, and `deposit_tent` (exact triangle deposit) |
| `methods.py` | the current reltrans camera re-implemented (`pixel_samples`); reproduces reltrans's kernels to 2e-7 |
| `contour.py` | vectorised contour solver: image radius at which a ray lands on disc radius R, for each image angle |
| `b1.py` | plain B1 (Gauss nodes between zone-edge contours) and B1 with coarse triangles |
| `b1up.py` | B1 contours + Lobatto nodes per zone, spectral upsampling (Fourier in angle, polynomial across the zone), exact triangle deposit; also the flat-space outer region |
| `b1c.py` | B1 line cache: rays traced once per (a, i) between the ISCO and r = 300 contours; zone edges and nodes for any rin by interpolation; then as `b1up.py`. This is the recommended variant. |
| `ctf.py` | per-radius contours with the Cunningham g* parameterisation and exact deposit |
| `study.py`, `run_study.sh` | the error-against-cost study (5 geometries × 3 frequency sets × all methods) |
| `analyse.py`, `plot_cost.py` | cost model and `figures/kernel_cost.png` |
| `plot_spectra.py`, `viz_methods.py` | `figures/spectra_*.png`, `figures/methods_viz.png` |
| `rtprof.py`, `profile_r0.py`, `run_case.py` | R0 profiling and the pixel-resolution ladder |
| `validate_harness.py` | checks the Python pixel sampler against reltrans's own kernel |
| `gaptest.py` | shows the lag error from the dropped annulus at the GR/flat seam |
| `experiments/` | exploratory runs from the session (convergence of the variants, line-cache settings, outer-region resolution at high frequency) |
| `bench/depbench.f90` | Fortran micro-benchmark behind the deposit cost constants |

## Running

Build reltrans as usual (`make`), with `RELTRANS_TABLES` set. From this
directory:

```
export HEADAS=$(python3 -c "import xspectrampoline_helpers as h; print(h.get_HEADAS())")
mkdir -p Output                 # ReIm = 7 writes Output/Total.dat unconditionally
python3 validate_harness.py     # pixel sampler vs reltrans kernel
python3 study.py default dc     # one geometry / frequency set -> results/default_dc.npz
python3 analyse.py              # tabulate errors and modelled costs
```

`RELTRANS_ROOT` (default: the parent directory) and `RELTRANS_BUILD`
(default: `build`) select the library. The study results used for the figures
(about 3 MB of `.npz`) are not committed.

## Caveats

- Costs for the B1 variants are modelled from measured Fortran unit costs
  (3.8 µs per ray, 0.41 µs per per-sample quantity evaluation, 15 ns per
  triangle and 9 ns per bin deposit per frequency); they were not timed as
  Fortran implementations.
- Only a single lamppost and `MU_ZONES=1` are supported by the kernel hooks.
- Thin disc only in the outer region (`honr = 0`).
