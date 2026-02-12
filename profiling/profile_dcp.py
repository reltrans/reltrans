"""
Structured profiling script for reltrans DCP evaluation.

This script focuses on steady-state performance by ignoring the first
evaluation (which may include FFTW initialization cost).
"""

import sys
import pathlib
import cProfile
import pstats
import io
import numpy as np

# Add project root to path
project_root = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

from test.wrapper import Reltrans, DCP_Parameters


def run_profile(n_runs=5):
    energy = np.logspace(-1, 2, 400).astype(np.float32)
    params = DCP_Parameters()

    try:
        rt = Reltrans()
    except RuntimeError as e:
        print("Native library not available.")
        print("Please build reltrans and ensure HEASoft is configured.")
        sys.exit(1)

    # Warmup (ignore FFTW initialization overhead)
    print("Running warmup evaluation...")
    rt.dcp(energy, params)

    print(f"Profiling steady-state over {n_runs} runs...")

    profiler = cProfile.Profile()
    profiler.enable()

    for _ in range(n_runs):
        rt.dcp(energy, params)

    profiler.disable()

    s = io.StringIO()
    stats = pstats.Stats(profiler, stream=s).sort_stats("cumtime")
    stats.print_stats(20)

    print("\nTop 20 functions by cumulative time:\n")
    print(s.getvalue())


if __name__ == "__main__":
    run_profile()

