# Performance Profiling Guide

This document describes how to profile reltrans model evaluations in a
structured and reproducible way.

## Why steady-state profiling?

The first call to `rt.dcp()` may include FFTW initialization overhead.
To measure realistic runtime performance, we:

1. Perform one warmup evaluation.
2. Profile only subsequent evaluations.
3. Focus on cumulative time (`cumtime`) to identify bottlenecks.

## Running the profiler

From the project root:

```bash
python profiling/profile_dcp.py

