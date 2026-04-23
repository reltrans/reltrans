import os
import numpy as np
import pytest
from wrapper import DCP_Parameters, Reltrans


@pytest.fixture(scope="module")
def rt():
    """Fresh Reltrans instance, independent of the session-scoped fixture."""
    return Reltrans()


@pytest.fixture(scope="module")
def energy():
    """Coarse energy grid, fast for unit tests."""
    return np.logspace(np.log10(0.3), np.log10(10.0), 41, dtype=np.float32)


@pytest.fixture(scope="module")
def reset_fn(rt):
    """
    Resolves reset_reltrans() from the compiled library under either symbol
    naming convention (standalone subroutine vs. gfortran module member).
    Skips dependent tests with a clear message if the symbol is absent.
    """
    for sym in ("reset_reltrans_", "__m_genreltrans_MOD_reset_reltrans"):
        try:
            fn = getattr(rt.lib_reltrans, sym)
            fn.argtypes = []
            fn.restype = None
            return fn
        except AttributeError:
            continue
    pytest.skip(
        "reset_reltrans symbol not found. Recompile: cd /workspace && make clean && make"
    )

# Tests

def test_reset_symbol_exists(rt):
    """reset_reltrans must be present in the compiled library."""
    found = any(
        hasattr(rt.lib_reltrans, sym)
        for sym in ("reset_reltrans_", "__m_genreltrans_MOD_reset_reltrans")
    )
    assert found, (
        "reset_reltrans not found in the library. "
        "Run: cd /workspace && make clean && make"
    )


def test_reset_produces_valid_output(rt, reset_fn, energy):
    """
    dcp() immediately after reset() must return finite, non-negative float32
    output of the expected shape. Covers: callable, idempotent (3 resets),
    and all basic output-validity properties in one pass.
    """
    for _ in range(3):
        reset_fn()

    out = rt.dcp(energy, DCP_Parameters(flo_hz=0.0, fhi_hz=0.0))

    assert out.dtype == np.float32,          f"Expected float32, got {out.dtype}"
    assert out.shape == (len(energy) - 1,),  f"Expected shape ({len(energy)-1},), got {out.shape}"
    assert np.all(np.isfinite(out)),         f"{(~np.isfinite(out)).sum()} NaN/Inf bins after reset()"
    assert np.all(out >= 0.0),               f"{(out < 0).sum()} negative bins in DC spectrum after reset()"


def test_aba_with_reset_matches_original(rt, reset_fn, energy):
    """
    Core Issue #124 regression: dcp(A) → dcp(B) → reset() → dcp(A)
    must reproduce A's output exactly.

    Also stress-tests with 5 different B parameters to confirm the
    cache-invalidation path is robust, not just lucky on one parameter set.
    """
    p_a = DCP_Parameters(h=4.0, a=0.3, inc=30.0, flo_hz=0.0, fhi_hz=0.0)
    reference = rt.dcp(energy, p_a).copy()

    for h_b in np.linspace(5.0, 50.0, 5):
        rt.dcp(energy, DCP_Parameters(h=float(h_b), a=0.9, flo_hz=0.0, fhi_hz=0.0))
        reset_fn()
        np.testing.assert_allclose(
            rt.dcp(energy, p_a), reference,
            rtol=1e-4, atol=1e-10,
            err_msg=f"ABA failed after B call with h={h_b:.1f}: reset() did not clear cached state.",
        )


def test_rev_nosav_aba_is_consistent(rt, energy):
    """
    REV_NOSAV=1 must force a full recompute on every call, making the ABA
    pattern consistent without any explicit reset() call.
    """
    os.environ["REV_NOSAV"] = "1"
    try:
        p_a = DCP_Parameters(h=4.0,  a=0.3, flo_hz=0.0, fhi_hz=0.0)
        p_b = DCP_Parameters(h=40.0, a=0.9, flo_hz=0.0, fhi_hz=0.0)

        out_a1 = rt.dcp(energy, p_a).copy()
        rt.dcp(energy, p_b)
        out_a2 = rt.dcp(energy, p_a).copy()

        np.testing.assert_allclose(
            out_a2, out_a1,
            rtol=1e-4, atol=1e-10,
            err_msg="ABA inconsistent with REV_NOSAV=1: force_recompute path not working.",
        )
    finally:
        os.environ.pop("REV_NOSAV", None)
