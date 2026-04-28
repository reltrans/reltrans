import numpy as np
import pytest
from wrapper import DCP_Parameters, Reltrans


@pytest.fixture(scope="module")
def reltrans_instance():
    """Fresh Reltrans instance, independent of the session-scoped fixture."""
    return Reltrans()


@pytest.fixture(scope="module")
def energy():
    """Coarse energy grid, fast for unit tests."""
    return np.logspace(np.log10(0.3), np.log10(10.0), 41, dtype=np.float32)


@pytest.fixture(scope="module")
def reset_fn(reltrans_instance):
    """
    Resolves reset_reltrans() from the compiled library under either symbol
    naming convention (standalone subroutine vs. gfortran module member).
    Skips dependent tests with a clear message if the symbol is absent.
    """
    for sym in ("reset_reltrans_", "__m_genreltrans_MOD_reset_reltrans"):
        try:
            fn = getattr(reltrans_instance.lib_reltrans, sym)
            fn.argtypes = []
            fn.restype = None
            return fn
        except AttributeError:
            continue
    pytest.skip("'reset_reltrans' symbol not found")


# Tests


def test_reset_produces_valid_output(reltrans_instance, reset_fn, energy, assert_snapshot):
    """
    dcp() immediately after a single reset() must return the same output as a
    known-good snapshot.
    """
    reset_fn()
    out = reltrans_instance.dcp(energy, DCP_Parameters(flo_hz=0.0, fhi_hz=0.0))
    assert_snapshot(out)


def test_aba_with_reset_matches_original(reltrans_instance, reset_fn, energy):
    """
    Core Issue #124 regression: dcp(A) → dcp(B) → reset() → dcp(A)
    must reproduce A's output exactly.

    After calling with parameters B (different from A), the cache is polluted.
    reset() should clear that cached state so the next call with A reproduces
    the original result.
    """
    p_a = DCP_Parameters(h=4.0, a=0.3, inc=30.0, flo_hz=0.0, fhi_hz=0.0)
    p_b = DCP_Parameters(h=40.0, a=0.9, inc=60.0, flo_hz=0.0, fhi_hz=0.0)

    reference = reltrans_instance.dcp(energy, p_a).copy()

    # Pollute cache with different parameters
    reltrans_instance.dcp(energy, p_b)

    # Reset and re-run with original parameters
    reset_fn()
    result = reltrans_instance.dcp(energy, p_a)

    np.testing.assert_allclose(
        result, reference,
        rtol=1e-4, atol=1e-10,
        err_msg="ABA failed: reset() did not clear cached state.",
    )


def test_rev_nosav_aba_is_consistent(reltrans_instance, energy, envars):
    """
    REV_NOSAV=1 must force a full recompute on every call, making the ABA
    pattern consistent without any explicit reset() call.

    Without REV_NOSAV, calling dcp(A) then dcp(B) then dcp(A) may return a
    stale cached result for the second A call. With REV_NOSAV=1, each call
    fully recomputes, so both A calls must agree.
    """
    envars["REV_NOSAV"] = "1"

    p_a = DCP_Parameters(h=4.0,  a=0.3, flo_hz=0.0, fhi_hz=0.0)
    p_b = DCP_Parameters(h=40.0, a=0.9, flo_hz=0.0, fhi_hz=0.0)

    out_a1 = reltrans_instance.dcp(energy, p_a).copy()
    reltrans_instance.dcp(energy, p_b)
    out_a2 = reltrans_instance.dcp(energy, p_a).copy()

    np.testing.assert_allclose(
        out_a2, out_a1,
        rtol=1e-4, atol=1e-10,
        err_msg="ABA inconsistent with REV_NOSAV=1: force_recompute path not working.",
    )
