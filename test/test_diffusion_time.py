"""Tests for reltransDT, the reltransDCp flavour with the diffusion
(reprocessing) time of the reflected photons."""
import dataclasses

import numpy as np
import pytest
from pyreltrans import DCP_Parameters, DT_Parameters

ENERGY = np.logspace(np.log10(0.1), np.log10(100), 301)
E_MID = 0.5 * (ENERGY[1:] + ENERGY[:-1])

# Single frequency bin centred at sqrt(1.0 * 1.2) Hz (nf = 1 for dlogf = 0.09)
XRB_NARROW = DT_Parameters(mass=10.0, flo_hz=1.0, fhi_hz=1.2)
F_CENTRE_NARROW = np.sqrt(1.2)
XRB_LAG = DT_Parameters(mass=10.0, flo_hz=1.0, fhi_hz=2.0, re_im=-4.0)


@pytest.fixture
def reference_band(envars):
    """Set the reference band used by the unfolded (negative re_im) modes."""
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"


def _max_rel_diff(a, b):
    return np.max(np.abs(a - b)) / np.max(np.abs(b))


@pytest.mark.parametrize("re_im", [-1.0, -2.0])
def test_diffusion_time_full_period(reltrans, reference_band, re_im):
    """A delay of exactly one period (tau = 1/f) is a 2 pi phase, so the
    per-frequency transform of the complex spectrum must reproduce the tau = 0
    result to round-off. This checks that the phase uses the kernel
    frequency."""
    no_delay = dataclasses.replace(XRB_NARROW, re_im=re_im)
    expected = reltrans.dt(ENERGY, no_delay)

    reltrans.reset()
    full_period = dataclasses.replace(no_delay, diff_time=1.0 / F_CENTRE_NARROW)
    output = reltrans.dt(ENERGY, full_period)
    assert _max_rel_diff(output, expected) < 1e-10

    reltrans.reset()
    quarter_period = dataclasses.replace(no_delay,
                                         diff_time=0.25 / F_CENTRE_NARROW)
    output = reltrans.dt(ENERGY, quarter_period)
    assert _max_rel_diff(output, expected) > 1e-2


@pytest.mark.parametrize("re_im", [-1.0, -2.0, -4.0])
def test_diffusion_time_negligible(reltrans, reference_band, re_im):
    """A negligible delay goes through the per-frequency (complex spectrum)
    path but must agree with tau = 0."""
    no_delay = dataclasses.replace(XRB_LAG, re_im=re_im)
    expected = reltrans.dt(ENERGY, no_delay)
    reltrans.reset()
    output = reltrans.dt(ENERGY, dataclasses.replace(no_delay,
                                                      diff_time=1e-12))
    assert _max_rel_diff(output, expected) < 1e-8


@pytest.mark.parametrize("flo_hz, fhi_hz, re_im", [(0.0, 0.0, 1.0),
                                                   (1.0, 2.0, -1.0),
                                                   (1.0, 2.0, -2.0),
                                                   (1.0, 2.0, -4.0)])
def test_reltransDT_zero_diffusion_equals_reltransDCp(reltrans, reference_band,
                                                      flo_hz, fhi_hz, re_im):
    """With diff_time = 0, reltransDT is reltransDCp."""
    dcp = DCP_Parameters(mass=10.0, flo_hz=flo_hz, fhi_hz=fhi_hz, re_im=re_im)
    expected = reltrans.dcp(ENERGY, dcp)
    reltrans.reset()
    output = reltrans.dt(ENERGY, DT_Parameters(**dataclasses.asdict(dcp)))
    np.testing.assert_array_equal(output, expected)


def test_diffusion_time_time_averaged_unchanged(reltrans):
    """The delay has no effect on the time-averaged spectrum (f = 0)."""
    expected = reltrans.dt(ENERGY, DT_Parameters(mass=10.0))
    reltrans.reset()
    output = reltrans.dt(ENERGY, DT_Parameters(mass=10.0, diff_time=5e-3))
    np.testing.assert_array_equal(output, expected)


def test_diffusion_time_cache(reltrans, reference_band):
    """Changing only diff_time between calls must redo the convolutions, and
    going back to the previous value must give the previous result."""
    first = reltrans.dt(ENERGY, XRB_LAG)
    delayed = reltrans.dt(ENERGY, dataclasses.replace(XRB_LAG,
                                                       diff_time=5e-3))
    back = reltrans.dt(ENERGY, XRB_LAG)
    assert not np.allclose(delayed, first)
    np.testing.assert_array_equal(back, first)

    reltrans.reset()
    fresh = reltrans.dt(ENERGY, dataclasses.replace(XRB_LAG, diff_time=5e-3))
    np.testing.assert_allclose(delayed, fresh, rtol=0, atol=1e-15)


@pytest.mark.parametrize("diff_time", [1e-3, 5e-3, 13e-3])
def test_diffusion_time_lag_sanity(reltrans, reference_band, diff_time):
    """Physical sanity of the lag-energy spectrum: the extra lag can never
    exceed the diffusion time (the reflection is diluted by the continuum),
    and the reflection-dominated Compton hump lags more than the
    continuum-dominated 2-3 keV band."""
    dE = np.diff(ENERGY)
    lag0 = reltrans.dt(ENERGY, XRB_LAG) / dE
    reltrans.reset()
    lag = reltrans.dt(ENERGY, dataclasses.replace(XRB_LAG,
                                                   diff_time=diff_time)) / dE
    delta_lag = lag - lag0
    assert np.all(np.isfinite(lag))
    assert np.max(np.abs(delta_lag)) <= diff_time
    hump = np.argmin(np.abs(E_MID - 20.0))
    continuum = np.argmin(np.abs(E_MID - 2.5))
    assert delta_lag[hump] > delta_lag[continuum]
