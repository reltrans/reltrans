import shutil
import pathlib
import dataclasses
from functools import partial

import pytest
import numpy as np

from pyreltrans import ReltransSimulator, Simrelt_Parameters

ENERGY = np.geomspace(0.2, 10.0, 50)
SIMULATION_FILES = ("sim_test.dat", "sim_cross_test.dat", "sim_xspec_test.dat")


@pytest.fixture
def simulator(telescope, tmp_path, monkeypatch, envars):
    monkeypatch.chdir(tmp_path)
    simulator = ReltransSimulator(
        telescope.rmf_path, telescope.arf_path, telescope.bkg_path
    )
    monkeypatch.setattr(
        simulator, "_setup_environ", partial(simulator._setup_environ, envars)
    )
    simulator._reltrans.reset()
    yield simulator
    simulator._reltrans.reset()


def test_basic_simulation(simulator):
    # Currently this is a smoke test to check that it does something vaguely
    # correct.
    result = simulator.simulate_lag(ENERGY, Simrelt_Parameters())

    assert result.background_count_rate == pytest.approx(0.03639, abs=1e-3)
    assert result.source_count_rate == pytest.approx(270.49288692412586, abs=1e-3)
    assert result.lag_true[0] == pytest.approx(-5.011609443771324, abs=1e-4)


def test_simulation_outputs(simulator):
    simulator.seed = 42
    result = simulator.simulate_lag(ENERGY, Simrelt_Parameters(), remove_files=False)
    for name in SIMULATION_FILES:
        assert pathlib.Path(name).is_file(), name

    widths = np.diff(ENERGY)
    np.testing.assert_allclose(result.energy, ENERGY[:-1] + 0.5 * widths)
    np.testing.assert_allclose(result.d_energy, 0.5 * widths)
    np.testing.assert_array_equal(result.std, result.d_lag)
    assert np.all(result.d_lag > 0)
    np.testing.assert_allclose(result.xspec_data[:, 0], ENERGY[:-1])
    np.testing.assert_allclose(result.xspec_data[:, 1], ENERGY[1:])
    np.testing.assert_allclose(result.xspec_data[:, 2], widths * result.lag, rtol=1e-5)
    np.testing.assert_allclose(result.xspec_data[:, 3], widths * result.d_lag, rtol=1e-5)

    # Reference lags in seconds for the seeded NICER simulation.
    np.testing.assert_allclose(
        result.lag[[0, 12, 24, 36, 48]],
        [-0.94530964, 19.64631736, -14.33964235, -6.51757300, 32.25010369],
        rtol=2e-4,
    )


def test_simulation_parameters(simulator):
    params = Simrelt_Parameters()
    base = simulator.simulate_lag(ENERGY, params)
    for name in SIMULATION_FILES:
        assert not pathlib.Path(name).exists(), name

    def simulate(**changes):
        return simulator.simulate_lag(ENERGY, dataclasses.replace(params, **changes))

    longer = simulate(exposure_time=4 * params.exposure_time)
    np.testing.assert_allclose(longer.d_lag, 0.5 * base.d_lag, rtol=1e-6)

    incoherent = simulate(coherence_squared=0.5)
    assert np.all(incoherent.d_lag > base.d_lag)

    stronger = simulate(power=2 * params.power)
    assert stronger.fractional_rms_sq_per_hz == pytest.approx(
        2 * base.fractional_rms_sq_per_hz
    )
    assert np.all(stronger.d_lag < base.d_lag)

    # None of these change the model or the count rates
    for result in (longer, incoherent, stronger):
        np.testing.assert_allclose(result.lag_true, base.lag_true, rtol=1e-4)
        assert result.source_count_rate == pytest.approx(base.source_count_rate)
        assert result.background_count_rate == pytest.approx(
            base.background_count_rate
        )

    brighter = simulate(Anorm=2 * params.Anorm)
    assert brighter.source_count_rate == pytest.approx(2 * base.source_count_rate)
    assert brighter.fractional_rms_sq_per_hz == pytest.approx(
        base.fractional_rms_sq_per_hz
    )

    simulator.set_refband_energies(0.5, 2.0)
    narrower = simulator.simulate_lag(ENERGY, params)
    assert narrower.source_count_rate < base.source_count_rate
    assert narrower.background_count_rate < base.background_count_rate


def test_simulation_seed(simulator):
    params = Simrelt_Parameters()
    simulator.seed = 123456
    first = simulator.simulate_lag(ENERGY, params)
    simulator.seed = 123457
    different = simulator.simulate_lag(ENERGY, params)
    simulator.seed = 123456
    repeated = simulator.simulate_lag(ENERGY, params)

    # An odd number of bins leaves a spare Gaussian deviate after each run.
    # Compare the noise separately from small cold-cache differences in lag_true.
    first_noise = (first.lag - first.lag_true) / first.d_lag
    different_noise = (different.lag - different.lag_true) / different.d_lag
    repeated_noise = (repeated.lag - repeated.lag_true) / repeated.d_lag
    np.testing.assert_allclose(repeated_noise, first_noise, rtol=0, atol=1e-10)
    assert np.all(np.isfinite(different_noise))
    assert not np.allclose(different_noise, first_noise)


@pytest.mark.skipif(shutil.which("flx2xsp") is None, reason="needs HEASoft flx2xsp")
def test_make_xspec(simulator):
    result = simulator.simulate_lag(ENERGY, Simrelt_Parameters())
    result.make_XSPEC("mock")
    assert pathlib.Path("mock.pha").is_file()
    assert pathlib.Path("mock.rsp").is_file()
