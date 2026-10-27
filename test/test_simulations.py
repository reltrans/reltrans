import os
import shutil
import pathlib
import dataclasses

import pytest
import numpy as np

from pyreltrans import ReltransSimulator, Simrelt_Parameters

ENERGY = np.geomspace(0.2, 10.0, 50)


@pytest.fixture
def simulator(telescope, tmp_path, monkeypatch):
    # The simulation writes its files to the working directory and sets
    # environment variables without restoring them.
    monkeypatch.chdir(tmp_path)
    environ = dict(os.environ)
    yield ReltransSimulator(
        telescope.rmf_path, telescope.arf_path, telescope.bkg_path
    )
    os.environ.clear()
    os.environ.update(environ)


def test_basic_simulation(telescope):
    # Currently this is a smoke test to check that it does something vaguely
    # correct.
    energy = np.geomspace(0.2, 10.0, 50)
    simulator = ReltransSimulator(
        telescope.rmf_path, telescope.arf_path, telescope.bkg_path
    )
    result = simulator.simulate_lag(energy, Simrelt_Parameters())

    assert result.background_count_rate == pytest.approx(0.03639, abs=1e-3)
    assert result.source_count_rate == pytest.approx(270.49288692412586, abs=1e-3)
    assert result.lag_true[0] == pytest.approx(-5.011609443771324, abs=1e-4)


def test_simulation_outputs(simulator):
    result = simulator.simulate_lag(ENERGY, Simrelt_Parameters(), remove_files=False)
    for name in ("sim_test.dat", "sim_cross_test.dat", "sim_xspec_test.dat"):
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

    # The simulated lags scatter around the true lag with width d_lag
    pull = (result.lag - result.lag_true) / result.d_lag
    assert abs(np.mean(pull)) < 0.5
    assert 0.6 < np.std(pull) < 1.4


def test_simulation_parameters(simulator):
    params = Simrelt_Parameters()
    base = simulator.simulate_lag(ENERGY, params)
    assert not any(pathlib.Path().iterdir())

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


@pytest.mark.xfail(strict=True, reason="simrelt_extra ignores SEED_SIM")
def test_simulation_seed(simulator):
    first = simulator.simulate_lag(ENERGY, Simrelt_Parameters())
    second = simulator.simulate_lag(ENERGY, Simrelt_Parameters())
    np.testing.assert_array_equal(first.lag, second.lag)


@pytest.mark.skipif(shutil.which("flx2xsp") is None, reason="needs HEASoft flx2xsp")
def test_make_xspec(simulator):
    result = simulator.simulate_lag(ENERGY, Simrelt_Parameters())
    result.make_XSPEC("mock")
    assert pathlib.Path("mock.pha").is_file()
    assert pathlib.Path("mock.rsp").is_file()
