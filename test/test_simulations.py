import pytest
import numpy as np

from pyreltrans import ReltransSimulator, Simrelt_Parameters


def test_basic_simulation(telescope):
    # Currently this is a smoke test to check that it does something vaguely
    # correct.
    energy = np.geomspace(0.2, 10.0, 50)
    simulator = ReltransSimulator(
        telescope.rmf_path, telescope.arf_path, telescope.bkg_path
    )
    result = simulator.simulate_lag(energy, Simrelt_Parameters())

    assert result.lag_true[0] == pytest.approx(-4.8397069, abs=1e-4)
