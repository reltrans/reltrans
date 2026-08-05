import os
import sys
import dataclasses
import subprocess

import numpy as np

from pyreltrans.wrappers import Reltrans, Simrelt_Parameters, DCP_Parameters


@dataclasses.dataclass
class SimulationResult:
    """
    The result of a simulation.
    """

    # Energy in keV.
    energy: np.ndarray
    # Energy bin width in keV.
    d_energy: np.ndarray
    # The simulated lag data point.
    lag: np.ndarray
    # The simulated lag error.
    d_lag: np.ndarray

    # The true lag without any scatter.
    lag_true: np.ndarray
    # The standard deviation of the lag.
    std: np.ndarray

    # The count rate.
    count_rate: np.ndarray

    # The real part of the cross spectrum.
    real: np.ndarray

    # The imaginary part of the cross spectrum.
    imag: np.ndarray

    # The data used to create an XSPEC mock spectrum.
    xspec_data: np.ndarray

    def make_XSPEC(self, name: str):
        """
        Run the requiste `flx2xsp` command to generate a mock XSPEC spectrum
        with a dummy response file. The files will be named `{name}.pha` and
        `{name}.rsp`.

        Requires HEASoft on the path.
        """
        # Write the XSPEC file:
        tmp_xspec_file = "tmp_xspec.dat"
        np.savetxt(tmp_xspec_file, self.xspec_data, delimiter="\t")
        subprocess.run(
            ["flx2xsp", tmp_xspec_file, f"{name}.pha", f"{name}.rsp"], check=True
        )
        # Remove the temporary file again:
        os.remove(tmp_xspec_file)

    @staticmethod
    def from_file(path: str, cross_path: str, xspec_path: str) -> "SimulationResult":
        """
        Read the simulation result from a set of paths.

        The `path` should point to the simulation data, and `cross_path` to the
        cross spectrum file. Finally `xspec_path` should point to the file
        containing the data for the XSPEC simulation.

        TODO: at some point simrelt should return the arrays to Python to avoid
        having to read them in from files.
        """
        data = np.genfromtxt(path, skip_header=2)
        cross_data = np.genfromtxt(cross_path)
        xspec_data = np.genfromtxt(xspec_path)

        E = data[:, 0]
        d_E = data[:, 1]
        lag_simulated = data[:, 2]
        d_lag_simulated = data[:, 3]
        lag_true = data[:, 4]

        return SimulationResult(
            energy=E,
            d_energy=d_E,
            lag=lag_simulated,
            d_lag=d_lag_simulated,
            lag_true=lag_true,
            count_rate=cross_data[:, 2],
            real=cross_data[:, 3],
            imag=cross_data[:, 4],
            std=cross_data[:, 5],
            xspec_data=xspec_data,
        )


class ReltransSimulator:
    """
    A simulation wrapper around the reltrans library that resets all the caches
    between calls. This **does not** make it thread safe but does mean seperate
    instances that simulate different instruments.
    """

    def __init__(self, rmf: str, arf: str, background: str):
        if os.environ.get("RELTRANS_TABLES", None) is None:
            raise Exception(
                "Must set the RELTRANS_TABLES environment variable before running the simulation."
            )

        self.arf = arf
        self.rmf = rmf
        self.background = background
        self.e_low = 0.2
        self.e_high = 10.0
        self.seed = 42
        self._reltrans = Reltrans()

    def _setup_environ(self):
        # Set the environment variables
        os.environ["SEED_SIM"] = f"{self.seed}"
        os.environ["ARF_SET"] = os.path.abspath(self.arf)
        os.environ["RMF_SET"] = os.path.abspath(self.rmf)
        os.environ["EMIN_REF"] = str(self.e_low)
        os.environ["EMAX_REF"] = str(self.e_high)
        os.environ["BKG_SET"] = os.path.abspath(self.background)
        os.environ["BACKSCL"] = "1"

    def set_refband_energies(self, e_low, e_high):
        """
        Set the reference band energies.
        """
        self.e_low = e_low
        self.e_high = e_high

    def simulate_lag(
        self,
        energy: np.ndarray,
        params: Simrelt_Parameters,
        remove_files=True,
    ) -> SimulationResult:
        """
        Simulate a reverberation lag.

        If `remove_files` is True, then the simulation files written during the
        simulation are deleted again.
        """
        self._setup_environ()
        self._reltrans.simrelt(energy, params)
        result = SimulationResult.from_file(
            "sim_test.dat", "sim_cross_test.dat", "sim_xspec_test.dat"
        )
        # Cleanup the temporary files
        if remove_files:
            os.remove("sim_test.dat")
            os.remove("sim_cross_test.dat")
            os.remove("sim_xspec_test.dat")

        return result
