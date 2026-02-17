import pytest
import numpy as np
try:
    from wrapper import DCP_Parameters
    WRAPPER_AVAILABLE = True
except ModuleNotFoundError:
    WRAPPER_AVAILABLE = False

pytestmark = [
    pytest.mark.integration,
    pytest.mark.regression,
]

if not WRAPPER_AVAILABLE:
    pytest.skip("Native wrapper not available.", allow_module_level=True)

def _debug_plot(energy, output, title=""):
    """Used when creating new tests for quickly looking at the data to make
    sure it is sensible."""
    import matplotlib.pyplot as plt

    plt.clf()
    plt.plot(energy[0:-1], output)
    plt.xscale("log")
    plt.show()


def test_basic_invocation(reltrans, assert_snapshot):
    """A smoke test to check whether the default values are working."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters())
    assert_snapshot(output)


def test_re_im_parameter(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    xrb1 = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="time_lag")

    xrb1.re_im = 1
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="real_part")

    xrb1.re_im = 2
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="imaginary_part")

    xrb1.re_im = 3
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="magnitude")
