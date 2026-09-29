import pytest
import numpy as np
from pyreltrans import PL_Parameters

from conftest import _get_snapshot

plot_spectral_kwargs = dict(
        yscale = "log", xlabel = "Energy", ylabel = "Flux" , xscale = "log", units = "ef"
)


@pytest.mark.rt
def test_basic_invocation_reltransPL(reltrans, assert_snapshot, envars):
    """A smoke test to check if reltransPL does NOT return NAN."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.pl(energy, PL_Parameters(h=10.0, boost=-1.0))
    np.testing.assert_equal(output, output)

    
def test_re_im_reltransPL(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    reltransPL are NOT NAN. This test does NOT compare with benchmarks values
    This test requires an RMF and ARF, which is provided by the `telescope` fixture."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)
    dE = (energy[1:] - energy[:-1])

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    xrb1 = PL_Parameters(mass=10.0, flo_hz=0.1, fhi_hz=0.2, re_im=4.0)
    output = reltrans.pl(energy, xrb1)
    np.testing.assert_equal(output, output)

    xrb1.re_im = 1
    output = reltrans.pl(energy, xrb1)
    np.testing.assert_equal(output, output)

    xrb1.re_im = 2
    output = reltrans.pl(energy, xrb1)
    np.testing.assert_equal(output, output)

    xrb1.re_im = 3
    output = reltrans.pl(energy, xrb1)
    np.testing.assert_equal(output, output)

    xrb1.re_im = 6.0
    output = reltrans.pl(energy, xrb1)/dE
    np.testing.assert_equal(output, output)

    xrb1.re_im = 5.0
    output = reltrans.pl(energy, xrb1)
    np.testing.assert_equal(output, output)
