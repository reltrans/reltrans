import pytest
import numpy as np
from wrapper import DCP_Parameters, rtdist_Parameters


def _debug_plot(energy, output, title="", xlabel="", ylabel="", yscale="linear", xscale="log"):
    """Used when creating new tests for quickly looking at the data to make
    sure it is sensible."""
    import matplotlib.pyplot as plt

    plt.clf()
    plt.plot(energy[0:-1], output)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

def test_basic_invocation(reltrans, assert_snapshot):
    """A smoke test to check whether the default values are working."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters())
    # _debug_plot(energy,output, "reltransDCp time-averaged spectrum [default parameters]")
    assert_snapshot(output)

def test_basic_invocation_rtdist(reltrans, assert_snapshot, envars):
    """A smoke test to check whether the default values are working."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.rtdist(energy, rtdist_Parameters())
    # _debug_plot(energy,output, "rtdist time-averaged spectrum [default parameters]")
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

def test_re_im_rtdist(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    RTDIST model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    agn1 = rtdist_Parameters(mass=4.7e6, flo_hz=1e-5, fhi_hz=1e-4, re_im=4.0)
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title = "rtdist lag spectrum", xlabel="Energy [keV]")
    assert_snapshot(output, name="time_lag", rtol=1e-3)

    agn1.re_im = 3
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist modulus spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="magnitude", rtol=1e-3)

    agn1.re_im = 1
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist real part spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="real_part", rtol=1e-3)

    agn1.re_im = 2
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist imaginary part spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="imaginary_part", rtol=1e-3)

    
