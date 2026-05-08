import pytest
import numpy as np
from wrapper import DCP_Parameters, Dbl_Parameters, rtdist_Parameters


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
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters())
    # _debug_plot(energy,output, "reltransDCp time-averaged spectrum [default parameters]")
    assert_snapshot(output)

def test_basic_absorption_invocation(reltrans, assert_snapshot):
    """A smoke test to check whether absorption is being correctly applied."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(nh = 0.2))
    # _debug_plot(energy,output, "reltransDCp time-averaged spectrum [default parameters]")
    assert_snapshot(output)
    
    
def test_re_im_parameter(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    reltrans.reset()
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

    
def test_basic_invocation_reltransDbl(reltrans, assert_snapshot, envars):
    """A smoke test to check whether the default values are working."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dbl_lamp(energy, Dbl_Parameters())
    # _debug_plot(energy,output, "rtdist time-averaged spectrum [default parameters]", xlabel="Energy [keV]", xscale='log', yscale='log')
    assert_snapshot(output)

    
def test_re_im_reltransDbl(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    reltransDbl model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    xrb2 = Dbl_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)
    output = reltrans.dbl_lamp(energy, xrb2)
    # _debug_plot(energy,output, title = "reltransDbl lag spectrum", xlabel="Energy [keV]")
    assert_snapshot(output, name="time_lag")

    xrb2.re_im = 3
    output = reltrans.dbl_lamp(energy, xrb2)
    # _debug_plot(energy,output, title="reltransDbl imaginary part spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="magnitude")

    xrb2.re_im = 1
    output = reltrans.dbl_lamp(energy, xrb2)
    # _debug_plot(energy,output, title="reltransDbl modulus spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="real_part")

    xrb2.re_im = 2
    output = reltrans.dbl_lamp(energy, xrb2)
    # _debug_plot(energy,output, title="reltransDbl real part spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="imaginary_part")

def test_reltransDCp_called_twice_no_resetting(reltrans, assert_snapshot):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_dcp = reltrans.dcp(energy, DCP_Parameters())
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters())
    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)
    
def test_reltransDCp_called_twice_after_resetting(reltrans, assert_snapshot):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_dcp = reltrans.dcp(energy, DCP_Parameters())
    reltrans.reset()
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters())
    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)

    
def test_resetting_between_flavours(reltrans):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_dcp = reltrans.dcp(energy, DCP_Parameters())
    reltrans.reset()
    output_rtdist = reltrans.rtdist(energy, rtdist_Parameters())
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters())
    reltrans.reset()
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters())
    reltrans.reset()
    output_rtdist2 = reltrans.rtdist(energy, rtdist_Parameters())
    reltrans.reset()
    output_dbl2 = reltrans.dbl_lamp(energy, Dbl_Parameters())

    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)
    np.testing.assert_allclose(output_rtdist, output_rtdist2, rtol=1e-4)
    np.testing.assert_allclose(output_dbl, output_dbl2, rtol=1e-4)


def test_strans_routines_getrgrid(reltrans, assert_snapshot):
    rnmin = 1.7266670949940284
    rnmax = 300.0
    mueff = 0.86602540378443871
    nro   = 200
    nphi  = 200
    rn , domega = reltrans.getrgrid(rnmin, rnmax, mueff, nro, nphi)
    assert_snapshot(rn, name="rn")
    assert_snapshot(domega, name="domega")


def test_strans_routines_grtrace_outputs(reltrans, assert_snapshot):
    reltrans.reset()

    rnmin = 1.7266670949940284
    rnmax = 300.0
    mueff = 0.86602540378443871
    nro   = 200
    nphi  = 200
    mu0 = mueff
    spin = 0.998
    rin  = 1.2369694567256766
    rout = 1000.0
    mudisk = 0.0
    distance = 18000000.0
    
    rn , domega = reltrans.getrgrid(rnmin, rnmax, mueff, nro, nphi)

    reltrans.set_re   (nro,nphi)
    reltrans.set_taudo(nro,nphi)
    reltrans.set_pem  (nro,nphi)
    reltrans.grtrace(
        nphi=nphi,
        rn=rn,
        mueff=mueff,
        mu0=mu0,
        spin=spin,
        rmin=rin,
        rout=rout,
        mudisk=mudisk,
        d=distance,
    )
    re    = reltrans.get_re   (nro,nphi)  
    taudo = reltrans.get_taudo(nro,nphi)  
    pem   = reltrans.get_pem  (nro,nphi)  
    assert_snapshot(re, name="re1", rtol=2e-4)
    assert_snapshot(taudo, name="taudo1", rtol=2e-4)
    assert_snapshot(pem, name="pem1", rtol=2e-4)


# def test_basic_invocation_rtdist(reltrans, assert_snapshot, envars):
#     """A smoke test to check whether the default values are working."""
#     reltrans.reset()
#     energy = np.logspace(np.log10(0.1), np.log10(100), 501)
#     output = reltrans.rtdist(energy, rtdist_Parameters())
#     # _debug_plot(energy,output, "rtdist time-averaged spectrum [default parameters]")
#     assert_snapshot(output)

    
# def test_re_im_rtdist(reltrans, assert_snapshot, telescope, envars):
#     """Test the re_im parameter to assert that all the different outputs of the
#     RTDIST model are working. This test requires an RMF and ARF, which is provided by
#     the `telescope` fixture."""
#     reltrans.reset()
#     energy = np.logspace(np.log10(0.1), np.log10(100), 101)

#     envars["RMF_SET"] = telescope.rmf_path
#     envars["ARF_SET"] = telescope.arf_path
#     envars["EMIN_REF"] = "0.3"
#     envars["EMAX_REF"] = "10.0"

#     agn1 = rtdist_Parameters(mass=4.7e6, flo_hz=1e-5, fhi_hz=1e-4, re_im=4.0)
#     output = reltrans.rtdist(energy, agn1)
#     # _debug_plot(energy,output, title = "rtdist lag spectrum", xlabel="Energy [keV]")
#     assert_snapshot(output, name="time_lag")

#     agn1.re_im = 3
#     output = reltrans.rtdist(energy, agn1)
#     # _debug_plot(energy,output, title="rtdist modulus spectrum", xlabel="Energy[keV]")
#     assert_snapshot(output, name="magnitude")

#     agn1.re_im = 1
#     output = reltrans.rtdist(energy, agn1)
#     # _debug_plot(energy,output, title="rtdist real part spectrum", xlabel="Energy[keV]")
#     assert_snapshot(output, name="real_part")

#     agn1.re_im = 2
#     output = reltrans.rtdist(energy, agn1)
#     # _debug_plot(energy,output, title="rtdist imaginary part spectrum", xlabel="Energy[keV]")
#     assert_snapshot(output, name="imaginary_part")


# def test_basic_invocation_model_iterations_rtdist(reltrans, assert_snapshot):
#     energy = np.logspace(np.log10(0.1), np.log10(100), 501)
#     reltrans.reset()
#     output_rtdist = reltrans.rtdist(energy, DCP_Parameters())
#     reltrans.reset()
#     output_rtdist2 = reltrans.rtdist(energy, DCP_Parameters())
#     np.testing.assert_allclose(output_rtdist, output_rtdist2, rtol=1e-4)


# def test_invocation_model_iterations_rtdist_lag(reltrans, telescope, envars):
#     reltrans.reset()
#     envars["RMF_SET"] = telescope.rmf_path
#     envars["ARF_SET"] = telescope.arf_path
#     envars["EMIN_REF"] = "0.3"
#     envars["EMAX_REF"] = "10.0"
#     energy = np.logspace(np.log10(0.1), np.log10(100), 501)
#     agn1 = rtdist_Parameters(mass=4.7e6, flo_hz=1e-5, fhi_hz=1e-4, re_im=4.0)
#     output_rtdist_lag = reltrans.rtdist(energy, agn1)
#     output_rtdist = reltrans.rtdist(energy, rtdist_Parameters())
#     reltrans.reset()
#     output_rtdist_lag2 = reltrans.rtdist(energy, agn1)
#     np.testing.assert_allclose(output_rtdist_lag, output_rtdist_lag2, rtol=1e-4)
