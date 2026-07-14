import pytest
import numpy as np
from pyreltrans import DCP_Parameters, Dbl_Parameters, rtdist_Parameters


plot_spectral_kwargs = dict(
        yscale = "log", xlabel = "Energy", ylabel = "Flux" , xscale = "log", units = "ef"
)

def test_dcp_to_dbl_continuum(reltrans, save_plot):
    """Test to compare the continuum output of reltransDCp and reltransDBL (boost = 0)"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = 0)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = 0))
    save_plot(
        energy[0:-1],
        output_dcp,
        output_dbl,
        label1="DCP",
        label2="DBL",
        rtol = 1e-4,
        **plot_spectral_kwargs,
    )
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_reflection_only(reltrans, save_plot):
    """Test to check if the reflection-only outputs (boost = -1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are the same"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = -1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = xrb.boost))
    save_plot(
        energy[0:-1],
        output_dcp,
        output_dbl,
        label1="DCP",
        label2="DBL",
        rtol = 1e-4,
        **plot_spectral_kwargs,
    )
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_eta0_reflection_only(reltrans, save_plot):
    """Test to check if the reflection-only outputs (boost = -1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are different, BUT eta_0 == 0 (which means that
    the second lamppost does NOT contribute)."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = -1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = 100.0, a = xrb.a, boost = xrb.boost, eta_0 = 0.0))
    save_plot(
        energy[0:-1],
        output_dcp,
        output_dbl,
        label1="DCP",
        label2="DBL",
        rtol = 1e-4,
        **plot_spectral_kwargs,
    )
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_full_spectrum(reltrans, save_plot):
    """Test to check if the full time-averaged spectra (boost = 1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are the same"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = 1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = xrb.boost))
    save_plot(
        energy[0:-1],
        output_dcp,
        output_dbl,
        label1="DCP",
        label2="DBL",
        rtol = 1e-4,
        **plot_spectral_kwargs,
    )
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_eta0_full_spectrum(reltrans, save_plot):
    """Test to check if the full time-averaged spectra (boost = 1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are different, BUT eta_0 == 0 (which means that
    the second lamppost does NOT contribute)."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = 1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = xrb.boost, eta_0 = 0))
    save_plot(
        energy[0:-1],
        output_dcp,
        output_dbl,
        label1="DCP",
        label2="DBL",
        rtol = 1e-4,
        **plot_spectral_kwargs,
    )
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_ReIm8_check_second_response(reltrans, telescope, envars):
    '''A test for checking if the second response is loaded when ReIm=8'''
    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"
    envars["RMF2SET"] = telescope.rmf_path
    envars["ARF2SET"] = telescope.arf_path

    energy = np.logspace(np.log10(0.1), np.log10(100), 101)
    reltrans.reset()
    xrb1 = DCP_Parameters(mass=10.0, flo_hz=1, fhi_hz=2, re_im=6.0)
    output = reltrans.dcp(energy, xrb1)
    resp2_needed = reltrans.get_needresp2()
    assert resp2_needed

    xrb1 = DCP_Parameters(mass=10.0, flo_hz=1, fhi_hz=2, re_im=8.0)
    output = reltrans.dcp(energy, xrb1)
    resp2_needed = reltrans.get_needresp2()
    assert not resp2_needed


def test_reltransDCp_called_twice_no_resetting(reltrans):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_dcp = reltrans.dcp(energy, DCP_Parameters())
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters())
    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)


def test_reltransDCp_called_twice_after_resetting(reltrans, assert_snapshot):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_dcp = reltrans.dcp(energy, DCP_Parameters(nh=0.5))
    reltrans.reset()
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters(nh=0.5))
    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)


def test_rtdist_called_twice(reltrans):
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    reltrans.reset()
    output_rtdist = reltrans.rtdist(energy, rtdist_Parameters())
    output_rtdist2 = reltrans.rtdist(energy, rtdist_Parameters())
    np.testing.assert_allclose(output_rtdist, output_rtdist2, rtol=1e-4)


def test_rtdist_called_twice_lag_and_time_averaged(reltrans, telescope, envars):
    reltrans.reset()
    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    agn1 = rtdist_Parameters(mass=4.7e6, flo_hz=1e-5, fhi_hz=1e-4, re_im=4.0)
    output_rtdist_lag = reltrans.rtdist(energy, agn1)
    output_rtdist = reltrans.rtdist(energy, rtdist_Parameters()) #call the time averaged spectrum
    output_rtdist_lag2 = reltrans.rtdist(energy, agn1)
    output_rtdist2 = reltrans.rtdist(energy, rtdist_Parameters())
    np.testing.assert_allclose(output_rtdist_lag, output_rtdist_lag2, rtol=1e-4)
    np.testing.assert_allclose(output_rtdist, output_rtdist2, rtol=1e-4)


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


def test_negative_re_im_called_twice(reltrans, telescope, envars):
    """Test """
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    envars["EMIN_REF"] = "2.0"
    envars["EMAX_REF"] = "10.0"

    xrb1 = DCP_Parameters(mass=10.0, flo_hz=0.1, fhi_hz=0.2, re_im=-4.0)
    output_propercrossNOmatrix = reltrans.dcp(energy, xrb1)
    reltrans.reset()
    xrb1 = DCP_Parameters(mass=10.0, flo_hz=1, fhi_hz=5, re_im=-3.0)
    output_propercross = reltrans.dcp(energy, xrb1)
    reltrans.reset()
    xrb1 = DCP_Parameters(mass=10.0, flo_hz=0.1, fhi_hz=0.2, re_im=-4.0)
    output_propercrossNOmatrix2 = reltrans.dcp(energy, xrb1)

    np.testing.assert_allclose(output_propercrossNOmatrix, output_propercrossNOmatrix2, rtol=1e-4)
