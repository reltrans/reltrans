import pytest
import numpy as np
from pyreltrans import DCP_Parameters, Dbl_Parameters, rtdist_Parameters


def _debug_plot(
        energy, output, title="", xlabel="", ylabel="", yscale="linear", xscale="log", ylim_min=None, ylim_max=None):
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
    if ylim_min is not None and ylim_max is not None:
        plt.ylim(ylim_min, ylim_max)
    plt.show()

    
def _debug_plot_multi_spec(energy, output, title="", xlabel="", ylabel="", yscale="linear", xscale="log"):
    """Used when creating new tests for quickly looking at the data to make
    sure it is sensible."""
    import matplotlib.pyplot as plt

    plt.clf()
    for spec in output:
        plt.plot(energy[0:-1], spec)
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
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(nh = 0.2))
    # _debug_plot(energy,output, "reltransDCp time-averaged spectrum [default parameters]")
    assert_snapshot(output)


def test_dcp_continuum(reltrans, assert_snapshot):
    """A test to check the output of the continum  (boost = 0)"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(boost = 0))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot(energy, output/dE*E**2, "reltransDCp time-averaged spectrum [boost = 0]", yscale="log", xscale="log")
    assert_snapshot(output)


def test_dcp_reflection(reltrans, assert_snapshot):
    """A test to check the output of boost = -1 (just reflection)"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters(boost = -1))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot(energy,output/dE*E**2, "reltransDCp time-averaged spectrum [boost = -1]", yscale="log", xscale="log")
    assert_snapshot(output)


def test_dcp_to_dbl_continuum(reltrans, assert_snapshot):
    """Test to compare the continuum output of reltransDCp and reltransDBL (boost = 0)"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = 0)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = 0))
    # ratio = output_dbl/output_dcp
    # outputs = np.stack((output_dcp, output_dbl))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot_multi_spec(energy, outputs, "reltrans time-averaged spectrum [boost = 0]", yscale="log", xscale="log")
    # _debug_plot(energy, ratio, "ratio dcp to dbl [boost = -1]", yscale="log", xscale="log")
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_reflection_only(reltrans, assert_snapshot):
    """Test to check if the reflection-only outputs (boost = -1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are the same"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = -1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = xrb.boost))
    # ratio = output_dbl/output_dcp
    # outputs = np.stack((output_dcp, output_dbl))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot_multi_spec(energy, outputs, "reltrans time-averaged spectrum [boost = -1]", yscale="log", xscale="log")
    # _debug_plot(energy, ratio, "ratio dcp to dbl [boost = 0]", yscale="log", xscale="log")
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_eta0_reflection_only(reltrans, assert_snapshot):
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
    # ratio = output_dbl/output_dcp
    # outputs = np.stack((output_dcp, output_dbl))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot_multi_spec(energy, outputs, "reltrans time-averaged spectrum [boost = -1]", yscale="log", xscale="log")
    # _debug_plot(energy, ratio, "ratio dcp to dbl [boost = 0]", yscale="log", xscale="log")
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_full_spectrum(reltrans, assert_snapshot):
    """Test to check if the full time-averaged spectra (boost = 1)
    of reltransDCp and reltransDBL are the same if the heights
    of the two lampposts are the same"""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    xrb = DCP_Parameters(h = 5.0, boost = 1)
    output_dcp = reltrans.dcp(energy, xrb)
    reltrans.reset()
    output_dbl = reltrans.dbl_lamp(energy, Dbl_Parameters(h1 = xrb.h, h2 = xrb.h, a = xrb.a, boost = xrb.boost))
    # ratio = output_dbl/output_dcp
    # outputs = np.stack((output_dcp, output_dbl))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot_multi_spec(energy, outputs, "reltrans time-averaged spectrum [boost = 1]", yscale="log", xscale="log")
    # _debug_plot(energy, ratio, "ratio dcp to dbl [boost = 1]", yscale="log", xscale="log")
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)


def test_dcp_to_dbl_eta0_full_spectrum(reltrans, assert_snapshot):
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
    # ratio = output_dbl/output_dcp
    # outputs = np.stack((output_dcp, output_dbl))
    # E = (energy[1:] + energy[:-1]) * 0.5
    # dE = (energy[1:] - energy[:-1])
    # _debug_plot_multi_spec(energy, outputs, "reltrans time-averaged spectrum [boost = 1]", yscale="log", xscale="log")
    # _debug_plot(energy, ratio, "ratio dcp to dbl [boost = 1]", yscale="log", xscale="log")
    np.testing.assert_allclose(output_dcp, output_dbl, rtol=1e-4)
    

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
    assert_snapshot(output, name="time_lag", atol=1e-9)

    xrb1.re_im = 1
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="real_part")

    xrb1.re_im = 2
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="imaginary_part", rtol=1e-3)

    xrb1.re_im = 3
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="magnitude")


def test_re_im_5_6(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(12), 101)
    dE = (energy[1:] - energy[:-1])

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    xrb1 = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=6.0)
    output = reltrans.dcp(energy, xrb1)/dE
    # _debug_plot(energy,output, title = "reltrans lag spectrum", xlabel="Energy [keV]", ylim_min=-1e-3, ylim_max=1e-3)
    assert_snapshot(output, name="time_lag")

    xrb1.re_im = 5
    output = reltrans.dcp(energy, xrb1)
    # _debug_plot(energy,output, title = "reltrans modulus spectrum", xlabel="Energy [keV]")
    assert_snapshot(output, name="real_part")


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
    assert_snapshot(output, name="time_lag", atol=1e-9)

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
    assert_snapshot(output, name="imaginary_part", rtol=1e-3)


def test_basic_invocation_rtdist(reltrans, assert_snapshot, envars):
    """A test to check whether the default values of rtdist model are working
    with density profile as SSD73."""
    reltrans.reset()
    envars["A_DENSITY"] = "1"
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.rtdist(energy, rtdist_Parameters())
    # _debug_plot(energy,output, "rtdist time-averaged spectrum [default parameters]")
    assert_snapshot(output)


def test_basic_invocation_rtdist_adens0(reltrans, assert_snapshot, envars):
    """A test to check whether the default values of rtdist are working
    with density profile constant."""
    reltrans.reset()
    envars["A_DENSITY"] = "0"
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.rtdist(energy, rtdist_Parameters())
    # _debug_plot(energy,output, "rtdist time-averaged spectrum [default parameters]")
    assert_snapshot(output)


def test_re_im_rtdist(reltrans, assert_snapshot, telescope, envars):
    """Test the re_im parameter to assert that all the different outputs of the
    RTDIST model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    reltrans.reset()
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    envars["A_DENSITY"] = "1"
    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"

    agn1 = rtdist_Parameters(mass=4.7e6, flo_hz=1e-5, fhi_hz=1e-4, re_im=4.0)
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title = "rtdist lag spectrum", xlabel="Energy [keV]")
    assert_snapshot(output, name="time_lag")

    agn1.re_im = 3
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist modulus spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="magnitude")

    agn1.re_im = 1
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist real part spectrum", xlabel="Energy[keV]")
    assert_snapshot(output, name="real_part")

    agn1.re_im = 2
    output = reltrans.rtdist(energy, agn1)
    # _debug_plot(energy,output, title="rtdist imaginary part spectrum", xlabel="Energy[keV]")
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
    output_dcp = reltrans.dcp(energy, DCP_Parameters(nh=0.5))
    reltrans.reset()
    output_dcp2 = reltrans.dcp(energy, DCP_Parameters(nh=0.5))
    np.testing.assert_allclose(output_dcp, output_dcp2, rtol=1e-4)


def test_rtdist_called_twice(reltrans, assert_snapshot):
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


def test_strans_routines_getrgrid(reltrans, assert_snapshot):
    rnmin = 1.7266670949940284
    rnmax = 300.0
    mueff = 0.86602540378443871
    nro = 200
    nphi = 200
    rn, domega = reltrans.getrgrid(rnmin, rnmax, mueff, nro, nphi)
    assert_snapshot(rn, name="rn")
    assert_snapshot(domega, name="domega")


def test_strans_routines_grtrace_outputs(reltrans, assert_snapshot):
    reltrans.reset()

    rnmin = 1.7266670949940284
    rnmax = 300.0
    mueff = 0.86602540378443871
    nro = 200
    nphi = 200
    mu0 = mueff
    spin = 0.998
    rin = 1.2369694567256766
    rout = 1000.0
    mudisk = 0.0
    distance = 18000000.0

    rn, domega = reltrans.getrgrid(rnmin, rnmax, mueff, nro, nphi)

    reltrans.set_re(nro, nphi)
    reltrans.set_taudo(nro, nphi)
    reltrans.set_pem(nro, nphi)
    reltrans.trace_disk_observer(
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
    re = reltrans.get_re(nro, nphi)
    taudo = reltrans.get_taudo(nro, nphi)
    pem = reltrans.get_pem(nro, nphi)
    assert_snapshot(re, name="re1", rtol=2e-4)
    assert_snapshot(taudo, name="taudo1", rtol=2e-4)
    assert_snapshot(pem, name="pem1", rtol=2e-4)


def test_raytrace_disk(reltrans, assert_snapshot):
    reltrans.reset()
    p = 0.16638448413345461
    f1234 = [0.97999999999999998, 0, -0.1989974874213242, 1]
    q = 1.1723701732639298
    lambda_ = 0.0
    sinobs = np.sin(1e-3)
    muobs = np.cos(1e-3)
    aspin = 0.998
    robs  = 6.0
    scal = 1.0
    radi, mu, phi, time, sigma = reltrans.raytrace_disk(p,f1234,lambda_,q,sinobs,muobs,aspin,robs, scal)
    # print(f'FROM THE TESTS: radi {radi}, mu {mu}, phi {phi}, time {time}, sigma {sigma}')
    assert radi.value  == pytest.approx(60095.57339645814, rel=1e-4) 
    assert mu.value    == pytest.approx(0.9701656408593914, rel=1e-4) 
    assert phi.value   == pytest.approx(6.247995161686957, rel=1e-6) 
    assert time.value  == pytest.approx(60108.777370445634, rel=1e-6) 
    assert sigma.value == pytest.approx(60089.55895566616, rel=1e-6)


def test_p_disk_crossing(reltrans, assert_snapshot):
    reltrans.reset()
    f1234 = [0.97999999999999998, 0, -0.1989974874213242, 1]
    q = 1.1723701732639298
    lambda_ = 0.0
    sins  = np.sin(10)
    mus   = np.cos(10)
    aspin   = 0.998
    h       = 6.0
    scal    = 1.0
    mudisk  = np.cos(0)
    r_max   = 300.0
    r_min   = 1.3
    
    p_out = reltrans.p_disk_crossing(f1234,lambda_,q,sins,mus,aspin,h,
          scal,mudisk,r_max,r_min)
    print(f'FROM THE TESTS: p_out = {p_out}')
    # p_ref = 0.16638448413345461
    assert p_out.value == -1


def test_p_coord_at_infinity(reltrans, assert_snapshot):
    reltrans.reset()
    f1234 = [0.97999999999999998, 0, -0.1989974874213242, 1]
    q = 1.1723701732639298
    lambda_ = 0.0
    sins  = np.sin(0)
    mus   = np.cos(0)
    aspin   = 0.998
    h       = 6.0
    scal    = 1.0
    p_out = reltrans.p_coord_at_infinity(f1234,lambda_,q,sins,mus,aspin,h,scal)
    print(f'FROM THE TESTS: p_out = {p_out}')
    assert 0 == -1

    
def test_ReIm8_check_second_response(reltrans,  assert_snapshot, telescope, envars):

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
