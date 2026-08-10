import pytest
import numpy as np
from pyreltrans import DCP_Parameters, Dbl_Parameters, rtdist_Parameters

def test_kerrz_observer_disk(reltrans):
    aspin = 0.998
    cos0  = np.cos(np.deg2rad(30))
    alpha = 3.0 #totally random
    beta  = 4.0 #totally random
    t, r, theta, phi = reltrans.kerrz_trace(aspin, cos0, alpha, beta)
    assert r == pytest.approx(3.3090221511440556, rel=1e-4)
    assert theta == pytest.approx(np.pi / 2.0, rel=1e-4)
    assert t == pytest.approx(18000034.946610235, rel=1e-4)

def test_kerrz_lamppost(reltrans):
    deltas = np.deg2rad(40)
    aspin  = 0.998
    h      = 6.0
    t, r, theta, phi = reltrans.kerrz_trace_lamppost(aspin, h, deltas)
    assert r  == pytest.approx(2.9239091166396736, rel=1e-4)
    assert theta == pytest.approx(np.pi / 2.0, rel=1e-4)
    assert t  == pytest.approx(10.236485162732807, rel=1e-4)

def test_kerrz_lensing(reltrans):
    lens, cosdelta, time = reltrans.kerrz_trace_lensing(0.998, 10.0, 1e5, np.cos(np.deg2rad(45)))
    assert lens == pytest.approx(0.8021903536100636, rel=1e-6)
    assert cosdelta == pytest.approx(-0.765057678774833, rel=1e-6)
    assert time == pytest.approx(100011.62907243619, rel=1e-6)

def test_getlens(reltrans):
    lens, del_t, cosdelta = reltrans.getlens(0.998, 10.0, np.cos(np.deg2rad(45)))
    assert lens == pytest.approx(0.8021903320813445, rel=1e-6)
    assert del_t == pytest.approx(11.629330638883403, rel=1e-6)
    assert cosdelta == pytest.approx(-0.7650576791076782, rel=1e-6)

    lens, del_t, cosdelta = reltrans.getlens(0.998, 10.0, np.cos(np.deg2rad(2)))
    assert lens == pytest.approx(0.8021328267886598, rel=1e-6)
    assert del_t == pytest.approx(8.716921714367345, rel=1e-6)
    assert cosdelta == pytest.approx(-0.9995113580864651, rel=1e-6)

    lens, del_t, cosdelta = reltrans.getlens(0.998, 10.0, np.cos(np.deg2rad(88)))
    assert lens == pytest.approx(0.8047144093640742, rel=1e-6)
    assert del_t == pytest.approx(18.706618177704513, rel=1e-6)
    assert cosdelta == pytest.approx(-0.22529823932104842, rel=1e-6)

def test_trace_observer_disk_single_photon(reltrans):
    '''This test computes the ray tracing from the observer to the disk for a single geodesic'''
    reltrans.reset()
    aspin = 0.998
    cos0  = np.cos(30.0/180.0 * np.pi)
    sin0  = np.sqrt(1.0 - cos0**2)
    dist  = 18000000.0
    scal  = 1.0
    alpha = 3.0 #totally random
    beta  = 4.0 #totally random
    #from the observer camera parameter to the Carter's constants of motion
    four_momentum, lambda_, q = reltrans.constants_of_motion(-alpha,-beta,dist,sin0,cos0,aspin,scal)
    mudisk  = 0.0
    r_max   = 1e8
    r_min   = 0.0
    #from the Carter's constant of motion to the affine parameter where the geodesic hit the disk
    p_out = reltrans.p_disk_crossing(four_momentum,lambda_.value,q.value,sin0,cos0,aspin,dist,scal,mudisk,r_max,r_min)
    #from the affine parameter and constant of motion to the interesting values
    radi, mu, phi, time, sigma = reltrans.get_raytrace_coords(p_out,four_momentum,lambda_,q,sin0,cos0,aspin,dist,scal)
    # print(f'FROM THE TESTS: radi {radi}, mu {mu}, phi {phi}, time {time}, sigma {sigma}')
    assert radi.value  == pytest.approx(3.3090221511440556, rel=1e-4)
    assert mu.value    == pytest.approx(0.0, rel=1e-4)
    assert time.value  == pytest.approx(18000034.946610235, rel=1e-4)
    assert sigma.value == pytest.approx(18000000.171032075, rel=1e-4)


def test_trace_source_disk_single_photon(reltrans):
    '''This test computes the ray-tracing from the lamppost source to the disk for a single geodesic'''
    reltrans.reset()
    deltas = 40.0/180.0 * np.pi #180 degree out from kerrz
    pr     = np.cos(deltas)           #cosdelta
    pp     = np.sqrt( 1.0 - pr**2 )   #sindelta
    pt     = 0.0
    mus    = 1.0   #Position of source: mus=1 means on-axis
    sins   = np.sqrt(1.0 - mus**2)   #sin of same angle
    aspin  = 0.998
    h      = 6.0
    scal   = 1.0
    #from the source paramter to the Carter's constants of motion
    four_momentum, lambda_, q = reltrans.initial_photon(pr,pt,pp,sins,mus,aspin,h)
    mudisk  = 0.0
    r_max   = 300.0
    r_min   = 1.3
    #from the Carter's constant of motion to the affine parameter where the geodedic hit the disk
    p_out = reltrans.p_disk_crossing(four_momentum,lambda_,q,sins,mus,aspin,h,
          scal,mudisk,r_max,r_min)
    #from the affine parameter and constant of motion to the interesting values
    radi, mu, phi, time, sigma = reltrans.get_raytrace_coords(p_out,four_momentum,lambda_.value,q.value,sins,mus,aspin,h,scal)
    # print(f'FROM THE TESTS: radi {radi}, mu {mu}, phi {phi}, time {time}, sigma {sigma}')
    assert radi.value  == pytest.approx(2.9239091166396736, rel=1e-4)
    assert mu.value    == pytest.approx(0.0, rel=1e-4)
    assert time.value  == pytest.approx(10.567334777173707, rel=1e-4)
    assert sigma.value == pytest.approx(5.442883301224947, rel=1e-4)
