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
    assert t == pytest.approx(100024.22630105753, rel=1e-4)

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

