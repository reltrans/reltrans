import pytest
import numpy as np
from pyreltrans import DCP_Parameters, Dbl_Parameters

# Files written by write_components, one per model component
COMPONENT_FILES = ("PivPL.dat", "Reverb.dat", "PivRef.dat", "IonVar.dat")

# re_im 1-4 rebin the components (crebin), 5-6 fold them through the response
# (cfoldandbin) and negative values skip the response matrix
# (propercross_NOmatrix). reltransDbl with beta_p = 0 takes the zero-coherence
# path (components_nocoh).
CASES = [("dcp", re_im) for re_im in (1, 2, 3, 4, 5, 6, -4)] + [("dbl", 4)]


@pytest.mark.parametrize("flavour, re_im", CASES)
def test_write_components(
    reltrans, assert_snapshot, telescope, envars, tmp_path, monkeypatch,
    flavour, re_im
):
    """With REV_VERB > 1, the lag-energy outputs (|ReIm| = 1..6) also write
    each model component to Output/*.dat through write_components. Check that
    every file is written with one finite row per energy bin, and snapshot the
    component values."""
    # write_components writes to Output/ relative to the working directory
    monkeypatch.chdir(tmp_path)
    (tmp_path / "Output").mkdir()

    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"
    # REV_VERB is only read on the first call after a reset
    envars["REV_VERB"] = "2"
    reltrans.reset()

    energy = np.logspace(np.log10(0.1), np.log10(100), 101)
    try:
        if flavour == "dcp":
            params = DCP_Parameters(
                mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=float(re_im),
                del_a=0.1, del_ab=0.2, g=0.5,
            )
            reltrans.dcp(energy, params)
        else:
            params = Dbl_Parameters(
                mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=float(re_im),
                beta_p=0.0, eta=0.8,
            )
            reltrans.dbl_lamp(energy, params)
    finally:
        # The verbosity is kept in the global config, so reset again to make
        # the next call re-read REV_VERB (restored by `envars`). Otherwise
        # later tests would keep writing to Output/ in their own directory.
        reltrans.reset()

    components = []
    for name in COMPONENT_FILES:
        data = np.loadtxt(tmp_path / "Output" / name)
        assert data.shape == (len(energy) - 1, 2), name
        assert np.all(np.isfinite(data)), name
        # First column is the bin-centre energy
        np.testing.assert_allclose(data[:, 0], 0.5 * (energy[1:] + energy[:-1]))
        components.append(data[:, 1])
    components = np.array(components)

    # Lags and imaginary parts cross zero, so a pure relative tolerance is not
    # meaningful there; scale the absolute tolerance to the data instead.
    assert_snapshot(
        components,
        name=f"{flavour}_reim{re_im}",
        atol=1e-9 * np.max(np.abs(components)),
    )
