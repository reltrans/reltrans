import numpy as np
import pytest
from pyreltrans import DCP_Parameters

COMPONENT_FILES = ("PivPL.dat", "Reverb.dat", "PivRef.dat", "IonVar.dat")


@pytest.mark.parametrize("level", [0, 1, 2, 3])
def test_verbose_levels(
    reltrans, telescope, envars, tmp_path, monkeypatch, capfd, level
):
    """Each REV_VERB level adds its own terminal output, and from level 2 the
    model components are written to Output/, which is created if missing."""
    monkeypatch.chdir(tmp_path)
    envars["REV_VERB"] = str(level)
    envars["RMF_SET"] = telescope.rmf_path
    envars["ARF_SET"] = telescope.arf_path
    envars["EMIN_REF"] = "0.3"
    envars["EMAX_REF"] = "10.0"
    # REV_VERB is only read on the first call after a reset
    reltrans.reset()

    energy = np.logspace(np.log10(0.1), np.log10(100), 101)
    params = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)
    try:
        reltrans.dcp(energy, params)
    finally:
        # Make the next call re-read REV_VERB (restored by `envars`)
        reltrans.reset()
    # conftest makes Fortran stdout unbuffered
    output = capfd.readouterr().out

    expected = {
        "reflection fraction for each source": level >= 1,
        "Peak ionisation from LP": level >= 2,
        "inside components": level >= 2,
        "Convolutions runtime": level >= 3,
    }
    for text, present in expected.items():
        assert (text in output) == present, text
    for name in COMPONENT_FILES:
        assert (tmp_path / "Output" / name).is_file() == (level >= 2), name
