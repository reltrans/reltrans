import os
import sys
import pathlib
import subprocess

import pytest

ROOT = pathlib.Path(__file__).parent.parent

COMPONENT_FILES = ("PivPL.dat", "Reverb.dat", "PivRef.dat", "IonVar.dat")

# Run in a subprocess: the Fortran output is only reliably flushed at exit,
# and an error writing the output files aborts the whole process.
SCRIPT = """
import numpy as np
from pyreltrans import Reltrans, DCP_Parameters
energy = np.logspace(np.log10(0.1), np.log10(100), 101)
params = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)
Reltrans().dcp(energy, params)
"""


@pytest.mark.parametrize("level", [0, 1, 2, 3])
def test_verbose_levels(telescope, tmp_path, level):
    """Each REV_VERB level adds its own terminal output, and from level 2 the
    model components are written to Output/, which is created if missing."""
    env = dict(
        os.environ,
        REV_VERB=str(level),
        RMF_SET=telescope.rmf_path,
        ARF_SET=telescope.arf_path,
        EMIN_REF="0.3",
        EMAX_REF="10.0",
        PYTHONPATH=os.pathsep.join(filter(None, [str(ROOT), os.environ.get("PYTHONPATH")])),
    )
    result = subprocess.run(
        [sys.executable, "-c", SCRIPT],
        cwd=tmp_path, env=env, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    expected = {
        "reflection fraction for each source": level >= 1,
        "Peak ionisation from LP": level >= 2,
        "inside components": level >= 2,
        "Convolutions runtime": level >= 3,
    }
    for text, present in expected.items():
        assert (text in result.stdout) == present, text
    for name in COMPONENT_FILES:
        assert (tmp_path / "Output" / name).is_file() == (level >= 2), name
