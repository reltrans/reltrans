import os
import inspect
import pathlib
import logging
import dataclasses

from typing import Any

import pytest
import wrapper

import numpy as np

logger = logging.getLogger(__name__)

# The `test` directory
ROOT_DIR = pathlib.Path(pathlib.Path(__file__).parent)
SNAPSHOT_DIR = ROOT_DIR / "_snapshots"


def _get_snapshot(name: str) -> None | np.ndarray:
    """Retrieve the snapshot by name or None if the file does not exist."""
    path = SNAPSHOT_DIR / name
    path = path.with_suffix("".join(path.suffixes) + ".npy")
    if not path.is_file():
        return None
    return np.load(str(path.absolute()))


def _create_snapshot(name: str, data: np.array):
    """Create a snapshot by name and store the data in the given numpy array."""
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    path = SNAPSHOT_DIR / name
    np.save(str(path.absolute()), data)


class EnvironmentVariables:
    def __init__(self):
        self.original_vars = {}

    def __setitem__(self, name: str, value: str):
        """Set an environment variable."""
        if name not in self.original_vars:
            self.original_vars[name] = os.environ.get(name, None)
        os.environ[name] = value

    def __getitem__(self, name: str) -> Any:
        """Get an environment variable."""
        return os.environ[name]

    def get(self, name: str, default=None) -> Any:
        """Get an environment variable with `default = None`."""
        return os.environ.get(name, default)

    def _restore(self):
        """Restore the original environment variables."""
        for name, v in self.original_vars.items():
            if v is None:
                del os.environ[name]
            else:
                os.environ[name] = v


@pytest.fixture(scope="function")
def envars() -> EnvironmentVariables:
    """
    Used to set environment variables that are only valid for the duration
    of a particular test.

    These may be used as a greatly simplified version of `os.environ`:

        def test_mytest(envars):
            envars["SOMETHING"] = "VALUE"
            assert envars["SOMETHING"] == "VALUE"
            assert envars.get("DIFFERENT", "5") == "5"

    All environment variables are unset or restored to their original values
    from before the test ran.

    Note this is not thread safe.
    """
    ev = EnvironmentVariables()
    yield ev
    ev._restore()

@pytest.fixture(scope="session")
def reltrans() -> wrapper.Reltrans:
    """
    Session-scoped fixture returning a single Reltrans instance.

    If the shared library is not available, tests are skipped gracefully
    with a clear diagnostic message.
    """
    import os
    import pytest
    import pathlib
    import platform

    if os.environ.get("RELTRANS_DEV"):
        pytest.skip("Skipping native library tests in dev mode")

    system = platform.system()
    lib_name = "libreltrans.so" if system == "Linux" else "libreltrans.dylib"

    project_root = pathlib.Path(__file__).resolve().parents[1]
    expected_path = project_root / "build" / "lib" / lib_name

    if not expected_path.exists():
        pytest.skip(
            "\nReltrans shared library not found.\n"
            f"Expected location: {expected_path}\n"
            f"HEADAS: {os.environ.get('HEADAS')}\n\n"
            "To enable native tests:\n"
            "  1. Ensure HEASOFT is installed and sourced.\n"
            "  2. Run `make` in the reltrans root directory.\n"
        )

    try:
        return wrapper.Reltrans()
    except OSError as e:
        pytest.skip(
            "\nReltrans library failed to load.\n"
            f"Error: {e}\n"
            f"HEADAS: {os.environ.get('HEADAS')}\n"
        )



@pytest.fixture
def assert_snapshot() -> callable:
    """
    Fixture used to assert that snapshots are reproduced to within specified
    tolerances.

    Returns a function which raises an assertion error if the sole argument
    provided to it does not match a cached snapshot in the
    `$ROOT_DIR/_snapshots` directory.

    If a snapshot does not exist, it creates a new file. Names are derived from
    the calling function (i.e. the test case) and may optionally have a `name`
    argument appended so that `assert_snapshot()` may be used multiple times
    with different data by the same test.
    """

    def _assert_snapshot_equal(data: np.ndarray, name="", rtol=1e-4) -> bool:
        calling_context = inspect.stack()[1][3]

        snapshot_name = calling_context
        if name:
            snapshot_name = f"{snapshot_name}.{name}"

        snapshot = _get_snapshot(snapshot_name)
        if snapshot is None:
            # TODO: print warning that no snapshot exists and that a new one
            # has been created
            logger.warning(
                "Snapshot %s does not exist. Creating new from supplied data",
                snapshot_name,
            )
            _create_snapshot(snapshot_name, data)
            return True

        np.testing.assert_allclose(data, snapshot, rtol=rtol)

    return _assert_snapshot_equal


@dataclasses.dataclass
class TelescopeData:
    arf_path: str
    rmf_path: str


@pytest.fixture
def telescope() -> TelescopeData:
    """Obtain paths to telescope files."""
    repo_root = pathlib.Path(ROOT_DIR.parent)
    rmf_path = (
        repo_root / "Benchmarks" / "resp_matrix" / "nicer-rmf6s-teamonly-array50.rmf"
    )
    arf_path = (
        repo_root
        / "Benchmarks"
        / "resp_matrix"
        / "nicer-consim135p-teamonly-array50.arf"
    )
    return TelescopeData(
        arf_path=str(arf_path.absolute()), rmf_path=str(rmf_path.absolute())
    )
