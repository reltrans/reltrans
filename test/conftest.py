import os
import sys
import inspect
import pathlib
import logging
import dataclasses

from typing import Any

# The `test` directory
ROOT_DIR = pathlib.Path(pathlib.Path(__file__).parent)
SNAPSHOT_DIR = ROOT_DIR / "_snapshots"

import pytest

try:
    import pyreltrans
except ModuleNotFoundError:
    sys.path.append(str(ROOT_DIR / ".."))
    import pyreltrans

import numpy as np

logger = logging.getLogger(__name__)


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
def reltrans() -> pyreltrans.Reltrans:
    """
    Obtain the reltrans library wrapper class.

    By returning a session scoped fixture, this class is essentially a
    singleton, and the same instance is used by all tests. This avoids having
    to load the reltrans library several times, and allows the reltrans cached
    to be reused between tests (eliding loading the xillver tables over and
    over).

    **Caveat**: this does mean values cached in one reltrans invocation will be
    potentially reachable by other tests.
    """
    # only initialise the library once
    return pyreltrans.Reltrans()


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

    The tolerances have been empirically determined during the refactoring. It
    turns out the order of operations incrues different rounding errors, which,
    as reltrans functions over a very large numerical range (that is, very many
    orders of magnitude), build up to the ~0.1% percent level.
    """

    def _assert_snapshot_equal(data: np.ndarray, name="", rtol=2e-4, atol=0) -> bool:
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

        np.testing.assert_allclose(data, snapshot, rtol=rtol, atol=atol)

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
