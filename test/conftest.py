import os
import sys
import inspect
import pathlib
import logging
import dataclasses
import enum

from typing import Any

DEFAULT_TOLERANCE_ABS = 0
DEFAULT_TOLERANCE_REL = 2e-4

# The `test` directory
ROOT_DIR = pathlib.Path(pathlib.Path(__file__).parent)
SNAPSHOT_DIR = ROOT_DIR / "_snapshots"
PLOT_DIR = ROOT_DIR / "_plots"

import pytest

try:
    import pyreltrans
except ModuleNotFoundError:
    sys.path.append(str(ROOT_DIR / ".."))
    import pyreltrans

import numpy as np

logger = logging.getLogger(__name__)


def _get_calling_function_name(name: str | None) -> str:
    """
    Get the calling function's name. A name may be provided that will be
    appended as a tail.
    """
    caller_name = inspect.currentframe().f_back.f_back.f_code.co_name
    if name:
        return f"{caller_name}.{name}"
    return caller_name


class PlotMode(enum.Enum):
    ALWAYS_PLOT = 1
    PLOT_ON_FAIL = 2
    NO_PLOT = 0


def pytest_addoption(parser):
    parser.addoption(
        "--plot",
        action="store",
        help="Whether to save debugging plots when tests fail.",
        default=True,
    )
    parser.addoption(
        "--plot-always",
        action="store",
        help="Whether to save debugging plots irrespective of whether tests fail.",
        default=False,
    )


@pytest.fixture
def plotting_mode(request):
    """
    This fixture can be used to get boolean that tells the test or fixture
    whether to enable debug plotting on fail.
    """
    if request.config.getoption("--plot-always"):
        return PlotMode.ALWAYS_PLOT
    if request.config.getoption("--plot"):
        return PlotMode.PLOT_ON_FAIL
    return PlotMode.NO_PLOT


def _plot_with_ratio(
    x: np.ndarray, y1: np.ndarray, y2: np.ndarray, title: str, **kwargs
):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(nrows=2, figsize=(10, 6))

    residuals = (y1 - y2) / y2
    axes[1].plot(x, residuals)

    units = kwargs.get("units", "f")
    if len(units) > 1:
        units_unknown = False

        y1_tmp = y1.copy()
        y2_tmp = y2.copy()

        for c in units[:-1]:
            if c == "e":
                y1_tmp *= x
                y2_tmp *= x
            else:
                units_unknown = True
                break

        if units[-1] != "f":
            units_unknown = True

        if units_unknown:
            logger.warn(f"Unknown plotting units '{units}'. Plotting as is.")
            y1_tmp = y1
            y2_tmp = y2

        y1 = y1_tmp
        y2 = y2_tmp

    axes[0].plot(x, y1, label=kwargs.get("label1", "actual"))
    axes[0].plot(x, y2, label=kwargs.get("label2", "expected"))

    axes[0].set_title(title)
    axes[0].legend()

    # The y-scale is separated for the two axes:
    axes[0].set_yscale(kwargs.get("yscale", "linear"))
    axes[1].set_yscale(kwargs.get("yscale_ratio", "linear"))
    axes[0].set_ylabel(kwargs.get("ylabel", "y") + f" (units: {units})")
    axes[1].set_ylabel("(y1 - y2) / y2")

    print(kwargs)

    for ax in axes:
        ax.set_xlabel(kwargs.get("xlabel", "x"))
        ax.set_xscale(kwargs.get("xscale", "linear"))

    fig.tight_layout()

    return fig


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
def save_plot(plotting_mode) -> callable:
    """
    This fixture can be used to save debugging plots depending on the `--plot`
    command line argument.
    """

    def plot_function(
        domain: np.ndarray | None,
        data: np.ndarray,
        snapshot: np.ndarray,
        name: str | None = None,
        derive_name_from_callstack=True,
        plotter_function: callable = _plot_with_ratio,
        atol=DEFAULT_TOLERANCE_ABS,
        rtol=DEFAULT_TOLERANCE_REL,
        **kwargs,
    ):
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            logger.warn("matplotlib missing, cannot create debug plot.")
            return

        if derive_name_from_callstack:
            name = _get_calling_function_name(name)
        else:
            name = name or "[name missing]"

        if (plotting_mode == PlotMode.ALWAYS_PLOT) or (
            (plotting_mode == PlotMode.PLOT_ON_FAIL)
            and not np.allclose(data, snapshot, rtol=rtol, atol=atol)
        ):
            domain = domain if domain is not None else np.array(range(len(data)))
            fig = plotter_function(domain, data, snapshot, name, **kwargs)
            # Then save the plot after ensuring the plot directory exists:
            PLOT_DIR.mkdir(parents=True, exist_ok=True)
            fig.savefig(PLOT_DIR / (name + ".jpg"))

    return plot_function


@pytest.fixture
def assert_snapshot(save_plot) -> callable:
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

    All other keyword arguments are passed to the debug plotting function.
    """

    def _assert_snapshot_equal(
        data: np.ndarray,
        name="",
        rtol=DEFAULT_TOLERANCE_REL,
        atol=DEFAULT_TOLERANCE_ABS,
        domain: np.ndarray | None = None,
        **kwargs,
    ) -> bool:
        snapshot_name = _get_calling_function_name(name)
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

        save_plot(
            domain,
            data,
            snapshot,
            name=snapshot_name,
            derive_name_from_callstack=False,
            **kwargs,
        )
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
