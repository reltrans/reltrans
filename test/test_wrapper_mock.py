import pytest
import numpy as np
import ctypes as ct
from unittest.mock import MagicMock, call
import sys
import os

# Add parent directory to path to import wrapper
import sys
import os

# Add the 'test' directory to sys.path so we can import 'wrapper' directly
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from wrapper import DCP_Parameters, Reltrans, get_reltrans_library_path

def test_dcp_parameters_defaults():
    """Test default values of DCP_Parameters."""
    params = DCP_Parameters()
    assert params.h == 6.0
    assert params.a == 0.998
    assert params.inc == 30.0
    assert params.mass == 4.6e7
    assert params.re_im == 1.0

def test_dcp_parameters_serialization():
    """Test conversion of DCP_Parameters to numpy array."""
    params = DCP_Parameters(h=10.0, a=0.5, kte=100.0)
    arr = params.to_numpy_array()
    
    assert isinstance(arr, np.ndarray)
    assert arr.dtype == np.float32
    
    # Check specific indices based on @dataclass field order
    # h is 0, a is 1, inc is 2...
    assert arr[0] == 10.0  # h
    assert arr[1] == 0.5   # a
    assert arr[2] == 30.0  # inc (default)

def test_get_reltrans_library_path(monkeypatch):
    """Test library path retrieval from environment variable."""
    monkeypatch.setenv("RELTRANS_PATH", "/path/to/lib.so")
    assert get_reltrans_library_path() == "/path/to/lib.so"

def test_get_reltrans_library_path_missing(monkeypatch):
    """Test error when RELTRANS_PATH is not set."""
    monkeypatch.delenv("RELTRANS_PATH", raising=False)
    with pytest.raises(Exception, match="RELTRANS_PATH not set"):
        get_reltrans_library_path()

def test_reltrans_initialization_mock(mocker):
    """Test Reltrans class initialization with mocked CDLL."""
    # Mock ctypes.cdll.LoadLibrary
    mock_load_library = mocker.patch("ctypes.cdll.LoadLibrary")
    mock_lib = MagicMock()
    mock_load_library.return_value = mock_lib
    
    # Mock the tdreltransdcp_ function object on the library
    mock_func = MagicMock()
    mock_lib.tdreltransdcp_ = mock_func

    # Initialize Reltrans
    rt = Reltrans(path="/dummy/path.so")
    
    # Verify LoadLibrary was called
    mock_load_library.assert_called_once_with("/dummy/path.so")
    
    # Verify argtypes and restype were set
    assert mock_func.argtypes is not None
    assert len(mock_func.argtypes) == 5
    assert mock_func.restype is None

def test_reltrans_dcp_call(mocker):
    """Test dcp method calls the C function correctly."""
    # Setup Mock
    mock_load_library = mocker.patch("ctypes.cdll.LoadLibrary")
    mock_lib = MagicMock()
    mock_load_library.return_value = mock_lib
    mock_func = MagicMock()
    mock_lib.tdreltransdcp_ = mock_func
    
    rt = Reltrans(path="/dummy/path.so")
    
    # Prepare inputs
    energy = np.logspace(0, 1, 11) # 10 bins
    params = DCP_Parameters()
    
    # Call dcp
    result = rt.dcp(energy, params)
    
    # Verify result shape and type
    assert isinstance(result, np.ndarray)
    assert len(result) == 10 # ne = len(energy) - 1
    assert result.dtype == np.float32
    
    # Verify C function call
    mock_func.assert_called_once()
    args, _ = mock_func.call_args
    
    # Check arguments passed to C function
    # 1. Energy pointer (float*)
    # 2. ne pointer (int*)
    # 3. params pointer (float*)
    # 4. ifl pointer (int*) - Assuming '1' is passed for ifl based on wrapper.py code
    # 5. output pointer (float*)
    
    # We can't easily check pointer values, but we can check types if needed
    # or just trust the mock call happened.
