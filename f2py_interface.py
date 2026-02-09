'''
Simple interface example

This script provides an interface to a fortran library
which you can import from anywhere in your system
'''
import ctypes as ct
import os.path
import numpy as np
from typing import Callable, Any

#######################################################################
# prepare a few pointer types for fortran
#
# these we'll need to deal with all we pass to and from fortran:
# everything is a pointer for fortran
#######################################################################

type_double_p = ct.POINTER(ct.c_double)
type_float_p = ct.POINTER(ct.c_float)
type_int_p    = ct.POINTER(ct.c_int)


#######################################################################
# connect to the library!
#
# using __FILE__  you can put this in a generic place
# the you can import this module from anywhere doing
# import sys
# sys.path.append(path_to_location_of_this_file)
#######################################################################

lib = ct.cdll.LoadLibrary(os.path.dirname(__file__) + "/lib_reltrans.so")

#######################################################################
# now the function(s)
#
# to see names of objects in library use: nm -gDC name_of_lib.so
#######################################################################

# this line is not strictly speaking needed if you are carefull in the call, but nice to do

wPL = lib.tdreltranspl_
wPL.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wPL.restype  = None

wDCp = lib.tdreltransdcp_
wDCp.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
# can actually get one return, even a pointer (a newly allocated array by the function)
wDCp.restype  = None

wx = lib.tdreltransx_
wx.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wx.restype  = None

wDbl = lib.tdreltransdbl_
wDbl.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wDbl.restype  = None

wdist = lib.tdrtdist_
wdist.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wdist.restype  = None

wsim_dist = lib.simrtdist_
wsim_dist.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wsim_dist.restype  = None

wsim_relt = lib.simrelt_
wsim_relt.argtypes = [type_float_p, type_int_p, type_float_p, type_int_p, type_float_p]
wsim_relt.restype  = None

def gen_wrap(ear: np.ndarray, params: np.ndarray, func: Callable[..., None]) -> np.ndarray:
    '''
    Generic wrapper for reltrans Fortran subroutines.

    Args:
        ear (np.ndarray): Array of energy bin edges (keV).
        params (np.ndarray): Array of model parameters (float32).
        func (Callable): The ctypes function object to call.

    Returns:
        np.ndarray: Calculated photon flux array (float32).
    '''

    # to be extra sure you could put the following
    # but it could slow down the code
    #
    # ear    = numpy.array(ear)
    # params = numpy.array(params)

    ne = len(ear) - 1

    photar = np.zeros(ne, dtype = np.float32)

    func(ear.ctypes.data_as(type_float_p),
               ct.byref(ct.c_int(ne)),
               params.ctypes.data_as(type_float_p),
               ct.byref(ct.c_int(1)),
               photar.ctypes.data_as(type_float_p))

    return photar

# def reltrans(ear, params):
#     return gen_wrap(ear, params, w)

def reltransPL(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Calculate transfer function for Power Law illumination.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params, wPL)

def reltransDCp(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Calculate transfer function for Disc-Corona (DCp) geometry.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params, wDCp)

def reltransDbl(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Calculate transfer function for Double Lamp Post geometry.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params,wDbl)

def reltransx(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Calculate transfer function for generic geometry X.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params, wx)

def rtdist(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Calculate Ray Tracing distance effects.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params, wdist)

def simrtdist(ear: np.ndarray, params: np.ndarray) -> np.ndarray:
    '''
    Simulate Ray Tracing distance effects.

    Args:
        ear (np.ndarray): Energy bin edges (keV).
        params (np.ndarray): Model parameters.

    Returns:
        np.ndarray: Photon flux.
    '''
    return gen_wrap(ear, params, wsim_dist)

