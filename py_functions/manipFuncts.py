"""
Functions for data manipulation
==============================

"""
import numpy as np

def nanhelper(y):
    """Helper to handle indices and logical indices of NaNs.
    Input:
    - y, 1d numpy array with possible NaNs
    Output:
    - nans, logical indices of NaNs
    - index, a function, with signature indices= index(logical_indices),
    to convert logical indices of NaNs to 'equivalent' indices
    Example:
    >>> # linear interpolation of NaNs
    >>> nans, x= nan_helper(y)
    >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """
    return np.isnan(y), lambda z: z.nonzero()[0]

def int2list(var):
    if not isinstance(var,list):
        var=[var]
    return var
