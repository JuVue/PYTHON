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
    if np.issubdtype(type(var),np.integer):
        var=[var]
    return var

def intersect_mtlb(a, b):
    a1, ia = np.unique(a, return_index=True)
    b1, ib = np.unique(b, return_index=True)
    aux = np.concatenate((a1, b1))
    aux.sort()
    c = aux[:-1][aux[1:] == aux[:-1]]
    return c, ia[np.isin(a1, c)], ib[np.isin(b1, c)]
