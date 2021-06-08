"""
Functions to convert timestamps
==============================

"""

from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from datetime import datetime as dt
from datetime import timedelta
import panda as pd

def date2datenum(d):
 ##  Convert python datetime to Matlab datenum
 ##:    input: Date as python datetime object
 ##:    return: corresponding datenum
    return 366 + d.toordinal() + (d - dt.fromordinal(d.toordinal())).total_seconds()/(24*60*60)

def datenum2date(datenum):
   ##Convert Matlab datenum into Python datetime.
   ##:input  datenum Date in datenum format
   ##:return:Datetime object corresponding to datenum.
   d=np.array(datenum)
   return pd.to_datetime(d-719529,unit='D')
    # if isinstance(datenum,list):
    #     for m in range(0,len(datenum)):
    #         d=datenum[m]
    #         d = np.array(d)
    #         dout[m] =np.append(dout, dt.fromordinal(int(d.astype(int)))+timedelta(days=d%1)-timedelta(days=366))
    # else:
    #     d=datenum
    #     d = np.array(d)
    #     dout[0] =  dt.fromordinal(int(d.astype(int)))+timedelta(days=d%1)-timedelta(days=366))
    # return dout


def calcTime_Mat2DOY(matlab_time):
        #### EXAMPLE OF USE:
        #### pytime = calcTime_Mat2DOY(matlab_time)

    print ('Converting MATLAB timesteps to DOY:')

    timestamps = pd.to_datetime(matlab_time-719529, unit='D')
    python_time = timestamps.dayofyear + (timestamps.hour / 24.0) + (timestamps.minute / 1440.0) + (timestamps.second / 86400.0)

    return python_time

def calcTime_Date2DOY(date):
    ####        date should be forematted as YYYYmmDD

    print ('Converting date to DOY:')

    mm = date[4:6]      #### month
    DD = date[6:8]      #### day
    refDateAug = 226       #### Aug reference date for drift: 14th Aug 2018
    refDateSep = 243       #### Sep reference date for drift: 1st Sep 2018

    if mm == '08':
        doy = (float(DD) - 14.0) + refDateAug
    elif mm == '09':
        doy = float(DD) + refDateSep
    else:
        print ('****Date not valid with this function****')

    print ('----')
    print ('Date = ', date)
    print ('DOY = ', doy)
    print ('')

    return doy
