### testing manual averaging of monc data


from  __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os


## import python functions
import sys
sys.path.insert(1, '/nfs/a26/lecimb/jutta/GITHUB/PYTHON/py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
#from physFuncts import calcThetaE, calcThetaVL
#from pyFixes import py3_FixNPLoad

monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'
m_out_dir = '3_control_20180913T0000Z/'


monc_filename= monc_root_dir + m_out_dir + 'moccha_casim_dg_72000.nc'

  #################################################################
print ('Loading MONC data:')
print ('')
###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)
monc_var_list =[['time_series_2_60','time_series_20_600' ,'z','rho', 'LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
              ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
              'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
              'graupel_mmr_mean'],
              ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp'],
              ['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]

ncm = {}
monc_data = {}
print(monc_filename)
  ncm = Dataset(monc_filename,'r')
  for c in range(0,len(monc_var_list)):
      for j in range(0,len(monc_var_list[c])):
          monc_data[monc_var_list[c][j]] = ncm.variables[monc_var_list[c][j]][:]

  monc_data['time2']=monc_data['time_series_20_600'] #2d data
  monc_data['time1']=monc_data['time_series_2_60'] #1d data
  monc_data.pop('time_series_2_60')
  monc_data.pop('time_series_20_600')

print (' Monc data Loaded!')

embed()
