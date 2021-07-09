### testing manual averaging of monc data


from  __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os

from IPython import embed

## import python functions
import sys
sys.path.insert(1, '/nfs/a26/lecimb/jutta/GITHUB/PYTHON/py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcAirDensity, calcThetaE, calcThetaVL, calcT
#from pyFixes import py3_FixNPLoad

monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'
m_out_dir = '3_control_20180913T0000Z/'


monc_filename= monc_root_dir + m_out_dir + 'moccha_casim_dg_72000.nc'
start = time.time()

  #################################################################
print ('Loading MONC data:')
print ('')
###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)
# monc_var_list =[['time_series_2_60','time_series_2_600','time_series_20_600' ,'z','zn', 'LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
#               ['rho','rhon','theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
#               'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
#               'graupel_mmr_mean'],
#               ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp']]
#               #['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]
monc_var_list =[['time_series_2_60','time_series_2_600','time_series_20_600' ,'z','zn','prefn','thref','thinit'],
              ['rho','rhon','theta_mean','w_wind_mean','u_wind_mean','v_wind_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
              'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
              'graupel_mmr_mean']]

#separate list for 4 d variables
monc_var_list_2  = ['q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']


ncm = {}
monc_data = {}
print(monc_filename)
ncm = Dataset(monc_filename,'r')
monc_data={}
zvar={}
tvar={}
for c in range(0,len(monc_var_list)):
    for j in range(0,len(monc_var_list[c])):
        var = monc_var_list[c][j]
        zvar[var]=[]
        tvar[var]=[]
        monc_data[var] = ncm.variables[var][:]
        ###getting right z and t dimensions
        tmp=ncm.variables[var].dimensions
        if "'z'" in str(tmp):
            zvar[var]='z'
        elif "'zn'" in str(tmp):
            zvar[var]='zn'
        else:
            zvar[var]=np.nan
        if "'time_series_2_60'" in str(tmp):
            tvar[var]='time1'
        elif "'time_series_2_600'" in str(tmp):
            tvar[var]='time2'
        elif "'time_series_20_600'" in str(tmp):
            tvar[var]='time3'
monc_data['zvar']=zvar
monc_data['tvar']=tvar

monc_data['time3']=monc_data['time_series_20_600'] #2d data
monc_data['time2']=monc_data['time_series_2_600'] #2d data
monc_data['time1']=monc_data['time_series_2_60'] #1d data
monc_data.pop('time_series_2_60')
monc_data.pop('time_series_2_600')
monc_data.pop('time_series_20_600')

## averaging 4D variables
# x = 50m
# y = 50m
th = ncm.variables['th'][:]+ ncm.variables['thinit'][0,:]
p=nncm.variables['prefn'][0,:]
p = np.array([[[p]*th.shape[2]]*th.shape[1]]*th.shape[0])
T = calcT(th,p)
monc_data['p_mean'] = np.nanmean(p,axis=(1,2))
monc_data['temp_mean'] = np.nanmean(T,axis=(1,2))
monc_data['th_mean'] = np.nanmean(th,axis=(1,2))
del(p,T,th)


for c in range(0,len(monc_var_list2)):
    var = monc_var_list2[c]
    zvar[var]=[]
    tvar[var]=[]
    ###getting right z and t dimensions
    tmp=ncm.variables[var].dimensions
    if "'z'" in str(tmp):
        zvar[var]='z'
    elif "'zn'" in str(tmp):
        zvar[var]='zn'
    else:
        zvar[var]=np.nan
    if "'time_series_2_60'" in str(tmp):
        tvar[var]='time1'
    elif "'time_series_2_600'" in str(tmp):
        tvar[var]='time2'
    elif "'time_series_20_600'" in str(tmp):
        tvar[var]='time3'
    #calculate mean values
    tmp = ncm.variables[var][:]
    tmp[tmp<=0.0]=np.nan
    monc_data[var +'_mean'] = np.nanmean(tmp,axis=(1,2))
    monc_data[var +'_mean'][np.isnan(monc_data[var +'_mean'])]=0.0
    del(tmp)




end = time.time()
print(end - start)
print (' Monc data Loaded!')
embed()

##################################
#### AVERAGING DATA MANUALLY
##################################
#
# ### 3D variables
# dnew={}
# vars=['lwp','iwp','swp','gwp']
#
# for c in vars:
#     monc_data[c][monc_data[c]<=0.0]=np.nan
#     dnew[c.upper() +'_mean'] = np.nanmean(monc_data[c],axis=(1,2))
#     dnew[c.upper() +'_mean'][np.isnan(dnew[c.upper() +'_mean'])]=0.0
#
# vars = list(dnew.keys())
#
# tvarid='time_series_2_60'
# m=1
# for c in vars:
#     plt.subplot(len(vars),1,m)
#     plt.plot(monc_data[tvarid],monc_data[c],'r')
#     plt.plot(monc_data['time_series_2_600'],dnew[c],'bo')
#     m+=1


### 4D variables
vars=['u','v','w','th']
embed()
for c in vars:
    dnew[c.upper() +'_mean'] = np.nanmean(monc_data[c],axis=(1,2))
    dnew[c.upper() +'_mean'][np.isnan(dnew[c.upper() +'_mean'])]=0.0
