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
from scipy.interpolate import interp1d

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

monc_var_list=[['time_series_2_60','time_series_2_600','time_series_20_600' ,'z','zn'],
                ['q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]

monc_var_list =[['time_series_2_60','time_series_2_600','time_series_20_600' ,'z','zn','prefn','thref','thinit'],
              ['rho','rhon','theta_mean','w_wind_mean','u_wind_mean','v_wind_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
              'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
              'graupel_mmr_mean'],
              ['q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]

ml2  =        ['liquid_mmr_mean','ice_mmr_mean','snow_mmr_mean','graupel_mmr_mean']
monc_var_avg= ['q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass']
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
# th = ncm.variables['th'][:]+ ncm.variables['thinit'][0,:]
# p=ncm.variables['prefn'][0,:]
# p = np.array([[[p]*th.shape[2]]*th.shape[1]]*th.shape[0])
# T = calcT(th,p)
# monc_data['p_mean'] = np.nanmean(p,axis=(1,2))
# monc_data['temp_mean'] = np.nanmean(T,axis=(1,2))
# monc_data['th_mean'] = np.nanmean(th,axis=(1,2))
# del(p,T,th)

monc_data['twc_tot']=monc_data['q_ice_mass']+monc_data['q_snow_mass']+monc_data['q_graupel_mass']+monc_data['q_cloud_liquid_mass']
#calculate mean values
var='q_ice_mass'
twc_thresh = np.zeros([np.size(monc_data[zvar[var]],0)])
monc_lt1km = np.where(monc_data[zvar[var]][:]<=1e3)
twc_thresh[monc_lt1km] =  1e-6
monc_lt1km = np.where(monc_data[zvar[var]][:]>=4e3)
twc_thresh[monc_lt1km] = 1e-7
monc_intZs = np.where(twc_thresh == 0.0)
x = [1e-6, 1e-7]
y = [1e3, 4e3]
f = interp1d(y, x)
twc_thresh[monc_intZs] = f(monc_data[zvar[var]][monc_intZs])
twc_thresh = np.squeeze(np.array([[[twc_thresh]*monc_data['twc_tot'].shape[2]]*monc_data['twc_tot'].shape[1]]))

for c in range(0,len(monc_var_avg)):
    var = monc_var_avg[c]
    #calculate mean values
    monc_data[var +'_mean']=np.empty((np.size(monc_data[tvar[var]],0),np.size(monc_data[zvar[var]],0)))
    monc_data[var +'_mean'][:]=np.nan
    for t in range(0,np.size(monc_data[tvar[var]],0)):
        #tt=monc_data[var][t,:]
        #twc=monc_data['twc_tot'][t,:]
        #tt[twc<twc_thresh] = np.nan
        #monc_data[var +'_mean'][t,:] = np.nanmean(tt,axis=(0,1))
        monc_data[var][t,monc_data['twc_tot'][t,:]<twc_thresh] = np.nan
        embed()
        monc_data[var +'_mean'][t,:] = np.nanmean(monc_data[var][t,:],axis=(0,1))
    monc_data[var +'_mean'][np.isnan(monc_data[var +'_mean'])]=0.0
    monc_data.pop(var)


end = time.time()
print(end - start)
print (' Monc data Loaded!')
import matplotlib.cm as mpl_cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,LogNorm

viridis = mpl_cm.get_cmap('viridis', 256)
newcolors = viridis(np.linspace(0, 1, 256))
greyclr = np.array([0.1, 0.1, 0.1, 0.1])
newcolors[:10, :] = greyclr
newcmp = ListedColormap(newcolors)

clevs=[1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3]
clevs=[1e-5,0.25e-4,0.5e-4,0.75e-4,
    1e-4,0.25e-3,0.5e-3,0.75e-3,
    1e-3,0.25e-2,0.5e-2,0.75e-2,
    1e-2,0.25e-1,0.5e-1,0.75e-1,
    1e-1,0.25e-0,0.5e-0,0.75e-0,1e-0]

    # fig=plt.figure(figsize=(6,9))
    # plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
    #     hspace = 0.4, wspace = 0.2)
    # plt.subplot(2,1,1)
    # img=plt.contourf(monc_data[tvar[ml2[c]]],monc_data[zvar[ml2[c]]],np.transpose(monc_data['q_cloud_liquid_mass_mean']),
    #     levels=clevs, norm = LogNorm(),
    #     cmap = newcmp)
    # cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    # cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    # plt.show()


for c in range(0,len(monc_var_avg)):
    var = monc_var_avg[c]
    fig=plt.figure(figsize=(6,9))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
        hspace = 0.4, wspace = 0.2)
    plt.subplot(2,1,1)
    ax = plt.gca()
    img=plt.contourf(monc_data[tvar[ml2[c]]],monc_data[zvar[ml2[c]]],np.transpose(monc_data[ml2[c]]),
        levels=clevs, norm = LogNorm(),
        cmap = newcmp)
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title(ml2[c])

    plt.subplot(2,1,2)
    plt.contourf(monc_data[tvar[var]],monc_data[zvar[var]],np.transpose(monc_data[var +'_mean']),
        levels=clevs, norm = LogNorm(),
        cmap = newcmp)
    plt.title(var)
    plt.show()
embed()

# fig=plt.figure(figsize=(6,9))
# plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
#     hspace = 0.4, wspace = 0.2)
# plt.subplot(2,1,1)
# ax = plt.gca()
# img=plt.contourf(monc_data[tvar[ml2[c]]],monc_data[zvar[ml2[c]]],np.transpose(monc_data['twc_tot_mean']),
#     levels=clevs, norm = LogNorm(),
#     cmap = newcmp)
# plt.show()
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
