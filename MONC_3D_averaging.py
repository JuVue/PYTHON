### manual averaging of 4D monc data and save it to output for future use
from  __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import glob
from IPython import embed
from scipy.interpolate import interp1d
import matplotlib.cm as mpl_cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,LogNorm
import pickle
## import python functions
import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcAirDensity, calcThetaE, calcThetaVL, calcT
#from pyFixes import py3_FixNPLoad
############################################################
############################################################
start = time.time()

machine='JASMIN'
#machine='LEEDS'
if machine=='LEEDS':
    monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'
    m_out_dir = '22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/'
elif machine=='JASMIN':
    monc_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
    m_out_dir = '22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/'
    monc_exp_dir= '/gws/nopw/j04/ncas_radar_vol1/jutta/MONC/output/'
monc_filename=glob.glob(monc_root_dir + m_out_dir +'*.nc')
start = time.time()

  #################################################################
print ('Loading MONC data:')
print ('')
###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)

# monc_var_list =[['time_series_30_600','time_series_30_60'],
#                 ['z','zn','prefn','thref','thinit','rho','rhon','theta_mean',
#                 'liquid_mmr_mean','ice_mmr_mean','snow_mmr_mean','graupel_mmr_mean'],
#                 ['u','v','w','q','th','p','q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]
monc_var_list =[['time_series_30_600','time_series_30_60'],
                ['z','zn','prefn','thref','th','p','q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass'],
                ['u','v','w','q_vapour']]

monc_var_avg= ['q_cloud_liquid_mass','q_ice_mass','q_snow_mass','q_graupel_mass','iwc_tot','lwc_tot','twc_tot'] #always have twc_toto at the end!

# ml2  =        ['liquid_mmr_mean','ice_mmr_mean','snow_mmr_mean','graupel_mmr_mean','model_twc']



ncm = {}
monc_data = {}
print(monc_filename)
ncm = Dataset(monc_filename[0],'r')
monc_data={}
zvar={}
tvar={}
for c in range(0,len(monc_var_list)):
    for j in range(0,len(monc_var_list[c])):
        var = monc_var_list[c][j]
        zvar[var]=[]
        tvar[var]=[]
        if c == len(monc_var_list)-1:
            monc_data[var+'_mean'] = np.nanmean(ncm.variables[var][:],axis=(1,2))
        else:
            monc_data[var] = ncm.variables[var][:]
        ###getting right z and t dimensions
        tmp=ncm.variables[var].dimensions
        if "'z'" in str(tmp):
            zvar[var]='z'
        elif "'zn'" in str(tmp):
            zvar[var]='zn'
        else:
            zvar[var]=np.nan
        if monc_var_list[0][0] in str(tmp):
            tvar[var]='time1'
        elif monc_var_list[0][1] in str(tmp):
            tvar[var]='time2'
        if len(monc_var_list[0])>2:
            if monc_var_list[0][2] in str(tmp):
                tvar[var]='time3'

monc_data['time1']=monc_data[monc_var_list[0][0]] #1d data
monc_data.pop(monc_var_list[0][0])
monc_data['time2']=monc_data[monc_var_list[0][1]] #2d data
monc_data.pop(monc_var_list[0][1])
if len(monc_var_list[0])>2:
    monc_data['time3']=monc_data[monc_var_list[0][2]] #2d data
    monc_data.pop(monc_var_list[0][2])

print('Loading done')
## averaging 4D variables
# x = 50m
# y = 50m

##################################
print('Calculating temperature and density')
th = ncm.variables['th'][:]+ ncm.variables['thref'][0,:]
p = ncm.variables['p'][:]+ ncm.variables['prefn'][0,:]

monc_data.pop('th')
monc_data.pop('thref')
monc_data.pop('p')
monc_data.pop('prefn')

T = calcT(th,p)
rho=calcAirDensity(T,p)

print('Calculating mean values:T,th,p,rho')
monc_data['p_mean'] = np.nanmean(p,axis=(1,2))
monc_data['T_mean'] = np.nanmean(T,axis=(1,2))
monc_data['th_mean'] = np.nanmean(th,axis=(1,2))
monc_data['rho_mean'] = np.nanmean(rho,axis=(1,2))

zvar['p_mean'] = zvar['p']
zvar['T_mean'] = zvar['p']
zvar['th_mean'] = zvar['th']
zvar['rho_mean'] = zvar['p']
tvar['p_mean'] = tvar['p']
tvar['T_mean'] = tvar['p']
tvar['th_mean'] = tvar['th']
tvar['rho_mean'] = tvar['p']
del(p,T,th)



print('done')
###############
print('')
print('calculating twc')
#monc_data['model_iwc']= (monc_data['ice_mmr_mean']+monc_data['graupel_mmr_mean']+monc_data['snow_mmr_mean'])*monc_data['rho_mean']
#monc_data['model_lwc']= monc_data['liquid_mmr_mean']*monc_data['rho_mean']
#monc_data['model_twc'] = monc_data['model_lwc'] +monc_data['model_iwc']
#tvar['model_twc']=tvar['ice_mmr_mean']
#zvar['model_twc']= zvar['ice_mmr_mean']
monc_data['iwc_tot']=monc_data['q_ice_mass']+monc_data['q_snow_mass']+monc_data['q_graupel_mass']
monc_data['iwc_tot']=monc_data['iwc_tot']*rho
monc_data['lwc_tot']=monc_data['q_cloud_liquid_mass']*rho
monc_data['twc_tot']=monc_data['iwc_tot']+monc_data['lwc_tot']
#monc_data['twc_tot']=monc_data['q_ice_mass']+monc_data['q_snow_mass']+monc_data['q_graupel_mass']+monc_data['q_cloud_liquid_mass']
#monc_data['twc_tot']=monc_data['twc_tot']*rho
tvar['twc_tot']=tvar['q_ice_mass']
zvar['twc_tot']= zvar['q_ice_mass']
tvar['iwc_tot']=tvar['q_ice_mass']
zvar['iwc_tot']= zvar['q_ice_mass']
tvar['lwc_tot']=tvar['q_ice_mass']
zvar['lwc_tot']= zvar['q_ice_mass']

print('setting up thresholds')
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
a=monc_data[zvar[var]][monc_intZs]
twc_thresh[monc_intZs] = f(a.tolist())
twc_thresh = np.squeeze(np.array([[[twc_thresh]*monc_data['twc_tot'].shape[2]]*monc_data['twc_tot'].shape[1]]))

print('averaging cloud variables')
for c in range(0,len(monc_var_avg)):
    var = monc_var_avg[c]
    #calculate mean values
    monc_data[var +'_mean']=np.empty((np.size(monc_data[tvar[var]],0),np.size(monc_data[zvar[var]],0)))
    monc_data[var +'_mean'][:]=np.nan
    for t in range(0,np.size(monc_data[tvar[var]],0)):
        monc_data[var][t,monc_data['twc_tot'][t,:]<twc_thresh] = np.nan
        monc_data[var +'_mean'][t,:] = np.nanmean(monc_data[var][t,:],axis=(0,1))
        zvar[var +'_mean']=zvar[var]
        tvar[var +'_mean']=tvar[var]
    monc_data[var +'_mean'][np.isnan(monc_data[var +'_mean'])]=0.0
    monc_data.pop(var)

monc_data['zvar']={}
monc_data['tvar']={}
for var in zvar.keys():
    if var in monc_data.keys():
        monc_data['zvar'][var]=zvar[var]
        monc_data['tvar'][var]=tvar[var]



end = time.time()
print(end - start)
print (' averaging complete!')
if not os.path.exists(glob.glob(monc_exp_dir + m_out_dir)):
        os.mkdir(glob.glob(monc_exp_dir + m_out_dir))
fname=(monc_exp_dir + m_out_dir +'3d_averages')

np.save(fname +'.npy',monc_data)
#
# afile=open(fname,'wb')
# pickle.dump(monc_data,afile)
# afile.close()


# ######################################################
# ### plots to check averaging worked
# ##checking calculations
# viridis = mpl_cm.get_cmap('viridis', 256)
# newcolors = viridis(np.linspace(0, 1, 256))
# greyclr = np.array([0.1, 0.1, 0.1, 0.1])
# newcolors[:10, :] = greyclr
# newcmp = ListedColormap(newcolors)
#
#
# ###################
# #### plot air density
# var = 'rho'
# fig=plt.figure(figsize=(6,9))
# plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
#     hspace = 0.4, wspace = 0.2)
# plt.subplot(311)
# ax = plt.gca()
# img=plt.contourf(monc_data[tvar[var]],monc_data[zvar[var]],np.transpose(monc_data[var]),
#     cmap = newcmp)
# cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
# cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
# plt.title(var)
# plt.subplot(312)
# plt.contourf(monc_data[tvar[var]],monc_data[zvar[var]],np.transpose(monc_data[var +'_mean']),
#     cmap = newcmp)
# plt.title(var +'_mean')
# plt.subplot(313)
# plt.contourf(monc_data[tvar[var]],monc_data[zvar[var]],np.transpose(monc_data[var ]-monc_data[var +'_mean']),
#     cmap = newcmp)
# plt.title(var +'_diff')
#
# plots_out_dir='./'
# fileout = plots_out_dir + var + '_comparison.png'
# plt.savefig(fileout)
#
#
# clevs=range(265,295,2)
#
# fig=plt.figure(figsize=(6,9))
# plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
#     hspace = 0.4, wspace = 0.2)
# plt.subplot(311)
# ax = plt.gca()
# img=plt.contourf(monc_data[tvar['theta_mean']],monc_data[zvar['theta_mean']],np.transpose(monc_data['th_mean']),
#     levels=clevs,cmap = newcmp)
# cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
# cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
# plt.title(' own th_mean calc')
# plt.subplot(312)
# plt.contourf(monc_data[tvar['theta_mean']],monc_data[zvar['theta_mean']],np.transpose(monc_data['theta_mean']),
#     levels=clevs,cmap = newcmp)
# plt.title(var )
#
# clevels=np.arange(-1,1,0.1)
# plt.subplot(313)
# img2=plt.contourf(monc_data[tvar['theta_mean']],monc_data[zvar['theta_mean']],np.transpose(monc_data['th_mean']-monc_data['theta_mean']),
#     levels=clevels,cmap = newcmp)
# cbaxes = fig.add_axes([0.225, 0.35, 0.6, 0.015])
# cb = plt.colorbar(img2, cax = cbaxes, orientation = 'horizontal')
# plt.title(var +'_diff')
# plots_out_dir='./'
# fileout = plots_out_dir +'theta_comparison.png'
# plt.savefig(fileout)
#
#
# viridis = mpl_cm.get_cmap('viridis', 256)
# newcolors = viridis(np.linspace(0, 1, 256))
# greyclr = np.array([0.1, 0.1, 0.1, 0.1])
# newcolors[:10, :] = greyclr
# newcmp = ListedColormap(newcolors)
#
# clevs=[1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
# # clevs=[1e-5,0.25e-4,0.5e-4,0.75e-4,
# #     1e-4,0.25e-3,0.5e-3,0.75e-3,
# #     1e-3,0.25e-2,0.5e-2,0.75e-2,
# #     1e-2,0.25e-1,0.5e-1,0.75e-1,
# #     1e-1,0.25e-0,0.5e-0,0.75e-0,1e-0]
#
#
# for c in range(0,len(monc_var_avg)):
#     var = monc_var_avg[c]
#     fig=plt.figure(figsize=(6,9))
#     plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
#         hspace = 0.4, wspace = 0.2)
#     plt.subplot(2,1,1)
#     ax = plt.gca()
#     img=plt.contourf(monc_data[tvar[ml2[c]]],monc_data[zvar[ml2[c]]],np.transpose(monc_data[ml2[c]]),
#         levels=clevs, norm = LogNorm(),
#         cmap = newcmp)
#     cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
#     cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
#     plt.title(ml2[c])
#
#     plt.subplot(2,1,2)
#     plt.contourf(monc_data[tvar[var]],monc_data[zvar[var]],np.transpose(monc_data[var +'_mean']),
#         levels=clevs, norm = LogNorm(),
#         cmap = newcmp)
#     plt.title(var)
#     plots_out_dir='./'
#     fileout = plots_out_dir + var + '_comparison.png'
#     plt.savefig(fileout)
