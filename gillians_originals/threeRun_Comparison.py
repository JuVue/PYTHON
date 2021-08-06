###
###
### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
###  loads in the raw model data (UM + IFS) for e.g., Fig 2 (radiation timeseries) in the paper
###

from __future__ import print_function
import time
import datetime
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import diags_MOCCHA as diags
import diags_varnames as varnames
import cartopy.crs as ccrs
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import os
import seaborn as sns

#### import python functions
import sys
sys.path.insert(1, '../py_functions/')
from time_functions import calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcThetaE, calcThetaVL
from pyFixes import py3_FixNPLoad

def readfile(filename):

    import pandas as pd

    # print '******'
    print( '')
    print( 'Reading .txt file with pandas')
    print( '')

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data

def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

def iceDrift(data):

    ###################################
    ## Define ice drift period
    ###################################

    Aug_drift_index = np.where(np.logical_and(data.values[:,2]>=14,data.values[:,1]==8))
    Sep_drift_index = np.where(np.logical_and(np.logical_and(data.values[:,2]<=14,data.values[:,1]==9),data.values[:,3]<=22))
    drift_index = range(Aug_drift_index[0][0],Sep_drift_index[0][-1])

    print ('******')
    print ('')
    # print 'Aug drift: ' + str(data.values[Aug_drift_index[0][0],0:3]) + ' - ' + str(data.values[Aug_drift_index[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_drift_index[0][0],0:3]) + ' - ' + str(data.values[Sep_drift_index[0][-1],0:3])
    print ('Whole drift: ' + str(data.values[drift_index[0],0:4]) + ' - ' + str(data.values[drift_index[-1],0:4]))
    print ('')

    return drift_index

def inIce(data):

    ###################################
    ## DEFINE IN ICE PERIOD
    ###################################
    Aug_inIce = np.where(np.logical_and(data.values[:,2]>=3,data.values[:,1]==8))
    Sep_inIce = np.where(np.logical_and(data.values[:,2]<20,data.values[:,1]==9))
    inIce_index = np.arange(Aug_inIce[0][0],Sep_inIce[0][-1])

    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    # Aug_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=12,data.values[:,1]==8),data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]>=13,data.values[:,1]==8),data.values[:,3]>=0))
    # # Sep_inIce = np.where(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9))
    # # Sep_inIce = np.where(np.logical_and(np.logical_and(data.values[:,2]<=20,data.values[:,1]==9),data.values[:,3]<=1))
    # inIce_index = range(Aug_inIce[0][0],Sep_inIce[0][-1])

    print( '******')
    print( '')
    # print 'Aug drift: ' + str(data.values[Aug_inIce[0][0],0:3]) + ' - ' + str(data.values[Aug_inIce[0][-1],0:3])
    # print 'Sep drift: ' + str(data.values[Sep_inIce[0][0],0:3]) + ' - ' + str(data.values[Sep_inIce[0][-1],0:3])
    # print 'In ice: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4])
    print( 'CloudNET: ' + str(data.values[inIce_index[0],0:4]) + ' - ' + str(data.values[inIce_index[-1],0:4]))
    print( '')
    print( 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')')
    print( 'Lon/lat of start point: (' + str(data.values[inIce_index[0],6]) + ', ' + str(data.values[inIce_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[inIce_index[-1],6]) + ', ' + str(data.values[inIce_index[-1],7]) + ')')
    print( 'Min/max longitude: ' + str(np.nanmin(data.values[inIce_index,6])) + ', ' + str(np.nanmax(data.values[inIce_index,6])))
    print( 'Min/max latitude: ' + str(np.nanmin(data.values[inIce_index,7])) + ', ' + str(np.nanmax(data.values[inIce_index,7])))
    print( '')

    return inIce_index

def trackShip(data):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==25,data.values[:,1]==8),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print( '******')
    print( '')
    print( 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    print( 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print( '')

    return trackShip_index

def plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined 1d timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    # datenums_obs['foremast'] = obs['foremast'].variables['time'][:] ### obs['foremast'] data on different timestep
    # time_obs['foremast'] = calcTime_Mat2DOY(datenums_obs['foremast'])

    # datenums_obs['deck7th'] = obs['deck7th'].variables['time'][:] ### 7th deck data on different timestep
    # time_obs['deck7th'] = calcTime_Mat2DOY(datenums_obs['deck7th'])

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)
    # ax = plt.gca()
    # plt.plot(data1['time'], data1['sfc_pressure'].data/1e2, color = 'darkblue', label = label1)
    # plt.plot(data2['time'], data2['sfc_pressure'].data/1e2, color = 'mediumseagreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'gold', label = label3)
    # else:
    #     plt.plot(data3['time'], data3['sfc_pressure'].data/1e2, color = 'gold',label = label3)
    # plt.plot(obs['foremast'].variables['doy'][:],obs['foremast'].variables['Psurf'][:], 'k')
    # plt.title('sfc_pressure [hPa]')
    # plt.ylim([980, 1010])
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'darkblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'mediumseagreen')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'gold', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'gold')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    # data1['surface_net_SW_radiation'].data[data1['surface_net_SW_radiation'].data == 0] = np.nan
    # data2['surface_net_SW_radiation'].data[data2['surface_net_SW_radiation'].data == 0] = np.nan
    # if out_dir3 != 'OUT_25H/': data3['surface_net_SW_radiation'].data[data3['surface_net_SW_radiation'].data == 0] = np.nan

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'mediumseagreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    # plt.subplot(3,3,6)
    # ax = plt.gca()
    # plt.plot(time_um, data['snowfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_snow'].data)
    # plt.title('sfc_snow_amount [kg/m2]')
    # plt.legend()
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    # plt.subplot(3,3,7)
    # ax = plt.gca()
    # plt.plot(time_um, data['rainfall_flux'].data)
    # plt.plot(time_um2, data2['sfc_ls_rain'].data)
    # plt.title('sfc_rain_amount [kg/m2]')
    # if month_flag == 8: ax.set_xlim([13.0, 31.0])
    # if month_flag == 9: ax.set_xlim([1.0, 15.0])
    # if month_flag == -1: ax.set_xlim([225.0, 258.0])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'mediumseagreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station_fluxes']['mday'],obs['ice_station_fluxes']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'mediumseagreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/244-256_oden_metum_ifs_casim-aeroprof_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined 1d timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 18

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(8,9))
    plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.9, left = 0.1,
            hspace = 0.4, wspace = 0.15)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    # time_um1 = data1['time'][:]
    # time_um2 = data2['time'][:]
    # time_um3 = data3['time'][:]

    #################################################################
    ## create figure and axes instances
    #################################################################

    plt.subplot(211)
    ax = plt.gca()
    plt.plot(data1['time'][:], data1['temp_1.5m'].data - 273.15, color = 'darkblue', label = label1)
    plt.plot(data2['time'][:], data2['temp_1.5m'].data - 273.15, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_temp_2m'].data - 273.15, color = 'gold', label =  label3)
    else:
        plt.plot(data3['time'][:], data3['temp_1.5m'].data - 273.15, color = 'gold')#, label = '2m')
    plt.plot(time_temp,obs['obs_temp'].variables['Tice'][:] - 273.15, color = 'black', label = 'Observations')
    plt.legend()
    plt.title('Temperature [$^{o}C$]')
    plt.ylim([263 - 273,275 - 273])
    # plt.grid('on')
    if month_flag == 8:
        ax.set_xlim([13.0, 31.0])
        plt.xlabel('Day of month [Aug]')
    if month_flag == 9:
        ax.set_xlim([1.0, 15.0])
        plt.xlabel('Day of month [Sep]')
    if month_flag == -1:
        ax.set_xlim([doy[0],doy[-1]])
        # plt.xlabel('Day of year')

    plt.subplot(2,1,2)
    ax = plt.gca()
    plt.plot(data1['time'][:], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'][:], data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'][:], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'][:], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.title('Net SW radiation [W/m2]')
    plt.ylim([0,100])
    # plt.grid('on')
    if month_flag == 8:
        ax.set_xlim([13.0, 31.0])
        plt.xlabel('Day of month [Aug]')
    if month_flag == 9:
        ax.set_xlim([1.0, 15.0])
        plt.xlabel('Day of month [Sep]')
    if month_flag == -1:
        ax.set_xlim([doy[0],doy[-1]])
        plt.xlabel('Day of year')

    # plt.subplot(3,1,3)
    # ax = plt.gca()
    # # data['surface_net_LW_radiation'].data[data['surface_net_LW_radiation'].data == 0] = np.nan
    # plt.plot(data1['time'][:], data1['surface_net_LW_radiation'].data, color = 'darkblue', label = label1)
    # plt.plot(data2['time'][:], data2['surface_net_LW_radiation'].data, color = 'mediumseagreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'][:], data3['sfc_net_lw'].data, color = 'gold')
    # else:
    #     plt.plot(data3['time'][:], data3['surface_net_LW_radiation'].data, color = 'gold')
    # plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'Observations')
    # # plt.legend()
    # plt.title('Net LW radiation [W/m2]')
    # # plt.ylim([260,275])
    # # plt.grid('on')
    # if month_flag == 8:
    #     ax.set_xlim([13.0, 31.0])
    #     plt.xlabel('Day of month [Aug]')
    # if month_flag == 9:
    #     ax.set_xlim([1.0, 15.0])
    #     plt.xlabel('Day of month [Sep]')
    # if month_flag == -1:
    #     ax.set_xlim([doy[0],doy[-1]])
    #     plt.xlabel('Day of year')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        if out_dir2[0:20] == '6_u-bm410_RA1M_CASIM':
            fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-200_tempoC_SW.png'
        elif out_dir2[0:20] == '5_u-bl661_RA1M_CASIM':
            if ifs_flag == True:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
            else:
                fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_casim-100_tempoC_SW.png'
        elif out_dir2[:18] == '4_u-bg610_RA2M_CON':
            fileout = '../FIGS/comparisons/' + out_dir1[:18] + '_oden_metum_tempoC_SW.png'
        else:
            fileout = '../FIGS/comparisons/' + out_dir2[:18] + '_oden_metum_tempoC_SW.png'
    # print 'Saving as: ' + fileout
    fileout = '../FIGS/comparisons/' + out_dir2[0:20] + '_oden_metum_metum-erai_IFS_tempoC_SW.svg'
    plt.savefig(fileout)#, dpi=400)
    plt.show()

def plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined 1d timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'steelblue', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'darkblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'steelblue')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'gold', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'gold')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'steelblue', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'steelblue')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'steelblue')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station_fluxes']['mday'],obs['ice_station_fluxes']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'steelblue')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/CASIM-NiceTest_C86-F62-M92_DOY243-249.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_CASIM_NdropTimeseries(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    from physFuncts import calcAirDensity
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Ndrop timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(13, 12))
    plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.88, left = 0.06,
            hspace = 0.3, wspace = 0.1)

    ### define axis instance
    # ax = plt.gca()

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    # ### convert radiation matlab timesteps to doy
    # datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    # time_radice_all = calcTime_Mat2DOY(datenums_radice)
    #
    # ### calculate net LW and SW from obs
    # time_radice = time_radice_all[:-1:2]
    # lwd = obs['obs_temp'].variables['LWdice'][:]
    # lwu = obs['obs_temp'].variables['LWuice'][:]
    # lwdmean = np.nanmean(lwd[:-1].reshape(-1, 2), axis=1)
    # lwumean = np.nanmean(lwu[:-1].reshape(-1, 2), axis=1)
    # swd = obs['obs_temp'].variables['SWdice'][:]
    # swu = obs['obs_temp'].variables['SWuice'][:]
    # swdmean = np.nanmean(swd[:-1].reshape(-1, 2), axis=1)
    # swumean = np.nanmean(swu[:-1].reshape(-1, 2), axis=1)
    # netLW = lwdmean - lwumean
    # netSW = swdmean - swumean

    ########            Radiation timeseries
    ########

    sw2 = data2['surface_downwelling_SW_radiation'][data2['hrly_flag']]# data2['surface_net_SW_radiation'][data2['hrly_flag']]
    lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    sw4 = data4['surface_downwelling_SW_radiation'][data4['hrly_flag']] # data4['surface_net_SW_radiation'][data4['hrly_flag']]
    lw4 = data4['surface_net_LW_radiation'][data4['hrly_flag']]
    crf4 = sw4# + lw4
    crf2 = sw2 #+ lw2

    ### ---------------------------------------------
    ### just downwelling radiation timeseries
    ### ---------------------------------------------
    # ax  = fig.add_axes([0.18,0.8,0.6,0.18])   # left, bottom, width, height
    # plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['SWd_ice'], color = 'grey', label = 'Ice_station')
    # plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['SWd_ship'], color = 'k', label = 'Ship')
    # plt.plot(data2['time'][data2['hrly_flag']], crf2, color = 'mediumseagreen', label = label2)
    # plt.plot(data4['time'][data4['hrly_flag']], crf4, color = 'purple', label = label4)
    # plt.xlim(doy[0], doy[-1])
    # plt.ylim([-10, 220])
    # plt.legend(bbox_to_anchor=(0.36, 0.79, 1., .102), loc=1, ncol=1)
    # plt.ylabel('SW$_{\downarrow}$ \n [W m$^{-2}$]')
    # plt.xlabel('Date')
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])

    ### ---------------------------------------------
    ### timeseries + pdf
    ### ---------------------------------------------
    ax  = fig.add_axes([0.1,0.81,0.62,0.18])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['SWd_ice'], color = 'grey', label = 'Ice_station')
    plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['SWd_ship'], color = 'k', label = 'Ship')
    plt.plot(data2['time'][data2['hrly_flag']], crf2, color = 'mediumseagreen', label = label2)
    plt.plot(data4['time'][data4['hrly_flag']], crf4, color = 'purple', label = label4)
    plt.xlim(doy[0], doy[-1])
    plt.ylim([-10, 240])
    plt.legend(bbox_to_anchor=(-0.05, 0.92, 1., .102), loc=1, ncol=2)
    plt.ylabel('SW$_{\downarrow}$ [W m$^{-2}$]')
    plt.xlabel('Date')
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])


    ax  = fig.add_axes([0.78,0.83,0.13,0.15])   # left, bottom, width, height
    yEmax = 0.02
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    f = sns.distplot(obs['fixed_radiation']['SWd_ice'], hist=False, color="grey", kde_kws={"linewidth": 2})
    f = sns.distplot(obs['fixed_radiation']['SWd_ship'], hist=False, color="black")
    f = sns.distplot(crf4, hist=False, color="purple", kde_kws={"shade": True})
    f = sns.distplot(crf2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # plt.annotate('Melt', xy=(87,0.087), xytext=(87,0.087), fontsize = 14)
    plt.xlim([-20,240])
    plt.ylim([0,yEmax])
    plt.yticks([0,0.01,0.02])
    plt.xlabel('SW$_{\downarrow}$ [W m$^{-2}$]')

    ########            Cloud fraction
    ########
    #### set up colourmaps to grey out zeros on figures
    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr
    newcmp = ListedColormap(newcolors)

    # plt.subplot(434)
    ax  = fig.add_axes([0.05,0.53,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    plt.contourf(data2['time'], data2['height'][:], np.transpose(data2['cloud_fraction']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,3.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.title(label2, fontsize = 14)# + '\n Cloud fraction')
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])

    # plt.subplot(435)
    ax  = fig.add_axes([0.38,0.53,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    img = plt.contourf(data4['time'], data4['height'][:], np.transpose(data4['cloud_fraction']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.title(label4, fontsize = 14)# + '\n Cloud fraction')
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # cbaxes = fig.add_axes([0.9, 0.56, 0.015, 0.15])
    cbaxes = fig.add_axes([0.71, 0.54, 0.015, 0.15])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    cbaxes.set_xlabel('C$_{V}$', rotation = 0, labelpad = 15, fontsize = 14)
    cbaxes.xaxis.set_label_position('top')

    #### set flagged data to nans
    data2['cloud_fraction'][data2['cloud_fraction'] < 0.0] = np.nan
    data4['cloud_fraction'][data4['cloud_fraction'] < 0.0] = np.nan

    ax  = fig.add_axes([0.83,0.55,0.13,0.15])   # left, bottom, width, height
    ax = plt.gca()
    ax.fill_betweenx(data2['height'],np.nanmean(data2['cloud_fraction'],0) - np.nanstd(data2['cloud_fraction'],0),
        np.nanmean(data2['cloud_fraction'],0) + np.nanstd(data2['cloud_fraction'],0), color = 'mediumaquamarine', alpha = 0.15)
    ax.fill_betweenx(data4['height'],np.nanmean(data4['cloud_fraction'],0) - np.nanstd(data4['cloud_fraction'],0),
        np.nanmean(data4['cloud_fraction'],0) + np.nanstd(data4['cloud_fraction'],0), color = 'violet', alpha = 0.15)
    plt.plot(np.nanmean(data2['cloud_fraction'],0),data2['height'], color = 'mediumseagreen', linewidth = 3, label = label2)
    plt.plot(np.nanmean(data4['cloud_fraction'],0),data4['height'], color = 'purple', linewidth = 3, label = label4)
    plt.xlabel('C$_{V}$')
    # ax.xaxis.set_label_position('top')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([0,1])
    plt.xticks(np.arange(0,1.01,0.25))
    ax.set_xticklabels([0,' ',0.5,' ',1.0])

    ########            NDROP
    ########
    #### calculate air density in kg/m3
    data2['rho'] = calcAirDensity(data2['temperature'].data, data2['pressure'].data / 1e2)
    data4['rho'] = calcAirDensity(data4['temperature'].data, data4['pressure'].data / 1e2)

    print(data2.keys())

    ### convert /kg to /m3
    data2['qnliq'] = data2['qnliq'] * data2['rho']
    data4['qnliq'] = data4['qnliq'] * data4['rho']

    #### set flagged um_data to nans
    data2['qnliq'][data2['qnliq'] < 0] = 0.0
    data4['qnliq'][data4['qnliq'] < 0] = 0.0

    #### set up colourmaps to grey out zeros on figures
    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:10, :] = greyclr
    newcmp = ListedColormap(newcolors)

    # plt.subplot(425)
    ax  = fig.add_axes([0.05,0.3,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    plt.contourf(data2['time'], data2['height'][:], np.transpose(data2['qnliq'])/1e6,
        [0, 10, 50, 100, 150, 200, 250],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])

    # plt.subplot(426)
    ax  = fig.add_axes([0.38,0.3,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    img = plt.contourf(data4['time'], data4['height'][:], np.transpose(data4['qnliq'])/1e6,
        [0, 10, 50, 100, 150, 200, 250],
        # vmin = 0, vmax = 150,
        cmap = newcmp #mpl_cm.Blues
        )
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    cbaxes = fig.add_axes([0.71, 0.31, 0.015, 0.15])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    cbaxes.set_xlabel('N$_{d}$ \n[cm$^{-3}$]', rotation = 0, labelpad = 10, fontsize = 14)
    cbaxes.xaxis.set_label_position('top')

    #### set flagged data to nans
    data2['qnliq'][data2['qnliq'] <= 0.0] = np.nan
    data4['qnliq'][data4['qnliq'] <= 0.0] = np.nan

    ax  = fig.add_axes([0.83,0.32,0.13,0.15])   # left, bottom, width, height
    ax = plt.gca()
    ax.fill_betweenx(data2['height'],np.nanmean(data2['qnliq']/1e6,0) - np.nanstd(data2['qnliq']/1e6,0),
        np.nanmean(data2['qnliq']/1e6,0) + np.nanstd(data2['qnliq']/1e6,0), color = 'mediumaquamarine', alpha = 0.15)
    ax.fill_betweenx(data4['height'],np.nanmean(data4['qnliq']/1e6,0) - np.nanstd(data4['qnliq']/1e6,0),
        np.nanmean(data4['qnliq']/1e6,0) + np.nanstd(data4['qnliq']/1e6,0), color = 'violet', alpha = 0.15)
    plt.plot(np.nanmean(data2['qnliq']/1e6,0),data2['height'], color = 'mediumseagreen', linewidth = 3, label = label2)
    plt.plot(np.nanmean(data4['qnliq']/1e6,0),data4['height'], color = 'purple', linewidth = 3, label = label4)
    plt.xlabel('N$_{d}$ [cm$^{-3}$]')
    # ax.xaxis.set_label_position('top')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([0,150])
    plt.xticks(np.arange(0,151,25))
    ax.set_xticklabels([0,' ',' ',75,' ',' ',150])
    # plt.legend()

    ########            QLIQ
    ########
    #### set flagged um_data to nans
    data2['qliq'][data2['qliq'] < 0] = 0.0
    data4['qliq'][data4['qliq'] < 0] = 0.0

    #### set up colourmaps to grey out zeros on figures
    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:10, :] = greyclr
    newcmp = ListedColormap(newcolors)

    # plt.subplot(427)
    ax  = fig.add_axes([0.05,0.07,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    plt.contourf(data2['time'], data2['height'][:], np.transpose(data2['qliq'])*1e3,
        [0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
        # vmin = 0, vmax = 0.35,
        cmap = newcmp
        )
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xlabel('Date')
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])

    # plt.subplot(428)
    ax  = fig.add_axes([0.38,0.07,0.3,0.18])   # left, bottom, width, height
    ax = plt.gca()
    img = plt.contourf(data4['time'], data4['height'][:], np.transpose(data4['qliq'])*1e3,
        [0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3],
        # vmin = 0, vmax = 0.35,
        cmap = newcmp
        )
    plt.xlabel('Date')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    cbaxes = fig.add_axes([0.71, 0.08, 0.015, 0.15])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    cbaxes.set_xlabel('q$_{liq}$ \n[g kg$^{-1}$]', rotation = 0, labelpad = 10, fontsize = 14)
    cbaxes.xaxis.set_label_position('top')

    #### set flagged data to nans
    data2['qliq'][data2['qliq'] <= 0.0] = np.nan
    data4['qliq'][data4['qliq'] <= 0.0] = np.nan

    ax  = fig.add_axes([0.83,0.09,0.13,0.15])   # left, bottom, width, height
    ax = plt.gca()
    ax.fill_betweenx(data2['height'],np.nanmean(data2['qliq']*1e3,0) - np.nanstd(data2['qliq']*1e3,0),
        np.nanmean(data2['qliq']*1e3,0) + np.nanstd(data2['qliq']*1e3,0), color = 'mediumaquamarine', alpha = 0.15)
    ax.fill_betweenx(data4['height'],np.nanmean(data4['qliq']*1e3,0) - np.nanstd(data4['qliq']*1e3,0),
        np.nanmean(data4['qliq']*1e3,0) + np.nanstd(data4['qliq']*1e3,0), color = 'violet', alpha = 0.15)
    plt.plot(np.nanmean(data2['qliq']*1e3,0),data2['height'], color = 'mediumseagreen', linewidth = 3, label = label2)
    plt.plot(np.nanmean(data4['qliq']*1e3,0),data4['height'], color = 'purple', linewidth = 3, label = label4)
    plt.xlabel('q$_{liq}$ [g kg$^{-1}$]')
    # ax.xaxis.set_label_position('top')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([0,0.1])
    plt.xticks(np.arange(0,0.11,0.025))
    ax.set_xticklabels([0,' ',0.05,' ',0.1])
    # plt.legend(bbox_to_anchor=(1.0, 0.0, 1., .102), loc=1, ncol=1)

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = '../FIGS/CASIM/CASIM-100_CASIM-AeroProf_SWdown-TS-PDFs_Cv_Ndrop_Qliq_hourlyObs_newColours_Dates_newRadiation.svg'
    # plt.savefig(fileout)
    plt.show()


    ###################################
    ## PLOT RADIOSONDE BIASES
    ###################################

    # print ('******')
    # print ('')
    # print ('Plotting radiosonde and model temperature profiles:')
    # print ('')
    #
    # #### change matlab time to doy
    # obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))
    #
    # ### set diagnostic naming flags for if IFS being used
    # if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
    #     ifs_flag = True
    # else:
    #     ifs_flag = False
    #
    # ### for reference in figures
    # zeros = np.zeros(len(data2['time']))
    #
    # #### set flagged values to nans
    # data1['temperature'][data1['temperature'] == -9999] = np.nan
    # data2['temperature'][data2['temperature'] == -9999] = np.nan
    # data3['temperature'][data3['temperature'] <= 0] = np.nan
    # data4['temperature'][data4['temperature'] == -9999] = np.nan
    # data5['temperature'][data5['temperature'] <= 0] = np.nan
    #
    # #### set flagged values to nans
    # data1['q'][data1['q'] == -9999] = np.nan
    # data2['q'][data2['q'] == -9999] = np.nan
    # data3['q'][data3['q'] <= 0] = np.nan
    # data4['q'][data4['q'] == -9999] = np.nan
    # data5['q'][data5['q'] <= 0] = np.nan
    #
    # #### ---------------------------------------------------------------
    # #### re-grid sonde and IFS data to UM vertical grid <10km
    # #### ---------------------------------------------------------------
    #
    # print ('...')
    # print ('Re-gridding sonde and ifs data...')
    # print ('')
    # data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'temp')
    # data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    # print ('')
    # print ('Done!')
    #
    # print ('')
    # print ('Starting radiosonde figure (quite slow!)...:')
    # print ('...')
    #
    # Tmin = -45
    # Tmax = 5
    # ymax = 9000
    # qmax = 4.0
    #
    # data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    #
    # data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    # data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    #
    # ### all model data share a timestamp
    # p3 = np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] < 230.0))
    # p4 = np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] < 240.0))
    # p5 = np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] < 247.0))
    # p6 = np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] < 251.0))
    #
    # ##################################################
    # ##################################################
    # #### MEDIAN PROFILES
    # ##################################################
    # ##################################################
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    #
    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(7.5,5.5))
    # plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
    #         hspace = 0.3, wspace = 0.2)
    #
    # ####        all model data share a timestamp
    # melt = np.where(data1['time_hrly'][::6] < 240.0)
    # freeze = np.where(data1['time_hrly'][::6] >= 240.0)
    #
    # ###-------------------------
    # plt.subplot(121)
    # ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
    # #     np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
    # #     color = 'blue', alpha = 0.05)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
    #     np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
    #     np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
    # #     np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
    # #     color = 'navajowhite', alpha = 0.35)
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    #
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'purple', label = label4, zorder = 2)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
    # #     np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
    # #     color = 'blue', alpha = 0.05)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
    # #     np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
    # #      color = 'lightblue', alpha = 0.15)
    # # plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'steelblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'steelblue', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
    # #     np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
    # #     color = 'mediumaquamarine', alpha = 0.15)
    # # plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
    # #     np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
    # #     color = 'navajowhite', alpha = 0.35)
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # #
    # # plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    # # plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # # plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4, zorder = 2)
    # # plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # plt.legend(bbox_to_anchor=(1.05, 1.03, 1., .102), loc=4, ncol=3)          ## 2 lines
    # # plt.legend(bbox_to_anchor=(1.25, 1.03, 1., .102), loc=4, ncol=3)            ## 5 lines
    # plt.ylabel('Z [km]')
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-2.0,1.0])
    # plt.xlabel('T bias [K]')
    #
    # ###-------------------------
    # plt.subplot(122)
    # ax1 = plt.gca()
    #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
    # #     np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
    # #     color = 'blue', alpha = 0.05)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
    #     np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
    #     np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
    # #     np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
    # #     color = 'navajowhite', alpha = 0.35)
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    #
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    # plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'purple', label = label4, zorder = 2)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
    # #     np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
    # #     color = 'blue', alpha = 0.05)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'darkblue', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
    # #     np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
    # #     color = 'lightblue', alpha = 0.15)
    # # plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'steelblue', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'steelblue', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
    # #     np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
    # #     color = 'mediumaquamarine', alpha = 0.15)
    # # plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # #
    # # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
    # #     np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
    # #     color = 'navajowhite', alpha = 0.35)
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
    # #     '--', color = 'gold', linewidth = 0.5)
    # #
    # # plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    # # plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # # plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4, zorder = 2)
    # # plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    # # plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)
    #
    # # plt.legend()
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax1.set_yticks(axminor, minor = True)
    # ax1.grid(which = 'major', alpha = 0.5)
    # # plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    # plt.grid('on')
    #
    # ###-------------------------
    # fileout = '../FIGS/CASIM/MedianProfiles_TandSpHum_casim-100_casim-aeroprof_wholeDrift.svg'
    # # plt.savefig(fileout)
    # plt.show()


    # ##################################################
    # ##################################################
    # #### 	P3 and P4
    # ##################################################
    # ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=LARGE_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.figure(figsize=(6.6,10))
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
    #         hspace = 0.22, wspace = 0.19)
    #
    # plt.subplot(221)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-5.0,2.5])
    # plt.xlabel('T bias [K]')
    #
    # plt.subplot(222)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    #
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])
    #
    # plt.subplot(223)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # plt.grid('on')
    # plt.xlim([-5.0,2.5])
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlabel('T bias [K]')
    #
    # plt.subplot(224)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # # plt.grid('on')
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])
    #
    # fileout = '../FIGS/CASIM/Temp-SpHumMedianProfiles_casim-100_casim-aeroprof_periodSelection-p3-p4_wSTDEV_newColours.svg'
    # # plt.savefig(fileout)
    # plt.show()
    #
    # ##################################################
    # ##################################################
    # #### 	P5 and P6
    # ##################################################
    # ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=LARGE_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.figure(figsize=(6.6,10))
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
    #         hspace = 0.22, wspace = 0.19)
    #
    # plt.subplot(221)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-5.0,2.5])
    # plt.xlabel('T bias [K]')
    #
    # plt.subplot(222)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    #
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])
    #
    # plt.subplot(223)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # plt.grid('on')
    # plt.xlim([-5.0,2.5])
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlabel('T bias [K]')
    #
    # plt.subplot(224)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1),
    #     color = 'violet', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'purple', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'purple', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    #
    # # plt.grid('on')
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])
    #
    # fileout = '../FIGS/CASIM/Temp-SpHumMedianProfiles_casim-100_casim-aeroprof_periodSelection-p5-p6_wSTDEV_newColours.svg'
    # # plt.savefig(fileout)
    # plt.show()



def plot_CASIM_QliqTimeseries(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Qliq timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(8,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    #### set flagged um_data to nans
    data2['qliq'][data2['qliq'] < 0] = np.nan
    data3['qliq'][data3['qliq'] < 0] = np.nan

    #### set up colourmaps to grey out zeros on figures
    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    plt.subplot(211)
    plt.pcolor(data2['time'], data2['height'][:], np.transpose(data2['qliq'])*1e3,
        vmin = 0, vmax = 0.35,
        cmap = mpl_cm.Blues)
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.title(label2 + ', Q_liq [g/kg]')
    plt.colorbar()

    plt.subplot(212)
    plt.pcolor(data3['time'], data3['height'][:], np.transpose(data3['qliq'])*1e3,
        vmin = 0, vmax = 0.35,
        cmap = mpl_cm.Blues)
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.xlabel('DOY')
    plt.title(label3 + ', Q_liq [g/kg]')
    plt.colorbar()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = '../FIGS/CASIM/CASIM-100_CASIM-AeroProf_QliqTimeseries_229-257DOY.png'
    # plt.savefig(fileout)
    plt.show()

def plot_CASIM_NiceTimeseries(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Ndrop timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(8,7))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    #### set flagged um_data to nans
    data2['qnice'][data2['qnice'] < 0] = np.nan
    data3['qnice'][data3['qnice'] < 0] = np.nan

    #### set up colourmaps to grey out zeros on figures
    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    plt.subplot(211)
    plt.pcolor(data2['time'], data2['height'][:], np.transpose(data2['qnice'])/1e6,
        vmin = 0, vmax = 0.1,
        cmap = mpl_cm.Blues)
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.title(label2 + ', N_ice [/kg]')
    plt.colorbar()

    plt.subplot(212)
    plt.pcolor(data3['time'], data3['height'][:], np.transpose(data3['qnice'])/1e6,
        vmin = 0, vmax = 0.1,
        cmap = mpl_cm.Blues)
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.xlabel('DOY')
    plt.title(label3 + ', N_ice [/kg]')
    plt.colorbar()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = '../FIGS/CASIM/CASIM-100_CASIM-AeroProf_NiceTimeseries_229-257DOY.png'
    plt.savefig(fileout)
    plt.show()

def plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined 1d timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(15,10))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_temp = obs['obs_temp'].variables['time'][:]
    time_temp = calcTime_Mat2DOY(datenums_temp)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))
    t1 = 11
    t3 = 45

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.subplot(3,2,1)

    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]

    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('CRF [W/m2]')
    # plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,2)
    ax1 = plt.gca()
    ax1.plot(time_tice,obs['obs_temp'].variables['Tice'][:] + 273.16, color = 'black', label = 'obs: ice')
    ax1.plot(data1['time'], data1['temp_1.5m'].data, color = 'darkblue', label = '1.5m')
    ax1.plot(data2['time'], data2['temp_1.5m'].data, color = 'purple')#, label = '2m')
    if ifs_flag == True:
        ax1.plot(data3['time'], data3['sfc_temp_2m'].data, color = 'gold', label = '2m')
    else:
        ax1.plot(data3['time'], data3['temp_1.5m'].data, color = 'gold')#, label = '2m')
    plt.title('near-sfc_temperature [K]')
    plt.legend()
    if month_flag == 8:  ax1.set_xlim([13.0, 31.0])
    if month_flag == 9:  ax1.set_xlim([1.0, 15.0])
    if month_flag == -1: ax1.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,3)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Observations')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    plt.legend()
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,4)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'purple')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    plt.title('surface_net_LW_radiation [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'purple')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1], obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1], 'k.')
    plt.plot(obs['ice_station_fluxes']['mday'],obs['ice_station_fluxes']['tafluxB'], 'r.')
    plt.ylim([-30, 30])
    plt.title('sensible_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,6)
    ax = plt.gca()
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'purple')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'gold')# * -1.0)
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1], obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1], 'k.')
    plt.title('latent_heat_flux [W/m2]')
    if month_flag == 8: ax.set_xlim([13.0, 31.0])
    if month_flag == 9: ax.set_xlim([1.0, 15.0])
    if month_flag == -1: ax.set_xlim([doy[0],doy[-1]])

    plt.subplot(3,2,5)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')
    plt.subplot(3,2,6)
    if month_flag == 8: plt.xlabel('Day of month [Aug]')
    if month_flag == 9: plt.xlabel('Day of month [Sep]')
    if month_flag == -1: plt.xlabel('Day of year')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/' + out_dir2[:9] + '_oden_metum_ra2t_ifs_TSa.svg'
    plt.savefig(fileout, dpi=300)
    plt.show()

def plot_Cv_RA2T(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, out_dir4, obs, doy, label1, label2, label3, label4): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    from physFuncts import calcAirDensity
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Cv for RA2T runs:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(11.5, 9))
    # plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
    #         hspace = 0.4, wspace = 0.13)

    ########            Cloud fraction
    ########
    #### set up colourmaps to grey out zeros on figures
    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr
    newcmp = ListedColormap(newcolors)

    #### load temp obs Cv data
    obs_data = np.load('../Cloudnet/obs_Cv_um-qf30.npy',allow_pickle=True, encoding = 'latin1').item()
    obs_data['Cv'][obs_data['Cv'] < 0.0] = 0.0

    period = np.where(np.logical_and(obs_data['time'] >= doy[0], obs_data['time'] < doy[-1]))

    ### build nanmask
    # nanmask = np.zeros([np.size(np.squeeze(obs_data['time'][period])), np.size(obs_data['Cv'],1)])
    # nanindex = np.zeros([np.size(np.squeeze(obs_data['time'][period]))])
    # print(nanindex.shape)
    # for i in range(len(np.squeeze(obs_data['time'][period]))):
    #     if np.isnan(np.nanmean(obs_data['Cv'][period[0][i],:])):       ## if there are only nans in the time profile
    #         nanmask[i,:] = 1.0
    #         nanmask[i-1,:] = 1.0
    #         nanmask[i+1,:] = 1.0
    #         nanindex[i] = 1
    # nanmask[nanmask == 0.0] = np.nan
    # nanind = np.where(nanindex == 1)
    #
    # obs_data['Cv'][nanind, :] = np.nan
    # data1['cloud_fraction'][nanind, :] = np.nan
    # data2['cloud_fraction'][nanind, :] = np.nan
    # data3['cloud_fraction'][nanind, :] = np.nan

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    ax0  = fig.add_axes([0.06,0.8, 0.48,0.17])   # left, bottom, width, height
    ax0 = plt.gca()
    # ax0.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.nanmean(obs_data['height'],0), np.transpose(obs_data['Cv']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    # plt.title('Obs-IFSgrid')# Cloud fraction')
    plt.xlim([doy[0], doy[-1]])
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,3.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax0.set_yticklabels([0,3,6,9])
    plt.xticks([244,245,246,247,248,249])
    ax0.set_xticklabels(['1 Sep','2 Sep','3 Sep','4 Sep','5 Sep','6 Sep'])
    axx = ax0.twinx()
    axx.set_ylabel('Obs-UMgrid', fontsize = 12, rotation = 270, labelpad = 13)
    axx.set_yticks([])
    cbaxes = fig.add_axes([0.59, 0.31, 0.015, 0.42])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'vertical')
    plt.xlabel('C$_{V}$', rotation = 0, labelpad = 15)
    cbaxes.xaxis.set_label_position('top')

    ax1  = fig.add_axes([0.06,0.56,0.48,0.17])   # left, bottom, width, height
    ax1 = plt.gca()
    plt.contourf(data2['time'], data2['height'][:], np.transpose(data2['cloud_fraction']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    # plt.title(label2)# + ' Cloud fraction')
    plt.xlim([doy[0], doy[-1]])
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,3,6,9])
    plt.xticks([244,245,246,247,248,249])
    ax1.set_xticklabels(['1 Sep','2 Sep','3 Sep','4 Sep','5 Sep','6 Sep'])
    axx = ax1.twinx()
    axx.set_ylabel('UM_RA2T \n noTurbMP',fontsize = 12,  rotation = 270, labelpad = 28)
    axx.set_yticks([])

    ax2  = fig.add_axes([0.06,0.32,0.48,0.17])   # left, bottom, width, height
    ax2 = plt.gca()
    plt.contourf(data4['time'], data4['height'][:], np.transpose(data4['cloud_fraction']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    # plt.title(label4[:-4])# + ' Cloud fraction')
    plt.xlim([doy[0], doy[-1]])
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax2.set_yticklabels([0,3,6,9])
    plt.xticks([244,245,246,247,248,249])
    ax2.set_xticklabels(['1 Sep','2 Sep','3 Sep','4 Sep','5 Sep','6 Sep'])
    axx = ax2.twinx()
    axx.set_ylabel(label4[:-4], fontsize = 12, rotation = 270, labelpad = 13)
    axx.set_yticks([])

    ax3 = fig.add_axes([0.06,0.08,0.48,0.17])   # left, bottom, width, height
    ax3 = plt.gca()
    plt.contourf(data1['time'], data1['height'][:], np.transpose(data1['cloud_fraction']),
        [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    # plt.title(label1)# + ' Cloud fraction')
    plt.xlabel('Date')
    # plt.colorbar()
    plt.xlim([doy[0], doy[-1]])
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax3.set_yticklabels([0,3,6,9])
    plt.xticks([244,245,246,247,248,249])
    ax3.set_xticklabels(['1 Sep','2 Sep','3 Sep','4 Sep','5 Sep','6 Sep'])
    axx = ax3.twinx()
    axx.set_ylabel(label1, fontsize = 12, rotation = 270, labelpad = 13)
    axx.set_yticks([])

    ax4  = fig.add_axes([0.72,0.22,0.25,0.5])   # left, bottom, width, height
    ax4 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][period,:]),0),np.nanmean(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs')
    ax4.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(np.squeeze(obs_data['Cv'][period,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][period,:]),0),
        np.nanmean(np.squeeze(obs_data['Cv'][period,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][period,:]),0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][period,:]),0) - np.nanstd(np.squeeze(obs_data['Cv'][period,:]),0), np.nanmean(obs_data['height'],0),
        'k--', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['Cv'][period,:]),0) + np.nanstd(np.squeeze(obs_data['Cv'][period,:]),0), np.nanmean(obs_data['height'],0),
        'k--', linewidth = 0.5)

    plt.plot(np.nanmean(data2['cloud_fraction'],0),data2['height'], color = 'darkorange', linewidth = 3, label = label2)
    ax4.fill_betweenx(data2['height'],np.nanmean(data2['cloud_fraction'],0) - np.nanstd(data2['cloud_fraction'],0),
        np.nanmean(data2['cloud_fraction'],0) + np.nanstd(data2['cloud_fraction'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(data2['cloud_fraction'],0) - np.nanstd(data2['cloud_fraction'],0), data2['height'],
        '--', color = 'darkorange', linewidth = 0.5)
    plt.plot(np.nanmean(data2['cloud_fraction'],0) + np.nanstd(data2['cloud_fraction'],0), data2['height'],
        '--', color = 'darkorange', linewidth = 0.5)

    plt.plot(np.nanmean(data4['cloud_fraction'],0),data4['height'], color = 'steelblue', linewidth = 3, label = label4[:-4])
    ax4.fill_betweenx(data4['height'],np.nanmean(data4['cloud_fraction'],0) - np.nanstd(data4['cloud_fraction'],0),
        np.nanmean(data4['cloud_fraction'],0) + np.nanstd(data4['cloud_fraction'],0), color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(data4['cloud_fraction'],0) - np.nanstd(data4['cloud_fraction'],0), data4['height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(data4['cloud_fraction'],0) + np.nanstd(data4['cloud_fraction'],0), data4['height'],
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(data1['cloud_fraction'],0),data1['height'], color = 'darkblue', linewidth = 3, label = label1)
    ax4.fill_betweenx(data1['height'],np.nanmean(data1['cloud_fraction'],0) - np.nanstd(data1['cloud_fraction'],0),
        np.nanmean(data1['cloud_fraction'],0) + np.nanstd(data1['cloud_fraction'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(data1['cloud_fraction'],0) - np.nanstd(data1['cloud_fraction'],0), data1['height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(data1['cloud_fraction'],0) + np.nanstd(data1['cloud_fraction'],0), data1['height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.xlabel('C$_{V}$')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax4.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax4.set_yticks(axminor, minor = True)
    ax4.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,1])
    plt.legend(bbox_to_anchor=(0.0, 1.02, 1., .102), loc=3, ncol=1)

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = '../Cloudnet/FIGS/CvTimeseries_Obs-all_' + label4 + '_' + label2 + '_' + label1 + '_1Sep-5Sep_Dates_FixedRA2T.svg'
    plt.savefig(fileout)
    # plt.show()
    plt.close()

def plot_CWC_RA2T(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    from physFuncts import calcAirDensity
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting cloud water contents for RA2T runs:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(13, 10))
    # plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
    #         hspace = 0.4, wspace = 0.13)

    ### choose var to plot [qliq / qice]
    var = 'qliq'

    #### set flagged data to zero
    data1[var][data1[var] < 0] = 0.0
    data2[var][data2[var] < 0] = 0.0
    # data3[var][data3[var] < 0] = 0.0
    data4[var][data4[var] < 0] = 0.0

    #### set nans to zero
    data1[var][np.isnan(data1[var])] = 0.0
    data2[var][np.isnan(data2[var])] = 0.0
    # data3[var][np.isnan(data3[var])] = 0.0
    data4[var][np.isnan(data4[var])] = 0.0

    ########
    ######## set up colourmaps to grey out small values on figures
    ########
    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr
    newcmp = ListedColormap(newcolors)

    ### define contour levels dependent on chosen var
    ###             define plot labels
    ###             define model parameter names from Cloudnet
    if var == 'qliq':
        crange = [0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
        lab = 'LWMR [g/kg]'
        modname = 'model_lwc'
        obsname = 'lwc'
    elif var == 'qice':
        crange = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
        lab = 'IWMR [g/kg]'
        modname = 'model_iwc_filtered'
        if out_dir3 == 'OUT_25H/':
            modname3 = 'model_snow_iwc_filtered'
        else:
            modname3 = modname
        obsname = 'iwc'

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    ax1  = fig.add_axes([0.08,0.7,0.58,0.22])   # left, bottom, width, height
    ax1 = plt.gca()
    plt.contourf(data1['time'], data1['height'][:], np.transpose(data1[var])*1e3,
        crange,
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.title(label1 + ' ' + lab)
    plt.colorbar()
    plt.xlim([doy[0], doy[-1]])

    ax2  = fig.add_axes([0.08,0.4,0.58,0.22])   # left, bottom, width, height
    ax2 = plt.gca()
    plt.contourf(data2['time'], data2['height'][:], np.transpose(data2[var])*1e3,
        crange,
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.title(label2 + ' ' + lab)
    plt.colorbar()
    plt.xlim([doy[0], doy[-1]])

    ax3 = fig.add_axes([0.08,0.1,0.58,0.22])   # left, bottom, width, height
    ax3 = plt.gca()
    # if out_dir3 == 'OUT_25H/':
    #     continue
    #     # ifs_data = np.load('../Cloudnet/ifs_Cv_data.npy').item()
    #     # plt.contourf(ifs_data['time'], np.nanmean(ifs_data['height'],0), np.transpose(ifs_data[modname3]),
    #     #     crange,
    #     #     # vmin = 0, vmax = 150,
    #     #     cmap = newcmp
    #     #     )
    # else:
    plt.contourf(data4['time'], data4['height'][:], np.transpose(data4[var])*1e3,
        crange,
        # vmin = 0, vmax = 150,
        cmap = newcmp
        )
    plt.ylabel('Height [m]')
    plt.ylim([0,8000])
    plt.title(label3 + ' ' + lab)
    plt.xlabel('Day of Year')
    plt.colorbar()
    plt.xlim([doy[0], doy[-1]])

    #### load temp obs Cv data
    obs_data = np.load('../Cloudnet/working_obs_data.npy').item()

    #### only include vals >1e-6kg/m3
    obs_data[obsname][obs_data[obsname] < 1e-6] = np.nan
    data1[var][data1[var] < 1e-6] = np.nan
    data2[var][data2[var] < 1e-6] = np.nan
    # data3[var][data3[var] < 1e-6] = np.nan
    data4[var][data4[var] < 1e-6] = np.nan

    ax4  = fig.add_axes([0.72,0.25,0.25,0.5])   # left, bottom, width, height
    ax4 = plt.gca()
    plt.plot(np.nanmean(obs_data[obsname]*1e3,0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs [g/m3]')
    ax4.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[obsname]*1e3,0) - np.nanstd(obs_data[obsname]*1e3,0),
        np.nanmean(obs_data[obsname]*1e3,0) + np.nanstd(obs_data[obsname]*1e3,0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(data1[var]*1e3,0),data1['height'], color = 'darkblue', linewidth = 3, label = label1)
    ax4.fill_betweenx(data1['height'],np.nanmean(data1[var]*1e3,0) - np.nanstd(data1[var]*1e3,0),
        np.nanmean(data1[var]*1e3,0) + np.nanstd(data1[var]*1e3,0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(data2[var]*1e3,0),data2['height'], color = 'purple', linewidth = 3, label = label2)
    ax4.fill_betweenx(data2['height'],np.nanmean(data2[var]*1e3,0) - np.nanstd(data2[var]*1e3,0),
        np.nanmean(data2[var]*1e3,0) + np.nanstd(data2[var]*1e3,0), color = 'thistle', alpha = 0.35)
    # if out_dir3 == 'OUT_25H/':
    #     continue
    #     # plt.plot(np.nanmean(ifs_data[modname3],0),np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS [g/m3]')
    #     # ax4.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data[modname3],0) - np.nanstd(ifs_data[modname3],0),
    #     #     np.nanmean(ifs_data[modname3],0) + np.nanstd(ifs_data[modname3],0), color = 'navajowhite', alpha = 0.15)
    # else:
    plt.plot(np.nanmean(data4[var]*1e3,0),data4['height'], color = 'mediumseagreen', linewidth = 3, label = label4)
    ax4.fill_betweenx(data4['height'],np.nanmean(data4[var]*1e3,0) - np.nanstd(data4[var]*1e3,0),
        np.nanmean(data4[var]*1e3,0) + np.nanstd(data4[var]*1e3,0), color = 'mediumaquamarine', alpha = 0.15)
    plt.xlabel(lab)
    plt.ylabel('Height [m]')
    plt.ylim([0,10000])
    if var == 'qliq': plt.xlim([0,0.15])
    if var == 'qice': plt.xlim([0,0.04])
    plt.legend()


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = '../Cloudnet/FIGS/' + var + 'Timeseries_Obs-all_' + label3 + '_' + label2 + '_' + label1 + '_14Aug-14Sep.png'
    # plt.savefig(fileout)
    plt.show()

def plot_line_subSect(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined timeseries and PDFs of radiation terms:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### add override for data2 to allow 24h data to be used for testing purposes
    if out_dir2[-4:] == '24h/':
        data2['surface_net_LW_radiation'][data2['surface_net_LW_radiation'] == 0] = np.nan
        data2['surface_net_SW_radiation'][data2['surface_net_SW_radiation'] == 0] = np.nan

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(12,10))

    ax  = fig.add_axes([0.07,0.7,0.55,0.22])   # left, bottom, width, height
    netLW = obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
    netSW = obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]
    ax = plt.gca()
    yA = [-65, 85]
    # plt.plot([240.0,240.0],[yA[0],yA[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice, netLW + netSW, color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('Net Radiation [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.legend(bbox_to_anchor=(-0.11, 0.65, 1., .102), loc=4, ncol=2)
    plt.ylim([-60,80])

    ax  = fig.add_axes([0.07,0.4,0.55,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    # plt.plot([240.0,240.0],[yB[0],yB[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['SWdice'][:] - obs['obs_temp'].variables['SWuice'][:]), color = 'black', label = 'Ice_station')
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'purple', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.title('surface_net_SW_radiation [W/m2]')
    # plt.legend()
    ax.set_xlim([doy[0],doy[-1]])
    plt.ylim([-3,120])

    ax  = fig.add_axes([0.07,0.1,0.55,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yC = [-90, 10]
    # plt.plot([240.0,240.0],[yC[0],yC[-1]],'--', color='red')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(time_radice,(obs['obs_temp'].variables['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]), color = 'black', label = 'obs: ice')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'purple')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    plt.title('surface_net_LW_radiation [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')
    plt.ylim([-90,5])

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)

    #### only compare obs over model dates available:
    ####        all model data share a timestamp

    subSect = np.where(np.logical_and(time_radice > data1['time_hrly'][0], time_radice <= data1['time_hrly'][-1]))

    sw1 = data1['surface_net_SW_radiation'][data1['hrly_flag']]
    lw1 = data1['surface_net_LW_radiation'][data1['hrly_flag']]
    sw2 = data2['surface_net_SW_radiation'][data2['hrly_flag']]
    lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    if ifs_flag == True:
        sw3 = data3['sfc_net_sw'][data3['hrly_flag']]
        lw3 = data3['sfc_net_lw'][data3['hrly_flag']]
    else:
        sw3 = data3['surface_net_SW_radiation'][data3['hrly_flag']]
        lw3 = data3['surface_net_LW_radiation'][data3['hrly_flag']]

    ax  = fig.add_axes([0.7,0.7,0.25,0.22])   # left, bottom, width, height
    yDmax = 0.12
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1 + lw1
    sns.distplot(crf1, hist=False, color="darkblue", kde_kws={"shade": True})
    crf3 = sw3 + lw3
    sns.distplot(crf3, hist=False, color="gold", kde_kws={"shade": True})
    crf2 = sw2 + lw2
    sns.distplot(crf2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netLW[subSect] + netSW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(55,0.07), xytext=(55,0.07), fontsize = 14)
    plt.xlabel('CRF [W/m2]')
    plt.xlim([-50,80])
    plt.ylim([0,yDmax])

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.4,0.25,0.22])   # left, bottom, width, height
    yEmax = 0.16
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(sw1, hist=False, color="darkblue", kde_kws={"shade": True})
    sns.distplot(sw3, hist=False, color="gold", kde_kws={"shade": True})
    sns.distplot(sw2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netSW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(87,0.07), xytext=(87,0.07), fontsize = 14)
    # plt.legend()
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('$SW_{net,surf}$ [W/m2]')

    # plt.subplot(212)
    ax  = fig.add_axes([0.7,0.1,0.25,0.22])   # left, bottom, width, height
    yFmax = 0.13
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    sns.distplot(lw1, hist=False, color="darkblue", kde_kws={"shade": True})
    sns.distplot(lw3, hist=False, color="gold", kde_kws={"shade": True})
    sns.distplot(lw2, hist=False, color="purple", kde_kws={"shade": True})
    sns.distplot(netLW[subSect], hist=False, color="black")
    # plt.title('Melt')
    # plt.annotate('Melt', xy=(0,0.14), xytext=(0,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('$LW_{net,surf}$ [W/m2]')


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/CRF_netSW_netLW_line+PDFS_oden_iceStation_' + label1[3:] + '_' + label3 + '_' + label2[3:] + '.png'
    plt.savefig(fileout)
    plt.show()

def plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting timeseries of turbulent fluxes:')
    print ('')

    ##################################################
    ##################################################
    #### 	SET AXES PROPERTIES
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 18
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### set diagnostic naming flags for if IFS being used
    ### -------------------------------
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### -------------------------------
    ### for reference in figures
    ### -------------------------------
    zeros = np.zeros(len(data2['time']))

    ### -------------------------------
    ### change matlab time for obs
    ### -------------------------------
    time_iceStation = calcTime_Mat2DOY(np.squeeze(obs['ice_station_fluxes']['mday'][:]))

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12))

    ax  = fig.add_axes([0.07,0.55,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1],
        obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    plt.plot(time_iceStation[np.squeeze(obs['ice_station_fluxes']['taflag'][:]==1)],
        np.squeeze(obs['ice_station_fluxes']['taflux'][obs['ice_station_fluxes']['taflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 7, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['sensible_heat_flux'].data, color = 'darkblue', label = label1)
    plt.plot(data2['time'], data2['sensible_heat_flux'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_sens_heat_flx'].data * -1.0, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['sensible_heat_flux'].data, color = 'gold')
    plt.ylim([-20, 40])
    plt.title('sensible_heat_flux [W/m2]')
    ax.set_xlim([doy[0],doy[-1]])

    ax  = fig.add_axes([0.07,0.1,0.56,0.35])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['foremast'].variables['doy'][obs['foremast'].variables['rflag'][:] == 1],
        obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1],
        'kv', markersize = 5, label = 'Foremast')
    # index = np.logical_and(obs['ice_station_fluxes']['lrflux']>=-30, obs['ice_station_fluxes']['lrflux']<=70)
    plt.plot(time_iceStation[np.squeeze(obs['ice_station_fluxes']['lrflag'][:]==1)],
        np.squeeze(obs['ice_station_fluxes']['lrflux'][obs['ice_station_fluxes']['lrflag'][:] == 1]),
        '^', color = 'darkgrey', markersize = 7, markeredgecolor = 'grey', label = 'Ice_station')
    plt.plot(data1['time'], data1['latent_heat_flux'].data, color = 'darkblue')
    plt.plot(data2['time'], data2['latent_heat_flux'].data, color = 'mediumseagreen')# * -1.0)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_down_lat_heat_flx'].data * -1.0, color = 'gold')
    else:
        plt.plot(data3['time'], data3['latent_heat_flux'].data, color = 'gold')# * -1.0)
    plt.title('latent_heat_flux [W/m2]')
    plt.ylim([-20, 60])
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)
    ax  = fig.add_axes([0.7,0.55,0.27,0.35])   # left, bottom, width, height
    # zerosC = np.zeros(len(data2['time']))
    yCmax = 0.16
    plt.plot([0,0],[0,yCmax],'--', color='lightgrey')
    ##---
    shf1 = data1['sensible_heat_flux'][data1['hrly_flag']].data
    indextaum = np.logical_and(shf1 >= -50, shf1 <= 50)
    sns.distplot(shf1[indextaum], hist=False, color="darkblue", kde_kws={"shade": True}, label = label1)
    ##---
    shf3 = data3['sfc_down_sens_heat_flx'][data3['hrly_flag']].data * -1.0
    indextaifs = np.logical_and(shf3 >= -50, shf3 <= 50)
    sns.distplot(shf3[indextaifs], hist=False, color="gold", kde_kws={"shade": True}, label = label3)
    ##---
    shf2 = data2['sensible_heat_flux'][data2['hrly_flag']].data
    indextacasim = np.logical_and(shf2 >= -50, shf2 <= 50)
    sns.distplot(shf2[indextacasim], hist=False, color="mediumseagreen", kde_kws={"shade": True}, label = label2)
    ##---
    fmst_taflux = obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1]
    indextafmst = np.logical_and(fmst_taflux>=-50, fmst_taflux<=50)
    sns.distplot(fmst_taflux[indextafmst], hist=False, color="black", label = 'Foremast')#, kde_kws={"shade": True}, label = 'Foremast')
    ##---
    taflux = np.squeeze(obs['ice_station_fluxes']['taflux'][obs['ice_station_fluxes']['taflag'][:] == 1])
    indexta = np.logical_and(taflux>=-50, taflux<=50)
    sns.distplot(taflux[indexta], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3}, label = 'Ice_station')
    plt.title('sensible_heat_flux [W/m2]')
    plt.legend()
    plt.xlim([-20,50])
    plt.ylim([0,yCmax])

    ax  = fig.add_axes([0.7,0.1,0.27,0.35])   # left, bottom, width, height
    yDmax = 0.12
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    ##---
    lhf1 = data1['latent_heat_flux'][data1['hrly_flag']].data
    indexlrum = np.logical_and(lhf1 >= -50, lhf1 <= 50)
    sns.distplot(lhf1[indexlrum], hist=False, color="darkblue", kde_kws={"shade": True})
    ##---
    lhf3 = data3['sfc_down_lat_heat_flx'][data3['hrly_flag']].data * -1.0
    indexlrifs = np.logical_and(lhf3 >= -50, lhf3 <= 50)
    sns.distplot(lhf3[indexlrifs], hist=False, color="gold", kde_kws={"shade": True})
    ##---
    lhf2 = data2['latent_heat_flux'][data2['hrly_flag']].data
    indexlrcasim = np.logical_and(lhf2 >= -50, lhf2 <= 50)
    sns.distplot(lhf2[indexlrcasim], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    ##---
    fmst_lrflux = obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1]
    indexlrfmst = np.logical_and(fmst_lrflux>=-50, fmst_lrflux<=50)
    sns.distplot(fmst_lrflux[indexlrfmst], hist=False, color="black")#, kde_kws={"shade": True})
    ##---
    lrflux = np.squeeze(obs['ice_station_fluxes']['lrflux'][obs['ice_station_fluxes']['lrflag'][:] == 1])
    indexlr = np.logical_and(lrflux>=-50, lrflux<=50)
    sns.distplot(lrflux[indexlr], hist=False, color="grey", kde_kws={'linestyle':'--','linewidth':3})
    plt.title('latent_heat_flux [W/m2]')
    plt.xlim([-20,50])
    plt.ylim([0,yDmax])

    fileout = '../FIGS/comparisons/SHF_LHF_line+PDFS_oden_foremast+iceStationQC_metum_ifs_casim-100.svg'
    plt.savefig(fileout)
    plt.show()

def table_Fluxes(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('calculating fluxes for chosen period selection:')
    print ('')

    #################################################################
    ## sort out observations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    # datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation (ice station) on different timestep
    # time_radice_all = calcTime_Mat2DOY(datenums_radice)
    #
    # datenums_radship = obs['obs_temp'].variables['time2'][:] ### radiation (ship) on different timestep
    # time_radship_all = calcTime_Mat2DOY(datenums_radship)

    # datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    # time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ship_time = obs['foremast'].variables['doy'][obs['foremast'].variables['taflag'][:] == 1]
    ship_shf = obs['foremast'].variables['taflux'][obs['foremast'].variables['taflag'][:] == 1]
    ship_lhf = obs['foremast'].variables['rflux'][obs['foremast'].variables['rflag'][:] == 1]

    ### update these to probe data['fixed_radiation']['time'] arrays instead
    p3mod = np.where(np.logical_and(data1['time_hrly'] >= doy[0], data1['time_hrly'] < 230.0))
    p4mod = np.where(np.logical_and(data1['time_hrly'] >= 230.0, data1['time_hrly'] < 240.0))
    p5mod = np.where(np.logical_and(data1['time_hrly'] >= 240.0, data1['time_hrly'] < 247.0))
    p6mod = np.where(np.logical_and(data1['time_hrly'] >= 247.0, data1['time_hrly'] < 251.0))

    p3ship = np.where(np.logical_and(ship_time >= doy[0], ship_time < 230.0))
    p4ship = np.where(np.logical_and(ship_time >= 230.0, ship_time < 240.0))
    p5ship = np.where(np.logical_and(ship_time >= 240.0, ship_time < 247.0))
    p6ship = np.where(np.logical_and(ship_time >= 247.0, ship_time < 251.0))

    # plt.plot(data1['time'])

    #####-------------------------------
    #####   net LW
    #####-------------------------------
    print (str(np.round(np.nanmean(np.squeeze(data3['sfc_net_sw'][p3mod])),2)))

    print ()
    print ('LHF (p3):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_lhf[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p3mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p3mod]*-1.0)) - np.nanmean(np.squeeze(ship_lhf[p3ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_lhf[p3ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_lhf[p3ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_lhf[p3ship])),2)) + ')')

    print ()
    print ('SHF (p3):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_shf[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p3mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p3mod]*-1.0)) - np.nanmean(np.squeeze(ship_shf[p3ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_shf[p3ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_shf[p3ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p3mod])) - np.nanmean(np.squeeze(ship_shf[p3ship])),2)) + ')')

    print ()
    print ('*******************')
    print ()

    print ()
    print ('LHF (p4):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_lhf[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p4mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p4mod]*-1.0)) - np.nanmean(np.squeeze(ship_lhf[p4ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_lhf[p4ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_lhf[p4ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_lhf[p4ship])),2)) + ')')

    print ()
    print ('SHF (p4):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_shf[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p4mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p4mod]*-1.0)) - np.nanmean(np.squeeze(ship_shf[p4ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_shf[p4ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_shf[p4ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p4mod])) - np.nanmean(np.squeeze(ship_shf[p4ship])),2)) + ')')

    print ()
    print ('*******************')
    print ()

    print ()
    print ('LHF (p5):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_lhf[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p5mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p5mod]*-1.0)) - np.nanmean(np.squeeze(ship_lhf[p5ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_lhf[p5ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_lhf[p5ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_lhf[p5ship])),2)) + ')')

    print ()
    print ('SHF (p5):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_shf[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p5mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p5mod]*-1.0)) - np.nanmean(np.squeeze(ship_shf[p5ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_shf[p5ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_shf[p5ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p5mod])) - np.nanmean(np.squeeze(ship_shf[p5ship])),2)) + ')')

    print ()
    print ('*******************')
    print ()

    print ()
    print ('LHF (p6):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_lhf[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p6mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_lat_heat_flx'][p6mod]*-1.0)) - np.nanmean(np.squeeze(ship_lhf[p6ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['latent_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_lhf[p6ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['latent_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_lhf[p6ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['latent_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_lhf[p6ship])),2)) + ')')

    print ()
    print ('SHF (p6):')
    # print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(ship_shf[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p6mod])*-1.0),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['sfc_down_sens_heat_flx'][p6mod]*-1.0)) - np.nanmean(np.squeeze(ship_shf[p6ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['sensible_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_shf[p6ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['sensible_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_shf[p6ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['sensible_heat_flux'][p6mod])) - np.nanmean(np.squeeze(ship_shf[p6ship])),2)) + ')')

def plot_paperRadiation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined timeseries and PDFs of radiation terms:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    # datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    # time_radice_all = calcTime_Mat2DOY(datenums_radice)
    #
    # datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    # time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ##############################################################################
    ####        need the following if NOT using check_Radiation
    ##############################################################################
    # datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    # time_radice_all = calcTime_Mat2DOY(datenums_radice)
    # time_radice = time_radice_all[:-1:2]
    # lwd = obs['obs_temp'].variables['LWdice'][:]
    # lwu = obs['obs_temp'].variables['LWuice'][:]
    # lwdmean = np.nanmean(lwd[:-1].reshape(-1, 2), axis=1)
    # lwumean = np.nanmean(lwu[:-1].reshape(-1, 2), axis=1)
    # swd = obs['obs_temp'].variables['SWdice'][:]
    # swu = obs['obs_temp'].variables['SWuice'][:]
    # swdmean = np.nanmean(swd[:-1].reshape(-1, 2), axis=1)
    # swumean = np.nanmean(swu[:-1].reshape(-1, 2), axis=1)
    # netLW = lwdmean - lwumean
    # netSW = swdmean - swumean
    #
    # obs['fixed_radiation'] = {}
    # obs['fixed_radiation']['time_ice'] = time_radice
    # obs['fixed_radiation']['SWnet_ice'] = netSW
    # obs['fixed_radiation']['LWnet_ice'] = netLW


    #########-------------------------------------------------------------------------------------------
    ####       DEFINE PERIODS
    ####               all model data share a timestamp
    p3 = np.where(np.logical_and(data1['time_hrly'] >= doy[0], data1['time_hrly'] <= 230.0))
    p4 = np.where(np.logical_and(data1['time_hrly'] >= 230.0, data1['time_hrly'] <= 240.0))
    p5 = np.where(np.logical_and(data1['time_hrly'] >= 240.0, data1['time_hrly'] <= 247.0))
    p6 = np.where(np.logical_and(data1['time_hrly'] >= 247.0, data1['time_hrly'] <= 251.0))
    p7 = np.where(np.logical_and(data1['time_hrly'] >= 251.0, data1['time_hrly'] <= 255.5))
    p8 = np.where(data1['time_hrly'] >= 255.5)

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12))

    ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height

    ax = plt.gca()
    yB = [-10, 120]
    plt.plot([240.0,240.0],[yB[0],yB[-1]],'--', color='grey')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # plt.plot(data1['time_hrly'][:-3], data1['fixed_radiation']['SWnet'].data, color = 'darkblue', label = label1)
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data4['time'], data4['surface_net_SW_radiation'].data, color = 'steelblue', label = label4[:-4])
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['SWnet_ice'], color = 'grey', linewidth = 3, label = 'Ice_station')
    plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['SWnet_ship'], color = 'k', linewidth = 2, label = 'Ship')
    plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.ylim([-3,120])

    ax  = fig.add_axes([0.07,0.4,0.53,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yC = [-90, 10]
    plt.plot([240.0,240.0],[yC[0],yC[-1]],'--', color='grey')
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    plt.plot(data4['time'], data4['surface_net_LW_radiation'].data, color = 'steelblue')
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'mediumseagreen')
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['LWnet_ice'], color = 'grey', linewidth = 3)
    plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['LWnet_ship'], color = 'k')
    plt.ylabel('LW$_{net}$ [W m$^{-2}$]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.ylim([-90,5])

    ax  = fig.add_axes([0.07,0.1,0.53,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yA = [-65, 85]
    #### plot periods
    ##  melt/freeze threshold
    plt.plot([240.0,240.0],[yA[0],yA[-1]],'--', color='grey')
    ## p3
    plt.plot([data1['time_hrly'][p3[0][0]], data1['time_hrly'][p3[0][-1]]], [-50, -50], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p3[0][0]], -50, '.', color = 'r')
    plt.plot(data1['time_hrly'][p3[0][-1]], -50, '.', color = 'r')
    plt.annotate('P3', xy=(227.5,-45), xytext=(227.5,-45), fontsize = 12, color = 'r')
    plt.plot([data1['time_hrly'][p3[0][-1]], data1['time_hrly'][p3[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    ## p4
    plt.plot([data1['time_hrly'][p4[0][0]], data1['time_hrly'][p4[0][-1]]], [-50, -50], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p4[0][0]], -50, '.', color = 'r')
    plt.plot(data1['time_hrly'][p4[0][-1]], -50, '.', color = 'r')
    plt.annotate('P4', xy=(234.5,-45), xytext=(234.5,-45), fontsize = 12, color = 'r')
    # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [yA[0],yA[-1]], '--', color = 'r', linewidth = 1)
    ## p5
    plt.plot([data1['time_hrly'][p5[0][0]], data1['time_hrly'][p5[0][-1]]], [70, 70], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p5[0][0]], 70, '.', color = 'r')
    plt.plot(data1['time_hrly'][p5[0][-1]], 70, '.', color = 'r')
    plt.annotate('P5', xy=(243,59), xytext=(243,59), fontsize = 12, color = 'r')
    plt.plot([data1['time_hrly'][p5[0][-1]], data1['time_hrly'][p5[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    ## p6
    plt.plot([data1['time_hrly'][p6[0][0]], data1['time_hrly'][p6[0][-1]]], [70, 70], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p6[0][0]], 70, '.', color = 'r')
    plt.plot(data1['time_hrly'][p6[0][-1]], 70, '.', color = 'r')
    plt.annotate('P6', xy=(248.75,59), xytext=(248.75,59), fontsize = 12, color = 'r')
    plt.plot([data1['time_hrly'][p6[0][-1]], data1['time_hrly'][p6[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    ## p7
    plt.plot([data1['time_hrly'][p7[0][0]], data1['time_hrly'][p7[0][-1]]], [70, 70], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p7[0][0]], 70, '.', color = 'r')
    plt.plot(data1['time_hrly'][p7[0][-1]], 70, '.', color = 'r')
    plt.annotate('P7', xy=(253,59), xytext=(253,59), fontsize = 12, color = 'r')
    plt.plot([data1['time_hrly'][p7[0][-1]], data1['time_hrly'][p7[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    ## p8
    plt.plot([data1['time_hrly'][p8[0][0]], data1['time_hrly'][p8[0][-1]]], [70, 70], '-.', color = 'r')
    plt.plot(data1['time_hrly'][p8[0][0]], 70, '.', color = 'r')
    plt.plot(data1['time_hrly'][p8[0][-1]], 70, '.', color = 'r')
    plt.annotate('P8', xy=(256.5,59), xytext=(256.5,59), fontsize = 12, color = 'r')
    ##
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data4['time'], data4['surface_net_LW_radiation'].data + data4['surface_net_SW_radiation'].data, color = 'steelblue', label = label4[:-4])
    plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    if ifs_flag == True:
        plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    else:
        plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['LWnet_ice'] + obs['fixed_radiation']['SWnet_ice'], color = 'grey',  linewidth = 3, label = 'Ice_station')
    plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['LWnet_ship'] + obs['fixed_radiation']['SWnet_ship'], color = 'k', label = 'Ship')
    plt.ylabel('Net Radiation [W m$^{-2}$]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.ylim([-60,80])
    plt.xlabel('Date')

    ### -------------------------------
    ### Build figure (PDFs)
    ### -------------------------------
    # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # fig = plt.figure(figsize=(7,9))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.15)
    # plt.subplot(211)

    #### only compare with observations where we have data:
    ####        all model data share a timestamp
    melt = np.where(np.logical_and(data1['time_hrly'][:-3] >= obs['fixed_radiation']['time_ice'][0], data1['time_hrly'][:-3] < 240.0))
    freeze = np.where(data1['time_hrly'][:-3] >= 240.0)

    icemelt = np.where(obs['fixed_radiation']['time_ice'] < 240.0)
    icefreeze = np.where(obs['fixed_radiation']['time_ice'] >= 240.0)
    shipmelt = np.where(obs['fixed_radiation']['time_ship'] < 240.0)
    shipfreeze = np.where(obs['fixed_radiation']['time_ship'] >= 240.0)

    # sw1 = data1['surface_net_SW_radiation'][data1['hrly_flag']]
    # lw1 = data1['surface_net_LW_radiation'][data1['hrly_flag']]
    # sw2 = data2['surface_net_SW_radiation'][data2['hrly_flag']]
    # lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    # sw3 = data3['sfc_net_sw'][data3['hrly_flag']]
    # lw3 = data3['sfc_net_lw'][data3['hrly_flag']]
    # sw4 = data4['surface_net_SW_radiation'][data4['hrly_flag']]
    # lw4 = data4['surface_net_LW_radiation'][data4['hrly_flag']]
    sw1 = data1['fixed_radiation']['SWnet']
    lw1 = data1['fixed_radiation']['LWnet']
    sw2 = data2['fixed_radiation']['SWnet']
    lw2 = data2['fixed_radiation']['LWnet']
    sw3 = data3['fixed_radiation']['SWnet']
    lw3 = data3['fixed_radiation']['LWnet']
    sw4 = data4['fixed_radiation']['SWnet']
    lw4 = data4['fixed_radiation']['LWnet']

    ax  = fig.add_axes([0.64,0.7,0.15,0.22])   # left, bottom, width, height
    yEmax = 0.1
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    f = sns.distplot(sw1[melt], hist=False, color="darkblue", kde_kws={"shade": True})
    f = sns.distplot(sw4[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    f = sns.distplot(sw2[melt], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    f = sns.distplot(sw3[melt], hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipmelt], hist=False, color="black")
    plt.annotate('Melt', xy=(87,0.087), xytext=(87,0.087), fontsize = 14)
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('SW$_{net}$ [W m$^{-2}$]')
    # ## Obs_ship peak SWnet = 18.22
    # ## Obs_ice peak SWnet = 13.30
    # ## ECMWF_IFS peak SWnet = 21.61
    # ## UM_CASIM-100 peak SWnet = 26.92
    # ## UM_RA2T peak SWnet = 41.85
    # ## UM_RA2M peak SWnet = 40.41
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('UM_RA2M peak SWnet = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    ax  = fig.add_axes([0.64,0.4,0.15,0.22])   # left, bottom, width, height
    yFmax = 0.16
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    f = sns.distplot(lw1[melt], hist=False, color="darkblue", kde_kws={"shade": True})
    f = sns.distplot(lw4[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    f = sns.distplot(lw2[melt], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    f = sns.distplot(lw3[melt], hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['LWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['LWnet_ship'][shipmelt], hist=False, color="black")
    plt.annotate('Melt', xy=(0,0.14), xytext=(0,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('LW$_{net}$ [W m$^{-2}$]')
    # # UM_RA2M peak LWnet = -3.40
    # # UM_RA2T peak LWnet = -4.01
    # # UM_CASIM-100 peak LWnet = -2.05
    # # ECMWF_IFS peak LWnet = -3.68
    # # Obs_ice peak LWnet = -3.71
    # # Obs_ship peak LWnet = -2.12
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('Obs_ship peak LWnet = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    ax  = fig.add_axes([0.64,0.1,0.15,0.22])   # left, bottom, width, height
    yDmax = 0.08
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1[melt] + lw1[melt]
    f = sns.distplot(crf1, hist=False, color="darkblue", kde_kws={"shade": True})
    crf4 = sw4[melt] + lw4[melt]
    f = sns.distplot(crf4, hist=False, color="steelblue", kde_kws={"shade": True})
    crf2 = sw2[melt] + lw2[melt]
    f = sns.distplot(crf2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    crf3 = sw3[melt] + lw3[melt]
    f = sns.distplot(crf3, hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icemelt] + obs['fixed_radiation']['LWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipmelt] + obs['fixed_radiation']['LWnet_ship'][shipmelt], hist=False, color="black")
    plt.annotate('Melt', xy=(47,0.07), xytext=(47,0.07), fontsize = 14)
    plt.xlabel('Net Radiation [W m$^{-2}$]')
    plt.xlim([-80,80])
    plt.ylim([0,yDmax])
    # # UM_RA2M peak netRad = 33.57
    # # UM_RA2T peak netRad = 37.46
    # # UM_CASIM-100 peak netRad = 22.47
    # # ECMWF_IFS peak netRad = 15.56
    # # Obs_ice peak netRad = 8.31
    # # Obs_ship peak netRad = 16.94
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('UM_RA2T peak netRad = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    ax  = fig.add_axes([0.83,0.7,0.15,0.22])   # left, bottom, width, height
    yEmax = 0.1
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    f = sns.distplot(sw1[freeze], hist=False, color="darkblue", kde_kws={"shade": True})
    f = sns.distplot(sw4[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    f = sns.distplot(sw2[freeze], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    f = sns.distplot(sw3[freeze], hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipfreeze], hist=False, color="black")
    plt.annotate('Freeze', xy=(77,0.087), xytext=(77,0.087), fontsize = 14)
    plt.xlim([-10,110])
    plt.ylim([0,yEmax])
    plt.xlabel('SW$_{net}$ [W m$^{-2}$]')
    # # UM_RA2M peak SWnet = 25.25
    # # UM_RA2T peak SWnet = 26.57
    # # UM_CASIM-100 peak SWnet = 10.65
    # # ECMWF_IFS peak SWnet = 10.00
    # # Obs_ice peak SWnet = 8.67
    # # Obs_ship peak SWnet = 7.87
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('UM_RA2T peak SWnet = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    ax  = fig.add_axes([0.83,0.4,0.15,0.22])   # left, bottom, width, height
    yFmax = 0.16
    plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    f = sns.distplot(lw1[freeze], hist=False, color="darkblue", kde_kws={"shade": True})
    f = sns.distplot(lw4[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    f = sns.distplot(lw2[freeze], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    f = sns.distplot(lw3[freeze], hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['LWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['LWnet_ship'][shipfreeze], hist=False, color="black")
    plt.annotate('Freeze', xy=(-8,0.14), xytext=(-8,0.14), fontsize = 14)
    plt.xlim([-80,20])
    plt.ylim([0,yFmax])
    plt.xlabel('LW$_{net}$ [W m$^{-2}$]')
    # # UM_RA2M peak LWnet = -9.31
    # # UM_RA2T peak LWnet = -11.49
    # # UM_CASIM-100 peak LWnet = -5.65
    # # ECMWF_IFS peak LWnet = -6.83
    # # Obs_ice peak LWnet = -7.57
    # # Obs_ship peak LWnet = -6.50
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('UM_RA2T peak LWnet = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    # plt.subplot(212)
    ax  = fig.add_axes([0.83,0.1,0.15,0.22])   # left, bottom, width, height
    yDmax = 0.08
    plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    crf1 = sw1[freeze] + lw1[freeze]
    f = sns.distplot(crf1, hist=False, color="darkblue", kde_kws={"shade": True})
    crf4 = sw4[freeze] + lw4[freeze]
    f = sns.distplot(crf4, hist=False, color="steelblue", kde_kws={"shade": True})
    crf2 = sw2[freeze] + lw2[freeze]
    f = sns.distplot(crf2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    crf3 = sw3[freeze] + lw3[freeze]
    f = sns.distplot(crf3, hist=False, color="gold", kde_kws={"shade": True})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icefreeze] + obs['fixed_radiation']['LWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipfreeze] + obs['fixed_radiation']['LWnet_ship'][shipfreeze], hist=False, color="black")
    plt.annotate('Freeze', xy=(35,0.07), xytext=(35,0.07), fontsize = 14)
    plt.xlim([-80,80])
    plt.ylim([0,yDmax])
    plt.xlabel('Net Radiation [W m$^{-2}$]')
    # # UM_RA2M peak netRad = 14.91
    # # UM_RA2T peak netRad = 18.18
    # # UM_CASIM-100 peak netRad = 6.10
    # # ECMWF_IFS peak netRad = 5.62
    # # Obs_ice peak netRad = 0.98
    # # Obs_ship peak netRad = 1.73
    # y = f.lines[0].get_ydata()
    # x = f.lines[0].get_xdata()
    # maxid = np.argmax(y)
    # print ('UM_RA2T peak netRad = ' + ('%.2f' % x[maxid]))
    # plt.plot(x[maxid],y[maxid],'o')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/netSW_netLW_netRad_line+PDFS-gt230DOY_oden_iceStation_metum_ifs_casim-aeroprof_ra2t_splitSeason_fixLabels_newColours_Dates_wPeriods_wShipRadiation_fixedRA2T.svg'
    # plt.savefig(fileout)
    plt.show()

def table_Radiation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('calculating radiation terms for chosen period selection:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(9,10))
    # # plt.rc('figure',titlesize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.08,
    #         hspace = 0.4, wspace = 0.13)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    # datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation (ice station) on different timestep
    # time_radice_all = calcTime_Mat2DOY(datenums_radice)
    #
    # datenums_radship = obs['obs_temp'].variables['time2'][:] ### radiation (ship) on different timestep
    # time_radship_all = calcTime_Mat2DOY(datenums_radship)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    time_radice = obs['fixed_radiation']['time_ice'][:]
    time_radship = obs['fixed_radiation']['time_ship'][:]
    lwdice = obs['fixed_radiation']['LWd_ice'][:]
    lwuice = obs['fixed_radiation']['LWu_ice'][:]

    netSW = obs['fixed_radiation']['SWnet_ship'][:]
    netLW = obs['fixed_radiation']['LWnet_ship'][:]

    lwdmean = obs['fixed_radiation']['LWd_ship'][:]
    swdmean = obs['fixed_radiation']['SWd_ship'][:]

    netSWice = obs['fixed_radiation']['SWnet_ice'][:]
    netLWice = obs['fixed_radiation']['LWnet_ice'][:]

    lwdmeanice = obs['fixed_radiation']['LWd_ice'][:]
    swdmeanice = obs['fixed_radiation']['SWd_ice'][:]

    ### update these to probe data['fixed_radiation']['time'] arrays instead
    p3mod = np.where(np.logical_and(data1['fixed_radiation']['time'] >= doy[0], data1['fixed_radiation']['time'] < 230.0))
    p4mod = np.where(np.logical_and(data1['fixed_radiation']['time'] >= 230.0, data1['fixed_radiation']['time'] < 240.0))
    p5mod = np.where(np.logical_and(data1['fixed_radiation']['time'] >= 240.0, data1['fixed_radiation']['time'] < 247.0))
    p6mod = np.where(np.logical_and(data1['fixed_radiation']['time'] >= 247.0, data1['fixed_radiation']['time'] < 251.0))
    # p3mod = np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] < 230.0))
    # p4mod = np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] < 240.0))
    # p5mod = np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] < 247.0))
    # p6mod = np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] < 251.0))

    p3ship = np.where(np.logical_and(time_radship >= doy[0], time_radship < 230.0))
    p4ship = np.where(np.logical_and(time_radship >= 230.0, time_radship < 240.0))
    p5ship = np.where(np.logical_and(time_radship >= 240.0, time_radship < 247.0))
    p6ship = np.where(np.logical_and(time_radship >= 247.0, time_radship < 251.0))

    p3ice = np.where(np.logical_and(time_radice >= doy[0], time_radice < 230.0))
    p4ice = np.where(np.logical_and(time_radice >= 230.0, time_radice < 240.0))
    p5ice = np.where(np.logical_and(time_radice >= 240.0, time_radice < 247.0))
    p6ice = np.where(np.logical_and(time_radice >= 247.0, time_radice < 251.0))

    # plt.plot(data1['time'])

    #####-------------------------------
    #####   net LW
    #####-------------------------------
    print (str(np.round(np.nanmean(np.squeeze(data3['sfc_net_sw'][p3mod])),2)))

    print ()
    print ('net SW (p3):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netSW[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSW[p3ship])),2)) + '')
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSWice[p3ice])),2)) + '')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSW[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSWice[p3ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSW[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSWice[p3ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSW[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p3mod])) - np.nanmean(np.squeeze(netSWice[p3ice])),2)))

    print ()
    print ('Downwelling SW (p3):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(swdmean[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmean[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmean[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmean[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p3mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmean[p3ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p3mod])) - np.nanmean(np.squeeze(swdmeanice[p3ice])),2)))

    print ()
    print ('net LW (p3):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netLWice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netLW[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p3mod])) - np.nanmean(np.squeeze(netLW[p3ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p3mod])) - np.nanmean(np.squeeze(netLW[p3ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p3mod])) - np.nanmean(np.squeeze(netLW[p3ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p3mod])) - np.nanmean(np.squeeze(netLW[p3ship])),2)) + ')')

    print ()
    print ('Downwelling LW (p3):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(lwdmeanice[p3ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(lwdmean[p3ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p3mod])) - np.nanmean(np.squeeze(lwdmean[p3ship])))) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p3mod])) - np.nanmean(np.squeeze(lwdmean[p3ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p3mod])) - np.nanmean(np.squeeze(lwdmean[p3ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p3mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p3mod])) - np.nanmean(np.squeeze(lwdmean[p3ship])),2)) + ')')


    print ()
    print ('*******************')
    print ()

    print ()
    print ('net SW (p4):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netSW[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSW[p4ship])),2)) + '')
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSWice[p4ice])),2)) + '')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSW[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSWice[p4ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSW[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSWice[p4ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSW[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p4mod])) - np.nanmean(np.squeeze(netSWice[p4ice])),2)))

    print ()
    print ('Downwelling SW (p4):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(swdmean[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmean[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmean[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmean[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p4mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmean[p4ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p4mod])) - np.nanmean(np.squeeze(swdmeanice[p4ice])),2)))

    print ()
    print ('net LW (p4):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netLWice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netLW[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p4mod])) - np.nanmean(np.squeeze(netLW[p4ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p4mod])) - np.nanmean(np.squeeze(netLW[p4ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p4mod])) - np.nanmean(np.squeeze(netLW[p4ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p4mod])) - np.nanmean(np.squeeze(netLW[p4ship])),2)) + ')')

    print ()
    print ('Downwelling LW (p4):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(lwdmeanice[p4ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(lwdmean[p4ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p4mod])) - np.nanmean(np.squeeze(lwdmean[p4ship])))) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p4mod])) - np.nanmean(np.squeeze(lwdmean[p4ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p4mod])) - np.nanmean(np.squeeze(lwdmean[p4ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p4mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p4mod])) - np.nanmean(np.squeeze(lwdmean[p4ship])),2)) + ')')


    print ()
    print ('*******************')
    print ()

    print ()
    print ('net SW (p5):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netSW[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSW[p5ship])),2)) + '')
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSWice[p5ice])),2)) + '')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSW[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSWice[p5ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSW[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSWice[p5ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSW[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p5mod])) - np.nanmean(np.squeeze(netSWice[p5ice])),2)))

    print ()
    print ('Downwelling SW (p5):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(swdmean[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmean[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmean[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmean[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p5mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmean[p5ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p5mod])) - np.nanmean(np.squeeze(swdmeanice[p5ice])),2)))

    print ()
    print ('net LW (p5):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netLWice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netLW[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p5mod])) - np.nanmean(np.squeeze(netLW[p5ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p5mod])) - np.nanmean(np.squeeze(netLW[p5ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p5mod])) - np.nanmean(np.squeeze(netLW[p5ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p5mod])) - np.nanmean(np.squeeze(netLW[p5ship])),2)) + ')')

    print ()
    print ('Downwelling LW (p5):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(lwdmeanice[p5ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(lwdmean[p5ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p5mod])) - np.nanmean(np.squeeze(lwdmean[p5ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p5mod])) - np.nanmean(np.squeeze(lwdmean[p5ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p5mod])) - np.nanmean(np.squeeze(lwdmean[p5ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p5mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p5mod])) - np.nanmean(np.squeeze(lwdmean[p5ship])),2)) + ')')


    print ()
    print ('*******************')
    print ()

    print ()
    print ('net SW (p6):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netSWice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netSW[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSW[p6ship])),2)) + '')
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSWice[p6ice])),2)) + '')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSW[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSWice[p6ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSW[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSWice[p6ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSW[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWnet'][p6mod])) - np.nanmean(np.squeeze(netSWice[p6ice])),2)))

    print ()
    print ('Downwelling SW (p6):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(swdmean[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmean[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmean[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmean[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p6mod])),2)))
    print ('    Ship =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmean[p6ship])),2)))
    print ('    Ice =  ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['SWd'][p6mod])) - np.nanmean(np.squeeze(swdmeanice[p6ice])),2)))

    print ()
    print ('net LW (p6):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(netLWice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(netLW[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWnet'][p6mod])) - np.nanmean(np.squeeze(netLW[p6ship])))) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWnet'][p6mod])) - np.nanmean(np.squeeze(netLW[p6ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWnet'][p6mod])) - np.nanmean(np.squeeze(netLW[p6ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWnet'][p6mod])) - np.nanmean(np.squeeze(netLW[p6ship])),2)) + ')')

    print ()
    print ('Downwelling LW (p6):')
    print ('Obs(ice) = ' + str(np.round(np.nanmean(np.squeeze(lwdmeanice[p6ice])),2)))
    print ('Obs(ship) = ' + str(np.round(np.nanmean(np.squeeze(lwdmean[p6ship])),2)))
    print ('ECMWF_IFS = ' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data3['fixed_radiation']['LWd'][p6mod])) - np.nanmean(np.squeeze(lwdmean[p6ship])),2)) + ')')
    print ('UM_CASIM-100 = ' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data2['fixed_radiation']['LWd'][p6mod])) - np.nanmean(np.squeeze(lwdmean[p6ship])),2)) + ')')
    print ('UM_RA2T = ' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data4['fixed_radiation']['LWd'][p6mod])) - np.nanmean(np.squeeze(lwdmean[p6ship])),2)) + ')')
    print ('UM_RA2M = ' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p6mod])),2)) + ' / (' + str(np.round(np.nanmean(np.squeeze(data1['fixed_radiation']['LWd'][p6mod])) - np.nanmean(np.squeeze(lwdmean[p6ship])),2)) + ')')

def plot_sfcAlbedo(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    ######################################################
    #### SIMPLE CALCULATION OF ALBEDO
    ######################################################
    p5 = np.where(np.logical_and(data3['time'] >= 240.0, data3['time'] < 247.0))
    p5ship = np.where(np.logical_and(obs['albedo']['DoY'] >= 240.0, obs['albedo']['DoY'] < 247.0))
    p5ice = np.where(np.logical_and(obs['fixed_radiation']['time_ice'] >= 240.0, obs['fixed_radiation']['time_ice'] < 247.0))

    albedo_ice = obs['fixed_radiation']['SWu_ice'] / obs['fixed_radiation']['SWd_ice']



    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(6,5))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    plt.plot(obs['fixed_radiation']['time_ice'][p5ice], albedo_ice[p5ice], color = 'grey', label = 'Ice_station')
    plt.plot(obs['albedo']['DoY'][p5ship], obs['albedo']['a'][p5ship], 'ko', label = 'Ship')
    plt.plot(data3['time'][p5],data3['sfc_albedo'][p5], color = 'gold', label = label3)
    # plt.plot(data1['time'],data1['surface_albedo'],color = 'darkblue')
    # plt.plot(data2['time'],data2['surface_albedo'],color = 'mediumseagreen')
    # plt.plot(data4['time'],data4['surface_albedo'],color = 'steelblue')
    plt.legend()
    # plt.ylim([0.4,1.0])
    # plt.grid('on')

    # ax1.set_xticklabels([0,' ',0.2,' ',0.4,' ',0.6,' ',0.8,' ',1.0])

    plt.show()

    # print (data1['surface_downwelling_SW_radiation'])
    # print (data1['surface_upwelling_SW_radiation'])
    # print (data1['surface_albedo'])


def plot_Precipitation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting combined timeseries and PDFs of precipitation terms:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(12,4.5))
    # plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
    #         hspace = 0.4, wspace = 0.1)

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation on different timestep
    time_radice = calcTime_Mat2DOY(datenums_radice)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. temp_1.5m -> sfc_temp_2m
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    data1['rainfall_flux'][data1['rainfall_flux'] < 0] = np.nan
    data2['rainfall_flux'][data2['rainfall_flux'] < 0] = np.nan
    data4['rainfall_flux'][data4['rainfall_flux'] < 0] = np.nan
    data1['snowfall_flux'][data1['snowfall_flux'] < 0] = np.nan
    data2['snowfall_flux'][data2['snowfall_flux'] < 0] = np.nan
    data4['snowfall_flux'][data4['snowfall_flux'] < 0] = np.nan

    # flx_ls_rain = np.nansum(data3['flx_ls_rain'],1)
    # flx_conv_rain = np.nansum(data3['flx_conv_rain'],1)
    # flx_conv_rain[flx_conv_rain < 0] = np.nan
    # flx_ls_snow = np.nansum(data3['flx_ls_snow'],1)

    flx_conv_rain = data3['flx_conv_rain'][:,0]         #### just take surface value
    flx_conv_snow = data3['flx_conv_snow'][:,0]
    flx_ls_rain = data3['flx_ls_rain'][:,0]         #### just take surface value
    flx_ls_snow = data3['flx_ls_snow'][:,0]             #### assumes all precip which forms at altitudes
                                                        #### above evaporates/sublimes before it reaches
                                                        #### the surface

    # flx_ls_rain = np.nanmean(data3['flx_ls_rain'],1)         #### take average rate
    # flx_ls_snow = np.nanmean(data3['flx_ls_snow'],1)

    # flx_ls_rain = np.nansum(data3['flx_ls_rain'],1)         #### take total
    # flx_ls_snow = np.nansum(data3['flx_ls_snow'],1)

    #### remove flagged values
    flx_conv_rain[flx_conv_rain < 0] = np.nan
    flx_conv_snow[flx_conv_snow < 0] = np.nan
    flx_ls_rain[flx_ls_rain < 0] = np.nan
    flx_ls_snow[flx_ls_snow < 0] = np.nan

    ### doy drift flag
    drift = np.squeeze(np.where(np.logical_and(obs['pws']['doy'] >= 226.0, obs['pws']['doy'] <= 258.0)))

    # print (drift.shape)

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    precip1 = data1['rainfall_flux'][data1['hrly_flag']].data*3600 + data1['snowfall_flux'][data1['hrly_flag']].data*3600
    precip2 = data2['rainfall_flux'][data2['hrly_flag']].data*3600 + data2['snowfall_flux'][data2['hrly_flag']].data*3600
    if ifs_flag: precip3 = (flx_ls_rain[data3['hrly_flag']] + flx_ls_snow[data3['hrly_flag']] + flx_conv_rain[data3['hrly_flag']] + flx_conv_snow[data3['hrly_flag']])*3600
    precip4 = data4['rainfall_flux'][data4['hrly_flag']].data*3600 + data4['snowfall_flux'][data4['hrly_flag']].data*3600

    precip3[precip3 < 0] == np.nan

    res = 3         ### hourly resolution to plot

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    # ax  = fig.add_axes([0.75,0.15,0.2,0.7])   # left, bottom, width, height
    fig = plt.figure(figsize=(12,3.5))
    plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    ax  = fig.add_axes([0.09,0.18,0.6,0.76])   # left, bottom, width, height
    plt.plot(data2['time'], zeros,'--', color='lightgrey')
    plt.plot(obs['pws']['doy'][drift[0]],obs['pws']['prec_int'][drift[0]], color = 'grey', label = 'Obs_PWS', zorder = 1)
    if ifs_flag == True:
        plt.plot(data3['time_hrly'][::res], precip3[::res],
            'v', color = 'gold', markeredgecolor = 'orange',zorder = 1)
    plt.plot(data2['time_hrly'][::res], precip2[::res],
        '<', color = 'mediumseagreen', markeredgecolor = 'darkgreen', zorder = 1)
    plt.plot(data4['time_hrly'][::res], precip4[::res],
        '>', color = 'steelblue', markeredgecolor = 'darkslategrey', zorder = 1)
    plt.plot(data1['time_hrly'][::res], precip1[::res],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', zorder = 1)
    plt.ylabel('Precipitation flux [mm hr$^{-1}$]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.ylim([0,2])
    plt.xlabel('Date')
    plt.yscale('log')
    plt.ylim([3e-2,2e0])
    ax = plt.gca()
    nans = ax.get_ylim()
    for file in missing_files:
        ax.fill_between(np.arange(file, file + 1, 1/24.0), nans[0], nans[-1],
            facecolor = 'white',
            # hatch = 'x',
            zorder = 2)
    plt.legend(bbox_to_anchor=(0.0, 0.81, 1., .102), loc=3, ncol=2)

    ax  = fig.add_axes([0.76,0.27,0.22,0.6])   # left, bottom, width, height
    sns.distplot(precip1, hist=False, color="darkblue", kde_kws={"shade": True})
    sns.distplot(precip4, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(precip2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    sns.distplot(precip3, hist=False, color="gold", kde_kws={"shade": True})
    sns.distplot(obs['pws']['prec_int'][drift[0]], hist=False, color="grey")
    plt.xlim([0,0.4])
    plt.yscale('log')
    plt.ylim([1e-1,3e1])
    plt.xlabel('Precipitation flux [mm hr$^{-1}$]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TotalPrecip_oden-pws_metum_ifs-z0_casim-aeroprof_ra2t_fixedRA2T_newColours_Date.svg'
    # plt.savefig(fileout)
    plt.show()


    ifs_convprecip = data3['flx_conv_rain'][data3['hrly_flag'],0]*3600 + data3['flx_conv_snow'][data3['hrly_flag'],0]*3600
    ifs_convprecip[ifs_convprecip < 0] == np.nan

    plt.figure()
    ax = plt.gca()
    plt.plot(data3['time_hrly'], precip3)
    plt.plot(data3['time_hrly'], np.squeeze(ifs_convprecip))

    ax.set_xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.show()

def plotWinds(data1, data2, data3, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT MODEL WIND FIELDS
    ###################################

    print ('******')
    print ('')
    print ('Plotting wind timeseries:')
    print ('')

    ##################################################
    ##################################################
    #### 	FIGURE
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(10,6))
    plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    ### remove flagged points
    data1['uwind'][data1['uwind'] == -9999] = np.nan
    data1['vwind'][data1['vwind'] == -9999] = np.nan
    data1['wwind'][data1['wwind'] == -9999] = np.nan
    data2['uwind'][data2['uwind'] == -9999] = np.nan
    data2['vwind'][data2['vwind'] == -9999] = np.nan
    data2['wwind'][data2['wwind'] == -9999] = np.nan
    data3['uwind'][data3['uwind'] == -9999] = np.nan
    data3['uwind'][data3['uwind'] < -400] = np.nan
    data3['vwind'][data3['vwind'] == -9999] = np.nan
    data3['vwind'][data3['vwind'] < -400] = np.nan
    data3['wwind'][data3['wwind'] == -9999] = np.nan
    data3['wwind'][data3['wwind'] < -400] = np.nan

    plt.subplot(331)
    plt.pcolor(data1['time'],data1['height'],np.transpose(data1['uwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('um_ra2m, u [m/s]')
    plt.ylabel('Z [m]')
    plt.colorbar()
    plt.subplot(332)
    plt.pcolor(data1['time'],data1['height'],np.transpose(data1['vwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('um_ra2m, v [m/s]')
    plt.colorbar()
    plt.subplot(333)
    plt.pcolor(data1['time'],data1['height'],np.transpose(data1['wwind']), vmin = -0.2, vmax = 0.2)
    plt.ylim([0,5e3])
    plt.title('um_ra2m, w [m/s]')
    plt.colorbar()

    plt.subplot(334)
    plt.pcolor(data2['time'],data2['height'],np.transpose(data2['uwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('um_casim-100, u [m/s]')
    plt.ylabel('Z [m]')
    plt.colorbar()
    plt.subplot(335)
    plt.pcolor(data2['time'],data2['height'],np.transpose(data2['vwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('um_casim-100, v [m/s]')
    plt.colorbar()
    plt.subplot(336)
    plt.pcolor(data2['time'],data2['height'],np.transpose(data2['wwind']), vmin = -0.2, vmax = 0.2)
    plt.ylim([0,5e3])
    plt.title('um_casim-100, w [m/s]')
    plt.colorbar()

    plt.subplot(337)
    plt.pcolor(data3['time'],np.nanmean(data3['height'],0),np.transpose(data3['uwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('ecmwf_ifs, u [m/s]')
    plt.ylabel('Z [m]')
    plt.colorbar()
    plt.subplot(338)
    plt.pcolor(data3['time'],np.nanmean(data3['height'],0),np.transpose(data3['vwind']), vmin = -30, vmax = 30)
    plt.ylim([0,5e3])
    plt.title('ecmwf_ifs, v [m/s]')
    plt.colorbar()
    plt.subplot(339)
    plt.pcolor(data3['time'],np.nanmean(data3['height'],0),np.transpose(data3['wwind']), vmin = -0.2, vmax = 0.2)
    plt.ylim([0,5e3])
    plt.title('ecmwf_ifs, w [m/s]')
    plt.colorbar()
    plt.xlabel('DOY')

    plt.show()

def plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Plotting BL depth timeseries:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ## save data into temp variables to allow subsampling
    #################################################################
    bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    if ifs_flag == True:
        bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    else:
        bldepth3 = data3['bl_depth'][data3['hrly_flag']]

    #################################################################
    ## convert model inversion timesteps
    #################################################################
    data1['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data1['inversions']['mday']))
    data2['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data2['inversions']['mday']))
    data3['inversions']['doy'] = calcTime_Mat2DOY(np.squeeze(data3['inversions']['mday']))

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['inversions']['doy'] >= 225.9, obs['inversions']['doy'] <= 258.0))

    ### save in dict for ease
    obs['inversions']['doy_drift'] = obs['inversions']['doy'][drift]

    #### split into melt and freeze in scatter plots:
    ####        all model data share a timestamp
    melt = np.where(np.logical_and(data1['time_hrly'] >= obs['inversions']['doy_drift'][0], data1['time_hrly'] < 240.0))
    freeze = np.where(data1['time_hrly'] >= 240.0)
    #### allow some leeway in radiosonde timesteps
    obsmelt = np.where(np.logical_and(obs['inversions']['doy'] >= data1['time_hrly'][0]-0.2, obs['inversions']['doy'] < 240.0))
    obsfreeze = np.where(np.logical_and(obs['inversions']['doy'] >= 240.0, obs['inversions']['doy'] <= data1['time_hrly'][-1]))

    # print ('obsmelt.shape = ', obsmelt.shape)

    #### make inversion tempvars to allow for easy subsampling
    inv1 = np.squeeze(data1['inversions']['invbase'][data1['hrly_flag'],0])
    inv2 = np.squeeze(data2['inversions']['invbase'][data2['hrly_flag'],0])
    inv3 = np.squeeze(data3['inversions']['invbase'][data3['hrly_flag'],0])

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(13,7))

    #################################################################
    ## create figure and axes instances
    #################################################################
    ax  = fig.add_axes([0.08,0.56,0.45,0.36])   # left, bottom, width, height
    plt.plot(np.squeeze(obs['inversions']['doy_drift']),np.squeeze(obs['inversions']['sfmlheight'][drift]),
        color = 'k', label = 'Obs_Radiosondes')
    plt.plot(data1['time_hrly'][::6], bldepth1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(data2['time_hrly'][::6], bldepth2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(data3['time_hrly'][::6], bldepth3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    plt.legend(bbox_to_anchor=(0.0, 0.67, 1., .102), loc=4, ncol=2)
    plt.title('BL_depth / sfmlheight [m]')
    ax.set_xlim([doy[0],doy[-1]])
    # plt.xlabel('Day of year')
    plt.ylabel('Z [m]')

    ax  = fig.add_axes([0.08,0.1,0.45,0.36])   # left, bottom, width, height
    plt.plot(np.squeeze(obs['inversions']['doy_drift']), np.squeeze(obs['inversions']['invbase'][drift]),
        color = 'k', label = 'Obs: main inversion')
    plt.plot(data1['time_hrly'][::6], inv1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(data2['time_hrly'][::6], inv2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(data1['time_hrly'][::6], inv3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    # plt.legend()
    plt.title('Main inversion height [m]')
    ax.set_xlim([doy[0],doy[-1]])
    plt.xlabel('Day of year')
    plt.ylabel('Z [m]')

    #### -----------------------------------------------------------------
    #### scatterplots
    #### -----------------------------------------------------------------
    ax  = fig.add_axes([0.6,0.64,0.15,0.22])   # left, bottom, width, height
    ymax1 = 750
    xmax1 = ymax1
    blmelt1 = bldepth1[melt]
    blmelt2 = bldepth2[melt]
    blmelt3 = bldepth3[melt]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsmelt]), blmelt3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{SML}$ [m]')
    plt.ylabel('Model$_{SML}$ [m]')
    plt.title('Melt')

    ax  = fig.add_axes([0.83,0.64,0.15,0.22])   # left, bottom, width, height
    ymax1 = 1500
    xmax1 = ymax1
    blfreeze1 = bldepth1[freeze]
    blfreeze2 = bldepth2[freeze]
    blfreeze3 = bldepth3[freeze]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['sfmlheight'][obsfreeze]), blfreeze3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{SML}$ [m]')
    plt.ylabel('Model$_{SML}$ [m]')
    plt.title('Freeze')

    ax  = fig.add_axes([0.6,0.18,0.15,0.22])   # left, bottom, width, height
    ymax1 = 2700
    xmax1 = ymax1
    invmelt1 = inv1[melt]
    invmelt2 = inv2[melt]
    invmelt3 = inv3[melt]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsmelt]), invmelt3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{inv}$ [m]')
    plt.ylabel('Model$_{inv}$ [m]')
    plt.title('Melt')

    ax  = fig.add_axes([0.83,0.18,0.15,0.22])   # left, bottom, width, height
    ymax1 = 2700
    xmax1 = ymax1
    invfreeze1 = inv1[freeze]
    invfreeze2 = inv2[freeze]
    invfreeze3 = inv3[freeze]
    plt.plot([0,xmax1],[0, ymax1], '--', color = 'lightgrey')
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze1[::6],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = label1)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze2[::6],
        'v', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = label2)
    plt.plot(np.squeeze(obs['inversions']['invbase'][obsfreeze]), invfreeze3[::6],
        'd', color = 'gold', markeredgecolor = 'saddlebrown',  label = label3)
    plt.xlim([0,xmax1])
    plt.ylim([0,ymax1])
    plt.xlabel('Obs$_{inv}$ [m]')
    plt.ylabel('Model$_{inv}$ [m]')
    plt.title('Freeze')


    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_oden_metum_ifs_casim-100.svg'
    # fileout = '../FIGS/comparisons/BLDepth_calcInvHeights_timeseries_wScatter-SplitSeason_metum.svg'
    plt.savefig(fileout)
    plt.show()

def plot_BLType(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    ###################################
    ## Plot occurrences of BL type from UM
    ###################################

    print ('******')
    print ('')
    print ('Plotting BL type histograms:')
    print ('')

    ##################################################
    ##################################################
    #### 	set up figure properties
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 18

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(9,8))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.5, left = 0.08,
            hspace = 0.4, wspace = 0.10)

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ## split data into relative contributions from each BL type
    #################################################################

    data1['histogram_n'], data1['histogram_bins'] = np.histogram(data1['bl_type'], bins = np.arange(1,9))
    data2['histogram_n'], data2['histogram_bins'] = np.histogram(data2['bl_type'], bins = np.arange(1,9))

    ### calculate total number of occurrences for normalised charts
    total1 = float(np.sum(data1['histogram_n']))
    total2 = float(np.sum(data2['histogram_n']))

    ### create arranys of instances for e.g. type I, type II etc. between each simulation
    types = {}
    for i in range(0,7): types[i+1] = [data1['histogram_n'][i]/total1, data2['histogram_n'][i]/total2]

    #### list of BL types from UM documentation
    doc = ['1: S BL', '2: Sc over stable SL', '3: Well-mixed BL', '4: Unstable BL, dSc not o/Cu',
        '5: dSc o/Cu', '6: Cu-capped BL', '7: Shear-dom unstable BL']

    # Type I: Stable boundary layer (with or without cloud)
    # Type II: Boundary layer with stratocumulus over a stable near-surface layer
    # Type III: Well mixed boundary layer
    # Type IV: Unstable boundary layer with a DSC layer not over cumulus
    # Type V: Boundary layer with a DSC layer over cumulus
    # Type VI: Cumulus-capped boundary layer
    # Type VII: Shear-dominated unstable layer

    #################################################################
    ## create figure and axes instances
    #################################################################
    ax = plt.gca()

    plt.bar([0,1], types[1], label = doc[0])
    plt.bar([0,1], types[2], bottom = types[1], label = doc[1]); bars = np.add(types[1], types[2]).tolist()
    for i in range(3,8):
        # print(i)
        if i == 4:
            plt.bar([0,1], types[i], bottom = bars, hatch = '/', label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
        else:
            plt.bar([0,1], types[i], bottom = bars, label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
    plt.xticks([0,1], [label1, label2])
    plt.legend(bbox_to_anchor=(1.1, 0.3, 1., .102), loc=4, ncol=1)
    plt.title('BL type occurrences (normalised)')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLType_oden_metum_casim-100.svg'
    plt.savefig(fileout)
    plt.show()


    #### ------------------------------------------------------------------
    #### ------------------------------------------------------------------
    ####        SEASONAL SPLIT
    #### ------------------------------------------------------------------
    #### ------------------------------------------------------------------


    melt = np.where(data1['time'] < 240.0)
    freeze = np.where(data1['time'] >= 240.0)

    #################################################################
    ## split data into relative contributions from each BL type
    #################################################################

    data1['histogram_nMelt'], data1['histogram_bins'] = np.histogram(data1['bl_type'][melt], bins = np.arange(1,9))
    data2['histogram_nMelt'], data2['histogram_bins'] = np.histogram(data2['bl_type'][melt], bins = np.arange(1,9))

    data1['histogram_nFreeze'], data1['histogram_bins'] = np.histogram(data1['bl_type'][freeze], bins = np.arange(1,9))
    data2['histogram_nFreeze'], data2['histogram_bins'] = np.histogram(data2['bl_type'][freeze], bins = np.arange(1,9))

    ### calculate total number of occurrences for normalised charts
    total1melt = float(np.sum(data1['histogram_nMelt']))
    total2melt = float(np.sum(data2['histogram_nMelt']))
    total1freeze = float(np.sum(data1['histogram_nFreeze']))
    total2freeze = float(np.sum(data2['histogram_nFreeze']))


    ### create arranys of instances for e.g. type I, type II etc. between each simulation
    types = {}
    for i in range(0,7): types[i+1] = [data1['histogram_nMelt'][i]/total1melt, data1['histogram_nFreeze'][i]/total1freeze,
                                        data2['histogram_nMelt'][i]/total2melt, data2['histogram_nFreeze'][i]/total2freeze]

    #### list of BL types from UM documentation
    doc = ['1: Stable BL', '2: Sc over stable SL', '3: Well-mixed BL', '4: Unstable BL, dSc not o/Cu',
        '5: dSc o/Cu', '6: Cu-capped BL', '7: Shear-dom unstable BL']

    # Type I: Stable boundary layer (with or without cloud)
    # Type II: Boundary layer with stratocumulus over a stable near-surface layer
    # Type III: Well mixed boundary layer
    # Type IV: Unstable boundary layer with a DSC layer not over cumulus
    # Type V: Boundary layer with a DSC layer over cumulus
    # Type VI: Cumulus-capped boundary layer
    # Type VII: Shear-dominated unstable layer

    #################################################################
    ## create figure and axes instances
    #################################################################
    plt.figure(figsize=(12,7))
    # plt.rc('figure',titlesize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.68, left = 0.05,
            hspace = 0.4, wspace = 0.10)
    ax = plt.gca()

    plt.bar(np.arange(0,4), types[1], label = doc[0])
    plt.bar(np.arange(0,4), types[2], bottom = types[1], label = doc[1]); bars = np.add(types[1], types[2]).tolist()
    for i in range(3,8):
        # print(i)
        if i == 4:
            plt.bar(np.arange(0,4), types[i], bottom = bars, hatch = '/', label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
        else:
            plt.bar(np.arange(0,4), types[i], bottom = bars, label = doc[i-1]); bars = np.add(bars, types[i]).tolist()
    plt.xticks(np.arange(0,4), [label1 + '\nMelt', label1 + '\nFreeze', label2 + '\nMelt', label2 + '\nFreeze'])
    plt.legend(bbox_to_anchor=(0.5, 0.3, 1., .102), loc=4, ncol=1)
    plt.title('BL type occurrences (normalised)')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/BLType_oden_metum_casim-100_splitSeason.svg'
    plt.savefig(fileout)
    plt.show()

def plot_paperGLMAnalysis(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data4['temperature'][data4['temperature'] == -9999] = np.nan
    data5['temperature'][data5['temperature'] <= 0] = np.nan

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan
    data4['q'][data4['q'] == -9999] = np.nan
    data5['q'][data5['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'temp')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    Tmin = -45
    Tmax = 5
    ymax = 9000
    qmax = 4.0

    data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)

    data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])

    ##################################################
    ##################################################
    #### MEDIAN PROFILES
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7.5,9))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.93, left = 0.1,
            hspace = 0.2, wspace = 0.2)

    ####        need different indices for ifs and um_glm since latter not cloudnet-ed
    m_um = np.where(data5['time_hrly'][::6] < 240.0)
    f_um = np.where(data5['time_hrly'][::6] >= 240.0)
    m_ifs = np.where(data3['time_hrly'][::6] < 240.0)
    f_ifs = np.where(data3['time_hrly'][::6] >= 240.0)

    ###-------------------------
    plt.subplot(221)
    ax1 = plt.gca()
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data5['temp_anomalies'][:,m_um]),1) - np.nanstd(np.squeeze(data5['temp_anomalies'][:,m_um]),1),
    #     np.nanmedian(np.squeeze(data5['temp_anomalies'][:,m_um]),1) + np.nanstd(np.squeeze(data5['temp_anomalies'][:,m_um]),1),
    #     color = 'lightgrey', alpha = 0.3)
    # plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,m_um]),1) - np.nanstd(np.squeeze(data5['temp_anomalies'][:,m_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,m_um]),1) + np.nanstd(np.squeeze(data5['temp_anomalies'][:,m_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,m_um]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend(bbox_to_anchor=(0.8, 1.03, 1., .102), loc=4, ncol=3)
    plt.legend(bbox_to_anchor=(1.32, 1.01, 1., .102), loc=4, ncol=3)
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-2.0,1.0])
    # plt.xlabel('T bias [K]')

    ###-------------------------
    plt.subplot(222)
    ax1 = plt.gca()
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data5['q_anomalies'][:,m_um]),1) - np.nanstd(np.squeeze(data5['q_anomalies'][:,m_um]),1),
    #     np.nanmedian(np.squeeze(data5['q_anomalies'][:,m_um]),1) + np.nanstd(np.squeeze(data5['q_anomalies'][:,m_um]),1),
    #     color = 'lightgrey', alpha = 0.3)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,m_um]),1) - np.nanstd(np.squeeze(data5['q_anomalies'][:,m_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,m_um]),1) + np.nanstd(np.squeeze(data5['q_anomalies'][:,m_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,m_ifs]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,m_ifs]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,m_ifs]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,m_ifs]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,m_ifs]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,m_ifs]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,m_ifs]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,m_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,m_ifs]),1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,m_um]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend()
    # plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.35,0.35])
    plt.grid('on')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Melt', rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ###-------------------------
    plt.subplot(223)
    ax1 = plt.gca()
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data5['temp_anomalies'][:,f_um]),1) - np.nanstd(np.squeeze(data5['temp_anomalies'][:,f_um]),1),
    #     np.nanmedian(np.squeeze(data5['temp_anomalies'][:,f_um]),1) + np.nanstd(np.squeeze(data5['temp_anomalies'][:,f_um]),1),
    #     color = 'lightgrey', alpha = 0.3)
    # plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,f_um]),1) - np.nanstd(np.squeeze(data5['temp_anomalies'][:,f_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,f_um]),1) + np.nanstd(np.squeeze(data5['temp_anomalies'][:,f_um]),1), data1['universal_height'],
    #     '--', color = 'grey', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'steelblue', label = label4, zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,f_um]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-2.0,1.0])
    plt.xlabel('T bias [K]')

    ###-------------------------
    plt.subplot(224)
    ax1 = plt.gca()
    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data5['q_anomalies'][:,f_um]),1) - np.nanstd(np.squeeze(data5['q_anomalies'][:,f_um]),1),
    # np.nanmedian(np.squeeze(data5['q_anomalies'][:,f_um]),1) + np.nanstd(np.squeeze(data5['q_anomalies'][:,f_um]),1),
    # color = 'lightgrey', alpha = 0.3)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,f_um]),1) - np.nanstd(np.squeeze(data5['q_anomalies'][:,f_um]),1), data1['universal_height'],
    # '--', color = 'grey', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,f_um]),1) + np.nanstd(np.squeeze(data5['q_anomalies'][:,f_um]),1), data1['universal_height'],
    # '--', color = 'grey', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,f_ifs]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,f_ifs]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,f_ifs]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,f_ifs]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,f_ifs]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,f_ifs]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,f_ifs]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,f_ifs]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'steelblue', label = label4, zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,f_ifs]),1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,f_um]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend()
    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.35,0.35])
    plt.grid('on')
    ax2 = ax1.twinx()
    ax2.set_ylabel('Freeze', rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ###-------------------------
    fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_ifs_UMGlobal_ra2m_ra2t_casim-aeroprof_fixedRA2T.svg'
    plt.savefig(fileout, dpi = 300)
    plt.show()

    Zindex1 = np.where(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1) < 0)
    Zindex3 = np.where(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1) < 0)
    print ('Z = ')
    print (data1['universal_height'][Zindex3])
    print (data1['universal_height'])

    print ('UM_RA2M = ')
    print (np.round(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,f_ifs]),1),2))

    print ('UM_CASIM-100 = ')
    print (np.round(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,f_ifs]),1),2))

    print ('ECMWF_IFS = ')
    print (np.round(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,f_ifs]),1),2))

    print ('UM_RA2T = ')
    print (np.round(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,f_ifs]),1),2))

    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,m_ifs]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,m_ifs]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,m_ifs]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,m_ifs]),1),2)))


def plot_paperRadiosondes(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data4['temperature'][data4['temperature'] == -9999] = np.nan
    data5['temperature'][data5['temperature'] <= 0] = np.nan

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan
    data4['q'][data4['q'] == -9999] = np.nan
    data5['q'][data5['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'temp')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    Tmin = -45
    Tmax = 5
    ymax = 9000
    qmax = 4.0

    data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)

    data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])

    ##################################################
    ##################################################
    #### MEDIAN PROFILES
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7.5,5.5))
    plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
            hspace = 0.3, wspace = 0.2)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)


    ###-------------------------
    plt.subplot(121)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
        np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
        np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
         color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
        np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
        np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend(bbox_to_anchor=(0.9, 1.03, 1., .102), loc=4, ncol=2)
    plt.legend(bbox_to_anchor=(1.25, 1.03, 1., .102), loc=4, ncol=3)
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-2.0,1.0])
    plt.xlabel('T bias [K]')

    ###-------------------------
    plt.subplot(122)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
        np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
        np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
        np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
        np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend()
    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    plt.grid('on')

    ###-------------------------
    fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_ifs_casim-100_ra2t_ra2m_wUMGlobal.svg'
    # plt.savefig(fileout)
    plt.show()

    ##################################################
    ##################################################
    #### PCOLOR TIMESERIES
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(12,11))

    Tmin = -45
    Tmax = 5
    qbiasmin = -1.0
    qbiasmax = 1.0
    ymax = 9000
    yticklist = [0, 3e3, 6e3, 9e3]

    #########-------------------------------------------------------------------------------------------
    ####       DEFINE PERIODS
    ####               all model data share a timestamp
    p3 = np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] <= 230.0))
    p4 = np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] <= 240.0))
    p5 = np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] <= 247.0))
    p6 = np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] <= 251.0))
    p7 = np.where(np.logical_and(data1['time_hrly'][::6] >= 251.0, data1['time_hrly'][::6] <= 255.5))
    p8 = np.where(data1['time_hrly'][::6] >= 255.5)
    #   p3: np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] < 230.0))
    #   p4: np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] < 240.0))
    #   p5: np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] < 247.0))
    #   p6: np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] < 251.0))

    #########-------------------------------------------------------------------------------------------

    ### -------------------------------
    ### model temp anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.05,0.81,0.4,0.12])   # left, bottom, width, height
    sfig1 = plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
        vmin = Tmin, vmax = Tmax)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][0]], data1['time_6hrly'][p3[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p3[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p3[0][-1]], 8000, 'w.')
    plt.annotate('P3', xy=(227,6800), xytext=(227.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][0]], data1['time_6hrly'][p4[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p4[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p4[0][-1]], 8000, 'w.')
    plt.annotate('P4', xy=(234,6800), xytext=(234.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][0]], data1['time_6hrly'][p5[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p5[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p5[0][-1]], 8000, 'w.')
    plt.annotate('P5', xy=(242.5,6800), xytext=(242.6,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][0]], data1['time_6hrly'][p6[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p6[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p6[0][-1]], 8000, 'w.')
    plt.annotate('P6', xy=(248,6800), xytext=(248.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][0]], data1['time_6hrly'][p7[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p7[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p7[0][-1]], 8000, 'w.')
    plt.annotate('P7', xy=(252.25,6800), xytext=(252.5,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p8
    plt.plot([data1['time_6hrly'][p8[0][0]], data1['time_6hrly'][p8[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p8[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p8[0][-1]], 8000, 'w.')
    plt.annotate('P8', xy=(256,6800), xytext=(256.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    plt.ylabel('Z [km]')
    plt.xlabel('Date')
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.title('T [$^{\circ}$C]')
    ax2 = ax.twinx()
    ax2.set_ylabel('Radiosondes', rotation = 270, labelpad = 15)
    ax2.set_yticks([])
    cbaxes1 = fig.add_axes([0.1, 0.98, 0.3, 0.015])
    cb1 = plt.colorbar(sfig1, cax = cbaxes1, orientation = 'horizontal')

    ax  = fig.add_axes([0.05,0.56,0.4,0.12])   # left, bottom, width, height
    # dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # data3['temp_anomalies'] = dat3
    sfig2 = plt.pcolor(data3['time_6hrly'], data1['universal_height'], data3['temp_anomalies'],
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.colorbar()
    # plt.set_cmap('seismic')
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label3, rotation = 270, labelpad = 15)
    ax2.set_yticks([])
    plt.title('T [K]')
    # plt.title(label3 + ' - Radiosondes, T[degC]')
    cbaxes2 = fig.add_axes([0.1, 0.73, 0.3, 0.015])
    cb2 = plt.colorbar(sfig2, cax = cbaxes2, orientation = 'horizontal')

    ax  = fig.add_axes([0.05,0.39,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], data2['temp_anomalies'],
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.colorbar()
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label2[:8], rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ax  = fig.add_axes([0.05,0.22,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data4['time_6hrly'],data1['universal_height'], data4['temp_anomalies'],
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.colorbar()
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label4[:-4], rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ax  = fig.add_axes([0.05,0.05,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], data1['temp_anomalies'],
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.xlabel('Date')
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label1, rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ###-------------------------------------
    ### -------------------------------
    ### model q anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.55,0.81,0.4,0.12])   # left, bottom, width, height
    sfig3 = plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_driftSondes_UM']),
        vmin = 0, vmax = qmax)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][0]], data1['time_6hrly'][p3[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p3[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p3[0][-1]], 8000, 'w.')
    plt.annotate('P3', xy=(227,6800), xytext=(227.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][0]], data1['time_6hrly'][p4[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p4[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p4[0][-1]], 8000, 'w.')
    plt.annotate('P4', xy=(234,6800), xytext=(234.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][0]], data1['time_6hrly'][p5[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p5[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p5[0][-1]], 8000, 'w.')
    plt.annotate('P5', xy=(242.5,6800), xytext=(242.6,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][0]], data1['time_6hrly'][p6[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p6[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p6[0][-1]], 8000, 'w.')
    plt.annotate('P6', xy=(248,6800), xytext=(248.1,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][0]], data1['time_6hrly'][p7[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p7[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p7[0][-1]], 8000, 'w.')
    plt.annotate('P7', xy=(252.25,6800), xytext=(252.5,6801), fontsize = 12, color = 'w')
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ## p8
    plt.plot([data1['time_6hrly'][p8[0][0]], data1['time_6hrly'][p8[0][-1]]], [8000, 8000], 'w--')
    plt.plot(data1['time_6hrly'][p8[0][0]], 8000, 'w.')
    plt.plot(data1['time_6hrly'][p8[0][-1]], 8000, 'w.')
    plt.annotate('P8', xy=(256,6800), xytext=(256.1,6801), fontsize = 12, color = 'w')
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'lightgrey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.colorbar()
    plt.ylabel('Z [km]')
    plt.xlabel('Date')
    plt.title('q [g kg$^{-1}$]')
    ax2 = ax.twinx()
    ax2.set_ylabel('Radiosondes', rotation = 270, labelpad = 15)
    ax2.set_yticks([])
    cbaxes3 = fig.add_axes([0.6, 0.98, 0.3, 0.015])
    cb3 = plt.colorbar(sfig3, cax = cbaxes3, orientation = 'horizontal')

    ax  = fig.add_axes([0.55,0.56,0.4,0.12])   # left, bottom, width, height
    sfig4 = plt.pcolor(data3['time_6hrly'], data1['universal_height'], data3['q_anomalies'],
        vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label3, rotation = 270, labelpad = 15)
    ax2.set_yticks([])
    plt.title('q [g kg$^{-1}$]')
    cbaxes4 = fig.add_axes([0.6, 0.73, 0.3, 0.015])
    cb4 = plt.colorbar(sfig4, cax = cbaxes4, orientation = 'horizontal')

    ax  = fig.add_axes([0.55,0.39,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'], data1['universal_height'], data2['q_anomalies'],
        vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label2[:8], rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ax  = fig.add_axes([0.55,0.22,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data4['time_6hrly'], data1['universal_height'], data4['q_anomalies'],
        vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.ylabel('Z [km]')
    ax2 = ax.twinx()
    ax2.set_ylabel(label4[:-4], rotation = 270, labelpad = 15)
    ax2.set_yticks([])

    ax  = fig.add_axes([0.55,0.05,0.4,0.12])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'], data1['universal_height'], data1['q_anomalies'],
        vmin = qbiasmin, vmax = qbiasmax, cmap=mpl_cm.BrBG)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    #### plot periods
    ## p3
    plt.plot([data1['time_6hrly'][p3[0][-1]], data1['time_6hrly'][p3[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p4
    plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p5
    plt.plot([data1['time_6hrly'][p5[0][-1]], data1['time_6hrly'][p5[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p6
    plt.plot([data1['time_6hrly'][p6[0][-1]], data1['time_6hrly'][p6[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p7
    plt.plot([data1['time_6hrly'][p7[0][-1]], data1['time_6hrly'][p7[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ## p8
    # plt.plot([data1['time_6hrly'][p8[0][-1]], data1['time_6hrly'][p8[0][-1]]], [0, 9000], '--', color = 'grey', linewidth = 1)
    ##
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    ax.set_yticklabels([0,3,6,9])
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    # plt.colorbar()
    plt.ylabel('Z [km]')
    plt.xlabel('Date')
    ax2 = ax.twinx()
    ax2.set_ylabel(label1, rotation = 270, labelpad = 15)
    ax2.set_yticks([])
    # plt.title(label1 + ' - Radiosondes, T[degC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TimeSeriesProfiles_TandSpHum-BrBG_ifs_casim-aeroprof_ra2t_ra2m_Dates_fixedRA2T.png'
    # plt.savefig(fileout, dpi=300)
    plt.show()
    # plt.close()

    doyindex = np.where(np.logical_and(data1['time_6hrly'] >= 247.0, data1['time_6hrly'] < 251.0))
    print ('Z = ')
    print (data1['universal_height'])

    # print ('UM_RA2M = ')
    # print (np.round(data1['q_anomalies'][:,doyindex[0]],2))
    #
    # print ('UM_CASIM-100 = ')
    # print (np.round(data2['q_anomalies'][:,doyindex[0]],2))
    #
    # print ('ECMWF_IFS = ')
    # print (np.round(data3['q_anomalies'][:,doyindex[0]],2))
    #
    # print ('UM_RA2T = ')
    # print (np.round(data4['q_anomalies'][:,doyindex[0]],2))
    #
    # print (np.nanmin(np.round(data1['q_anomalies'][:,doyindex[0]],2)))
    # print (np.nanmin(np.round(data2['q_anomalies'][:,doyindex[0]],2)))
    # print (np.nanmin(np.round(data3['q_anomalies'][:,doyindex[0]],2)))
    # print (np.nanmin(np.round(data4['q_anomalies'][:,doyindex[0]],2)))

    print ('UM_CASIM-100 = ')
    print (np.round(np.nanmean(data2['temp_6hrly'][doyindex[0],:],0),2) - 273.15)

    print ('ECMWF_IFS = ')
    print (np.round(np.nanmean(data3['temp_6hrly_UM'][doyindex[0],:],0),2) - 273.15)


def plot_paperERAIProfiles(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data4['temperature'][data4['temperature'] == -9999] = np.nan
    data5['temperature'][data5['temperature'] <= 0] = np.nan

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan
    data4['q'][data4['q'] == -9999] = np.nan
    data5['q'][data5['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'temp')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    Tmin = -45
    Tmax = 5
    ymax = 9000
    qmax = 4.0

    data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)

    data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])

    ##################################################
    ##################################################
    #### MEDIAN PROFILES
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7.5,5.5))
    plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
            hspace = 0.3, wspace = 0.2)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)


    ###-------------------------
    plt.subplot(121)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
        np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
    #     np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
    #      color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
        np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
        color = 'salmon', alpha = 0.15)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'firebrick', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'firebrick', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
        np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'firebrick', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend(bbox_to_anchor=(0.9, 1.03, 1., .102), loc=4, ncol=2)
    plt.legend(bbox_to_anchor=(1.0, 1.03, 1., .102), loc=4, ncol=2)
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-2.0,1.5])
    plt.xlabel('T bias [K]')

    ###-------------------------
    plt.subplot(122)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
        np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
    #     np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
        np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
        color = 'salmon', alpha = 0.15)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'firebrick', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'firebrick', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
        np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'gold', label = label3, zorder = 4)
    plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'firebrick', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', label = label1, zorder = 3)
    plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend()
    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.4,0.4])#plt.xlim([-0.05,0.45])
    plt.grid('on')

    ###-------------------------
    fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_ifs_ra2m_ra2m-erai_wUMGlobal.svg'
    # plt.savefig(fileout)
    plt.show()

    Zindex1 = np.where(np.nanmedian(np.squeeze(data1['temp_anomalies']),1) < 0)
    Zindex3 = np.where(np.nanmedian(np.squeeze(data3['temp_anomalies']),1) < 0)
    print ('Z(universal) = ')
    # print (data1['universal_height'][Zindex3])
    print (data1['universal_height'])

    print ('UM_RA2M = ')
    print (np.round(np.nanmedian(np.squeeze(data1['temp_anomalies']),1),2))

    print ('UM_CASIM-100 = ')
    print (np.round(np.nanmedian(np.squeeze(data2['temp_anomalies']),1),2))

    print ('ECMWF_IFS = ')
    print (np.round(np.nanmedian(np.squeeze(data3['temp_anomalies']),1),2))

    print ('UM_RA2T = ')
    print (np.round(np.nanmedian(np.squeeze(data4['temp_anomalies']),1),2))

    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data1['temp_anomalies']),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data2['temp_anomalies']),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data3['temp_anomalies']),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data4['temp_anomalies']),1),2)))


def plot_paperCASIMNiceProfiles(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data4['temperature'][data4['temperature'] == -9999] = np.nan
    data5['temperature'][data5['temperature'] <= 0] = np.nan

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan
    data4['q'][data4['q'] == -9999] = np.nan
    data5['q'][data5['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'temp')
    data1, data2, data3, data4, data5, obs, drift = reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    Tmin = -45
    Tmax = 5
    ymax = 9000
    qmax = 4.0

    data3['temp_anomalies'] = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data1['temp_anomalies'] = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data2['temp_anomalies'] = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data4['temp_anomalies'] = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)
    data5['temp_anomalies'] = np.transpose(data5['temp_6hrly_UM'][:]) - np.transpose(obs['sondes']['temp_subsetSondes_UM'] + 273.15)

    data3['q_anomalies'] = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data1['q_anomalies'] = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data2['q_anomalies'] = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data4['q_anomalies'] = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])
    data5['q_anomalies'] = np.transpose(data5['q_6hrly_UM'][:])*1e3 - np.transpose(obs['sondes']['q_subsetSondes_UM'])

    ##################################################
    ##################################################
    #### MEDIAN PROFILES
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7.5,5.5))
    plt.subplots_adjust(top = 0.8, bottom = 0.15, right = 0.97, left = 0.1,
            hspace = 0.3, wspace = 0.2)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)


    ###-------------------------
    plt.subplot(121)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
        np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
        color = 'lavender', alpha = 0.3)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'slateblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'slateblue', linewidth = 0.5)

    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
    #     np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
    #      color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
        np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
        np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'goldenrod', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'goldenrod', linewidth = 0.5)

    plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2 + '_Cooper', zorder = 1)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'goldenrod', label = label4, zorder = 4)
    # plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'slateblue', label = label1, zorder = 3)
    # plt.plot(np.nanmedian(data5['temp_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend(bbox_to_anchor=(0.9, 1.03, 1., .102), loc=4, ncol=2)
    plt.legend(bbox_to_anchor=(1.25, 1.03, 1., .102), loc=4, ncol=2)
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-2.0,1.5])
    plt.xlabel('T bias [K]')

    ###-------------------------
    plt.subplot(122)
    ax1 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
        np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
        color = 'lavender', alpha = 0.3)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'slateblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'slateblue', linewidth = 0.5)

    # ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
    #     np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
        np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
        np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'goldenrod', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'goldenrod', linewidth = 0.5)

    plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'goldenrod', label = label4, zorder = 4)
    plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'mediumseagreen', label = label2, zorder = 1)
    # plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4], zorder = 2)
    plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],'.-' ,color = 'slateblue', label = label1, zorder = 3)
    # plt.plot(np.nanmedian(data5['q_anomalies'],1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.legend()
    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax1.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax1.set_yticks(axminor, minor = True)
    ax1.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.35,0.35])#plt.xlim([-0.05,0.45])
    plt.grid('on')

    ###-------------------------
    fileout = '../FIGS/comparisons/MedianProfiles_TandSpHum_casim-cooper_meyers_fletcher.svg'
    # plt.savefig(fileout)
    plt.show()

    ###------------------------
    ###     Analysis for paper
    ###------------------------

    print ('Cooper = ')
    print (np.nanmedian(data2['q_anomalies'],1))
    print ('Cooper[16] = ')
    print (np.nanmedian(data2['q_anomalies'],1)[16])
    print ('Meyers = ')
    print (np.nanmedian(data4['q_anomalies'],1))
    print ('Meyers[16] = ')
    print (np.nanmedian(data4['q_anomalies'],1)[16])
    print ('Cooper[16] - Meyers[16] = ')
    print (np.nanmedian(data2['q_anomalies'],1)[16] - np.nanmedian(data4['q_anomalies'],1)[16])
    print ('Cooper[16] - Meyers[16]) at Z = ')
    print (data1['universal_height'][16])

def plot_RadiosondesTemperature(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model temperature profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data4['temperature'][data4['temperature'] == -9999] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, obs, drift = reGrid_Sondes(data1, data2, data3, data4, obs, doy, ifs_flag, 'temp')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')
    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    # fig = plt.figure(figsize=(19,10))
    #
    # Tmin = -45
    # Tmax = 5
    # ymax = 8000
    #
    # ### -------------------------------
    # ### original data
    # ### ------------------------------
    # ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['temperature'][:,drift[0]],
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.set_cmap('viridis')
    # plt.ylabel('Z [m]')
    # plt.title('Sondes, T[degC]')
    #
    # ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['temp_6hrly'])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label3 + ', T[degC]')
    #
    # ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['temp_6hrly'])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label1 + ', T[degC]')
    #
    # ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['temp_6hrly'])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.xlabel('Day of year')
    # plt.title(label2 + ', T[degC]')
    #
    # ### -------------------------------
    # ### sonde and ifs data interpolated to um grid (Z<10km)
    # ### ------------------------------
    # ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title('Sondes(REGRID), T[degC]')
    #
    # ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['temp_hrly_UM'][::6])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label3 + '(REGRID), T[degC]')
    #
    # ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label1 + ', T[degC]')
    #
    # ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']])-273.15,
    #     vmin = Tmin, vmax = Tmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.xlabel('Day of year')
    # # plt.ylabel('Z [m]')
    # plt.title(label2 + ', T[degC]')
    #
    # ### -------------------------------
    # ### model anomalies wrt radiosondes
    # ### ------------------------------
    # ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
    #     vmin = Tmin, vmax = Tmax)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title('Sondes(REGRID), T[degC]')
    #
    # ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    # dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.set_cmap('seismic')
    # # plt.ylabel('Z [m]')
    # plt.title(label3 + '(REGRID) - Sondes(REGRID), T[K]')
    #
    # ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    # dat1 = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label1 + ' - Sondes(REGRID), T[K]')
    #
    # ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    # dat2 = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    # plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
    #     vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.xlabel('Day of year')
    # # plt.ylabel('Z [m]')
    # plt.title(label2 + ' - Sondes(REGRID), T[K]')
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # fileout = '../FIGS/comparisons/TemperatureProfiles_REGRID_10km_sondes_metum_ifs_casim-100.png'
    # # plt.savefig(fileout, dpi = 300)
    # plt.show()
    # # plt.close()


    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7,10))

    Tmin = -45
    Tmax = 5
    ymax = 9000
    yticklist = [0, 3e3, 6e3, 9e3]

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.15,0.82,0.85,0.13])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['temp_driftSondes_UM']),
        vmin = Tmin, vmax = Tmax)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title('Radiosondes, T[degC]')

    ax  = fig.add_axes([0.15,0.63,0.85,0.13])   # left, bottom, width, height
    dat3 = np.transpose(data3['temp_hrly_UM'][::6]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data3['temp_anomalies'] = dat3
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    plt.ylabel('Z [m]')
    plt.title(label3 + ' - Radiosondes, T[degC]')

    ax  = fig.add_axes([0.15,0.44,0.85,0.13])   # left, bottom, width, height
    dat1 = np.transpose(data1['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data1['temp_anomalies'] = dat1
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ' - Radiosondes, T[degC]')

    ax  = fig.add_axes([0.15,0.25,0.85,0.13])   # left, bottom, width, height
    dat2 = np.transpose(data2['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data2['temp_anomalies'] = dat2
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ' - Radiosondes, T[degC]')

    ax  = fig.add_axes([0.15,0.06,0.85,0.13])   # left, bottom, width, height
    dat4 = np.transpose(data4['temp_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['temp_driftSondes_UM'] + 273.15)
    data4['temp_anomalies'] = dat4
    plt.pcolor(data4['time_6hrly'],data1['universal_height'], dat4,
        vmin = -5.0, vmax = 5.0, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    plt.ylabel('Z [m]')
    plt.title(label4[:-4] + ' - Radiosondes, T[degC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/TemperatureProfiles_REGRID_anomOnly_sondes_ifs_ra2m_casim-100_ra2t.png'
    # plt.savefig(fileout, dpi=300)
    plt.show()
    # plt.close()

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(11,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)

    plt.subplot(131)
    ax1 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'darkblue', label = label1)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1),
        np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) - np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1) + np.nanstd(data1['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(data2['temp_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'mediumseagreen', label = label2)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1),
        np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) - np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['temp_anomalies'],1) + np.nanstd(data2['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['temp_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'gold', label = label3)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1),
        np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) - np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['temp_anomalies'],1) + np.nanstd(data3['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data4['temp_anomalies'],1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4])
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1),
        np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) - np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['temp_anomalies'],1) + np.nanstd(data4['temp_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.legend()
    plt.ylim([0,1e4])
    plt.xlim([-2.0,1.0])#plt.xlim([-1.6,1.0])
    plt.ylabel('Z [m]')
    plt.xlabel('Temperature bias [K]')
    plt.grid('on')
    plt.title('Total drift')

    plt.subplot(132)
    ax2 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'darkblue', label = label1 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,melt]),1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'mediumseagreen', label = label2 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,melt]),1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'gold', label = label3 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,melt]),1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4] + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,melt]),1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-2.0,1.0])#plt.xlim([-1.6,1.0])
    plt.xlabel('Temperature bias [K]')
    plt.title('Melt')

    plt.subplot(133)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'darkblue', label = label1 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,freeze]),1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'mediumseagreen', label = label2 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,freeze]),1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'gold', label = label3 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,freeze]),1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4] + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,freeze]),1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    plt.grid('on')
    plt.ylim([0,1e4])
    plt.xlim([-2.0,1.0])#plt.xlim([-1.6,1.0])
    plt.xlabel('Temperature bias [K]')
    plt.title('Freeze')

    fileout = '../FIGS/comparisons/TemperatureMedianProfiles_metum_ifs_casim-100_erai-glm.svg'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesQ(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model thetaE profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan
    data4['q'][data4['q'] == -9999] = np.nan

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, data4, obs, drift = reGrid_Sondes(data1, data2, data3, data4, obs, doy, ifs_flag, 'q')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    # ### -------------------------------
    # ### Build figure (timeseries)
    # ### -------------------------------
    # fig = plt.figure(figsize=(19,10))
    #
    # qmax = 4.0
    # ymax = 8000
    #
    # ### -------------------------------
    # ### original data
    # ### ------------------------------
    # ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['mr'][:,drift[0]],
    #     vmin = 0.0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.set_cmap('viridis')
    # plt.ylabel('Z [m]')
    # plt.title('Sondes, q [g/kg]')
    #
    # ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['q_6hrly'])*1e3,
    #     vmin = 0.0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label3 + ', q [g/kg]')
    #
    # ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['q_6hrly'])*1e3,
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label1 + ', q [g/kg]')
    #
    # ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['q_6hrly'])*1e3,
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.title(label2 + ', q [g/kg]')
    # plt.xlabel('Day of year')
    #
    # ### -------------------------------
    # ### sonde and ifs data interpolated to um grid (Z<10km)
    # ### ------------------------------
    # ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:]),
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title('Sondes(REGRID), q [g/kg]')
    #
    # ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['q_hrly_UM'][::6])*1e3,
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label3 + '(REGRID), q [g/kg]')
    #
    # ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3,
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label1 + ', q [g/kg]')
    #
    # ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3,
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.xlabel('Day of year')
    # # plt.ylabel('Z [m]')
    # plt.title(label2 + ', q [g/kg]')
    #
    # ### -------------------------------
    # ### model anomalies wrt radiosondes
    # ### ------------------------------
    # ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    # plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:]),
    #     vmin = 0, vmax = qmax)
    # plt.ylim([0,ymax])
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title('Sondes(REGRID),  q [g/kg]')
    #
    # ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    # dat3 = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    # plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
    #     vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.set_cmap('seismic')
    # # plt.ylabel('Z [m]')
    # plt.title(label3 + '(REGRID) - Sondes(REGRID),  q [g/kg]')
    #
    # ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    # dat1 = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    # plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
    #     vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.title(label1 + ' - Sondes(REGRID), q [g/kg]')
    #
    # ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    # dat2 = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_allSondes_UM'][drift[0],:])
    # plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
    #     vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    # plt.ylim([0,ymax])
    # plt.xlim([doy[0],doy[-1]])
    # plt.colorbar()
    # plt.xlabel('Day of year')
    # # plt.ylabel('Z [m]')
    # plt.title(label2 + ' - Sondes(REGRID), q [g/kg]')
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # fileout = '../FIGS/comparisons/QProfiles_REGRID_10km_sondes_metum_ifs_casim-100.png'
    # # plt.savefig(fileout, dpi = 300)
    # plt.show()

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    #### anomaly plot on it's own

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(7,10))

    ymax = 8000
    qmax = 4.0
    yticklist = [0, 3e3, 6e3, 9e3]

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.15,0.82,0.85,0.13])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['q_driftSondes_UM']),
        vmin = 0, vmax = qmax)
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Radiosondes,  q [g/kg]')

    ax  = fig.add_axes([0.15,0.63,0.85,0.13])   # left, bottom, width, height
    dat3 = np.transpose(data3['q_hrly_UM'][::6])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data3['q_anomalies'] = dat3
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + ' - Radiosondes,  q [g/kg]')

    ax  = fig.add_axes([0.15,0.44,0.85,0.13])   # left, bottom, width, height
    dat1 = np.transpose(data1['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data1['q_anomalies'] = dat1
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Radiosondes, q [g/kg]')

    ax  = fig.add_axes([0.15,0.25,0.85,0.13])   # left, bottom, width, height
    dat2 = np.transpose(data2['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data2['q_anomalies'] = dat2
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.yticks(yticklist)    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Radiosondes, q [g/kg]')

    ax  = fig.add_axes([0.15,0.06,0.85,0.13])   # left, bottom, width, height
    dat4 = np.transpose(data4['q_6hrly'][:,data1['universal_height_UMindex']])*1e3 - np.transpose(obs['sondes']['q_driftSondes_UM'])
    data4['q_anomalies'] = dat4
    plt.pcolor(data4['time_6hrly'],data1['universal_height'], dat4,
        vmin = -1.5, vmax = 1.5, cmap=mpl_cm.RdBu_r)
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k')
    plt.ylim([0,ymax])
    plt.yticks(yticklist)
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label4[:-4] + ' - Radiosondes, q [g/kg]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/QProfiles_REGRID_anomOnly_sondes_ifs_ra2m_casim-100_ra2t.png'
    # plt.savefig(fileout, dpi=300)
    plt.show()
    # plt.close()

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################
    SMALL_SIZE = 12
    MED_SIZE = 16
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(11,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)

    ####        all model data share a timestamp
    melt = np.where(data1['time_hrly'][::6] < 240.0)
    freeze = np.where(data1['time_hrly'][::6] >= 240.0)


    plt.subplot(131)
    ax1 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(data1['q_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'darkblue', label = label1)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1),
        np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) - np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data1['q_anomalies'],1) + np.nanstd(data1['q_anomalies'],1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(data2['q_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'mediumseagreen', label = label2)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1),
        np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) - np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(data2['q_anomalies'],1) + np.nanstd(data2['q_anomalies'],1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(data3['q_anomalies'],1),data1['universal_height'],
        '.-' ,color = 'gold', label = label3)
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1),
        np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) - np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(data3['q_anomalies'],1) + np.nanstd(data3['q_anomalies'],1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(data4['q_anomalies'],1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4])
    ax1.fill_betweenx(data1['universal_height'], np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1),
        np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) - np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(data4['q_anomalies'],1) + np.nanstd(data4['q_anomalies'],1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.legend()
    plt.ylim([0,1e4])
    plt.xlabel('Q bias [g/kg]')
    plt.ylabel('Z [m]')
    plt.ylim([0,1e4])
    plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    plt.grid('on')
    plt.title('Total drift')

    plt.subplot(132)
    ax2 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'darkblue', label = label1 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,melt]),1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'mediumseagreen', label = label2 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,melt]),1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'gold', label = label3 + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,melt]),1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,melt]),1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4] + ' median')
    ax2.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,melt]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,melt]),1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,melt]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,melt]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,melt]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.grid('on')
    plt.xlabel('Q bias [g/kg]')
    plt.ylim([0,1e4])
    plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    plt.title('Melt')

    plt.subplot(133)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'darkblue', label = label1 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,freeze]),1),
        color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'mediumseagreen', label = label2 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,freeze]),1),
        color = 'mediumaquamarine', alpha = 0.2)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'gold', label = label3 + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,freeze]),1),
        color = 'navajowhite', alpha = 0.3)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,freeze]),1),data1['universal_height'],
        '.-', color = 'steelblue', label = label4[:-4] + ' median')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,freeze]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,freeze]),1),
        color = 'salmon', alpha = 0.1)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,freeze]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,freeze]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,freeze]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    plt.grid('on')
    plt.xlabel('Q bias [g/kg]')
    plt.ylim([0,1e4])
    plt.xlim([-0.25,0.45])#plt.xlim([-0.05,0.45])
    plt.title('Freeze')

    fileout = '../FIGS/comparisons/QMedianProfiles_metum_ifs_casim-100_erai-glm.svg'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model thetaE profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data1['pressure'][data1['pressure'] == -9999] = np.nan
    data2['pressure'][data2['pressure'] == -9999] = np.nan
    data3['pressure'][data3['pressure'] <= 0] = np.nan
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### calculate equivalent potential temperature
    #### ---------------------------------------------------------------
    data1['theta'], data1['thetaE'] = calcThetaE(data1['temperature'], data1['pressure'], data1['q'])
    data2['theta'], data2['thetaE'] = calcThetaE(data2['temperature'], data2['pressure'], data2['q'])
    data3['theta'], data3['thetaE'] = calcThetaE(data3['temperature'], data3['pressure'], data3['q'])

    obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
        np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3)

    obs['sondes']['thetaE'] = np.transpose(obs['sondes']['thetaE'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'thetaE')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(19,10))

    Tmin = -5
    Tmax = 45
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['epottemp'][:,drift[0]],
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['thetaE_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta_{E}$ [degC]')
    plt.xlabel('Day of year')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['thetaE_hrly_UM'][::6])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['thetaE_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['thetaE_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta_{E}$ [degC]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta_{E}$ [degC]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['thetaE_hrly_UM'][::6]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), $\Theta_{E}$ [K]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['thetaE_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), $\Theta_{E}$ [K]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['thetaE_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['thetaE_allSondes_UM'][drift[0],:])
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -8.0, vmax = 8.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), $\Theta_{E}$ [K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/ThetaEProfiles_REGRID_10km_sondes-calculated_metum_ifs_casim-100.png'
    plt.savefig(fileout)
    plt.show()

def plot_RadiosondesTheta(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    # from scipy.interpolate import interp1d

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting radiosonde and model theta profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data1['pressure'][data1['pressure'] == -9999] = np.nan
    data2['pressure'][data2['pressure'] == -9999] = np.nan
    data3['pressure'][data3['pressure'] <= 0] = np.nan
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### calculate equivalent potential temperature
    #### ---------------------------------------------------------------
    data1['theta'], data1['thetaE'] = calcThetaE(data1['temperature'], data1['pressure'], data1['q'])
    data2['theta'], data2['thetaE'] = calcThetaE(data2['temperature'], data2['pressure'], data2['q'])
    data3['theta'], data3['thetaE'] = calcThetaE(data3['temperature'], data3['pressure'], data3['q'])

    obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
        np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3)

    obs['sondes']['theta'] = np.transpose(obs['sondes']['theta'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'theta')
    print ('')
    print ('Done!')

    print ('')
    print ('Starting radiosonde figure (quite slow!)...:')
    print ('...')

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(19,10))

    Tmin = -5
    Tmax = 45
    ymax = 8000

    ### -------------------------------
    ### original data
    ### ------------------------------
    ax  = fig.add_axes([0.06,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],obs['sondes']['gpsaltitude'][:,drift[0][0]],obs['sondes']['pottemp'][:,drift[0]],
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.set_cmap('viridis')
    plt.ylabel('Z [m]')
    plt.title('Sondes, $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],np.nanmean(data3['height'],0),np.transpose(data3['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label3 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['height'],np.transpose(data1['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.06,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data2['height'],np.transpose(data2['theta_6hrly'])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta$ [degC]')
    plt.xlabel('Day of year')

    ### -------------------------------
    ### sonde and ifs data interpolated to um grid (Z<10km)
    ### ------------------------------
    ax  = fig.add_axes([0.38,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.54,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data3['time_6hrly'],data1['universal_height'],np.transpose(data3['theta_hrly_UM'][::6])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.3,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data1['time_6hrly'],data1['universal_height'],np.transpose(data1['theta_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ', $\Theta$ [degC]')

    ax  = fig.add_axes([0.38,0.06,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(data2['time_6hrly'],data1['universal_height'],np.transpose(data2['theta_6hrly'][:,data1['universal_height_UMindex']])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ', $\Theta$ [degC]')

    ### -------------------------------
    ### model anomalies wrt radiosondes
    ### ------------------------------
    ax  = fig.add_axes([0.7,0.78,0.3,0.17])   # left, bottom, width, height
    plt.pcolor(obs['sondes']['doy_drift'],data1['universal_height'],np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])-273.15,
        vmin = Tmin, vmax = Tmax)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title('Sondes(REGRID), $\Theta$ [degC]')

    ax  = fig.add_axes([0.7,0.54,0.3,0.17])   # left, bottom, width, height
    dat3 = np.transpose(data3['theta_hrly_UM'][::6]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data3['time_6hrly'], data1['universal_height'], dat3,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.set_cmap('seismic')
    # plt.ylabel('Z [m]')
    plt.title(label3 + '(REGRID) - Sondes(REGRID), $\Theta$ [K]')

    ax  = fig.add_axes([0.7,0.3,0.3,0.17])   # left, bottom, width, height
    dat1 = np.transpose(data1['theta_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data1['time_6hrly'],data1['universal_height'], dat1,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    # plt.ylabel('Z [m]')
    plt.title(label1 + ' - Sondes(REGRID), $\Theta$ [K]')

    ax  = fig.add_axes([0.7,0.06,0.3,0.17])   # left, bottom, width, height
    dat2 = np.transpose(data2['theta_6hrly'][:,data1['universal_height_UMindex']]) - np.transpose(obs['sondes']['theta_allSondes_UM'][drift[0],:])
    plt.pcolor(data2['time_6hrly'],data1['universal_height'], dat2,
        vmin = -6.0, vmax = 6.0, cmap=mpl_cm.RdBu_r)
    plt.ylim([0,ymax])
    plt.xlim([doy[0],doy[-1]])
    plt.colorbar()
    plt.xlabel('Day of year')
    # plt.ylabel('Z [m]')
    plt.title(label2 + ' - Sondes(REGRID), $\Theta$ [K]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = '../FIGS/comparisons/ThetaProfiles_REGRID_10km_sondes-calculated_metum_ifs_casim-100.png'
    plt.savefig(fileout)
    plt.show()

def period_Selection(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5):

    '''
    Must be run AFTER plot_RadiosondesQ and plot_RadiosondesTemperature
    '''

    ### all model data share a timestamp
    p3 = np.where(np.logical_and(data1['time_hrly'][::6] >= doy[0], data1['time_hrly'][::6] < 230.0))
    p4 = np.where(np.logical_and(data1['time_hrly'][::6] >= 230.0, data1['time_hrly'][::6] < 240.0))
    p5 = np.where(np.logical_and(data1['time_hrly'][::6] >= 240.0, data1['time_hrly'][::6] < 247.0))
    p6 = np.where(np.logical_and(data1['time_hrly'][::6] >= 247.0, data1['time_hrly'][::6] < 251.0))
    # p7 = np.where(np.logical_and(data1['time_hrly'][::6] >= 251.0, data1['time_hrly'][::6] < 255.5))
    # p8 = np.where(data1['time_hrly'][::6] >= 255.5)

    ##################################################
    ##################################################
    #### 	P3 and P4
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.figure(figsize=(6.6,10))
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
            hspace = 0.22, wspace = 0.19)

    plt.subplot(221)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p3]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p3]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)

    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax3.set_yticklabels([])
    ax3.set_yticks(axminor, minor = True)
    ax3.grid(which = 'major', alpha = 0.5)
    plt.xlim([-5.0,2.5])
    plt.xlabel('T bias [K]')

    plt.subplot(222)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p3]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p3]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p3]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p3]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p3]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p3]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,p3]),1),data1['universal_height'],'.-', color = 'grey', label = label1 + ' median', zorder = 3)

    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9e3])
    plt.yticks(axmajor)
    ax3.set_yticklabels([])
    ax3.set_yticks(axminor, minor = True)
    ax3.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.6,0.6])
    plt.xticks([-0.6,-0.3,0,0.3,0.6])

    # plt.subplot(223)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p4]),1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p4]),1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    #
    # plt.grid('on')
    # plt.xlim([-5.0,2.5])
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlabel('T bias [K]')

    # plt.subplot(224)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data1['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p4]),1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p4]),1),
    #     np.nanmedian(np.squeeze(data3['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p4]),1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p4]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p4]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p4]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p4]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    #
    # # plt.grid('on')
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])
    #
    # fileout = '../FIGS/comparisons/Temp-SpHumMedianProfiles_metum_ifs_casim-100_ra2t_periodSelection-p3-p4_wSTDEV_newColours_fixedRA2T.svg'
    # # plt.savefig(fileout)
    # plt.show()

    ##################################################
    ##################################################
    #### 	P5 and P6
    ##################################################
    ##################################################

    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=LARGE_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.figure(figsize=(6.6,10))
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
    #         hspace = 0.22, wspace = 0.19)
    #
    # plt.subplot(221)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p5]),1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p5]),1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    # plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    #
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-5.0,2.5])
    # plt.xlabel('T bias [K]')
    #
    # plt.subplot(222)
    # ax3 = plt.gca()
    # plt.plot([0,0], [0,1e4], '--', color='grey')
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data1['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p5]),1),
    #     color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1),
    #     color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1),
    #     color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p5]),1),
    #     np.nanmedian(np.squeeze(data3['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p5]),1),
    #     color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p5]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p5]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p5]),1), data1['universal_height'],
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    # plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    # plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    # plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p5]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    #
    # plt.xlabel('q bias [g kg$^{-1}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax3.set_yticklabels([])
    # ax3.set_yticks(axminor, minor = True)
    # ax3.grid(which = 'major', alpha = 0.5)
    # plt.xlim([-0.6,0.6])
    # plt.xticks([-0.6,-0.3,0,0.3,0.6])

    plt.subplot(223)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p6]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data1['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data1['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p6]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data3['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data3['temp_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    # plt.plot(np.nanmedian(np.squeeze(data5['temp_anomalies'][:,p6]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    plt.grid('on')
    plt.xlim([-5.0,2.5])
    plt.ylim([0,9e3])
    plt.yticks(axmajor)
    ax3.set_yticklabels([])
    ax3.set_yticks(axminor, minor = True)
    ax3.grid(which = 'major', alpha = 0.5)
    plt.xlabel('T bias [K]')

    plt.subplot(224)
    ax3 = plt.gca()
    plt.plot([0,0], [0,1e4], '--', color='grey')
    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data1['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data1['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p6]),1),
        color = 'blue', alpha = 0.05)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data1['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data1['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'darkblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1),
        color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data4['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'steelblue', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1),
        color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data2['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax3.fill_betweenx(data1['universal_height'], np.nanmedian(np.squeeze(data3['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p6]),1),
        np.nanmedian(np.squeeze(data3['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p6]),1),
        color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p6]),1) - np.nanstd(np.squeeze(data3['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p6]),1) + np.nanstd(np.squeeze(data3['q_anomalies'][:,p6]),1), data1['universal_height'],
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmedian(np.squeeze(data1['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'darkblue', label = label1 + ' median', zorder = 3)
    plt.plot(np.nanmedian(np.squeeze(data4['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'steelblue', label = label4[:-4] + ' median', zorder = 2)
    plt.plot(np.nanmedian(np.squeeze(data2['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'mediumseagreen', label = label2 + ' median', zorder = 1)
    plt.plot(np.nanmedian(np.squeeze(data3['q_anomalies'][:,p6]),1),data1['universal_height'],'.-', color = 'gold', label = label3 + ' median', zorder = 4)
    # plt.plot(np.nanmedian(np.squeeze(data5['q_anomalies'][:,p6]),1),data1['universal_height'],'.-' ,linewidth = 3, markersize = 8, color = 'grey', label = label5, zorder = 1)

    # plt.grid('on')
    plt.xlabel('q bias [g kg$^{-1}$]')
    plt.ylim([0,9e3])
    plt.yticks(axmajor)
    ax3.set_yticklabels([])
    ax3.set_yticks(axminor, minor = True)
    ax3.grid(which = 'major', alpha = 0.5)
    plt.xlim([-0.6,0.6])
    plt.xticks([-0.6,-0.3,0,0.3,0.6])

    fileout = '../FIGS/comparisons/Temp-SpHumMedianProfiles_metum_ifs_casim-aeroprof_ra2t_periodSelection-p3-p6_wSTDEV_newColours_fixedRA2T.svg'
    # plt.savefig(fileout)
    plt.show()

    Zindex1 = np.where(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1) < 0)
    Zindex3 = np.where(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1) < 0)
    print ('Z(universal) = ')
    # print (data1['universal_height'][Zindex3])
    print (data1['universal_height'])

    print ('UM_RA2M = ')
    print (np.round(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1),2))

    print ('UM_CASIM-100 = ')
    print (np.round(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1),2))

    print ('ECMWF_IFS = ')
    print (np.round(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1),2))

    print ('UM_RA2T = ')
    print (np.round(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1),2))

    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data1['temp_anomalies'][:,p3]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data2['temp_anomalies'][:,p3]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data3['temp_anomalies'][:,p3]),1),2)))
    print (np.nanmin(np.round(np.nanmedian(np.squeeze(data4['temp_anomalies'][:,p3]),1),2)))


def reGrid_Sondes(data1, data2, data3, data4, data5, obs, doy, ifs_flag, var):

    from scipy.interpolate import interp1d

    ### 6-hourly time binning for model
    ### um['time'][:24:6].data
    ###     BUT there is a problem since we have 25 timesteps (i.e. [24] == [25])
    ###     need to pick out where we have a repeated time value, then remove it so
    ###     that the time array can be indexed easily

    #### ---------------------------------------------------------------
    ### build list of variables names wrt input data [OBS, UM, CASIM, IFS]
    #### ---------------------------------------------------------------
    if var == 'temp':
        varlist = ['temperature','temperature','temperature','temperature', 'temperature', 'temperature']
    elif var == 'thetaE':
        # varlist = ['epottemp','thetaE','thetaE','thetaE']     # use sonde file's epottemp
        varlist = ['thetaE','thetaE','thetaE','thetaE', 'thetaE', 'thetaE']         # use sonde calculated thetaE
    elif var == 'theta':
        # varlist = ['pottemp','theta','theta','theta']     # use sonde file's pottemp
        varlist = ['theta','theta','theta','theta','theta','theta']         # use sonde calculated theta
    elif var == 'q':
        varlist = ['sphum','q','q','q','q','q']

    # ### stop double counting of 0000 and 2400 from model data
    # temp = np.zeros([len(data1['time'])])
    # for i in range(0, len(temp)-1):
    #     if data1['time'][i] == data1['time'][i+1]:
    #         continue
    #     else:
    #         temp[i] = data1['time'][i]
    # ii = np.where(temp != 0.0)      ### picks out where data are non-zero

    #### ---------------------------------------------------------------
    #### save hourly temperature model profiles (using the ii index defined by the time indices)
    #### ---------------------------------------------------------------
    data1[var + '_hrly'] = np.squeeze(data1[varlist[1]][data1['hrly_flag'],:])
    data2[var + '_hrly'] = np.squeeze(data2[varlist[2]][data2['hrly_flag'],:])
    data3[var + '_hrly'] = np.squeeze(data3[varlist[3]][data3['hrly_flag'],:])
    data4[var + '_hrly'] = np.squeeze(data4[varlist[4]][data4['hrly_flag'],:])
    data5[var + '_hrly'] = np.squeeze(data5[varlist[5]][data5['hrly_flag'],:])

    # print(data2[var + '_hrly'].shape)

    #### ---------------------------------------------------------------
    #### explicitly save 6-hourly temperature model profiles and time binning for ease
    #### ---------------------------------------------------------------
    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_6hrly'] = data1['time_hrly'][::6]
    data2['time_6hrly'] = data2['time_hrly'][::6]
    data3['time_6hrly'] = data3['time_hrly'][::6]
    data4['time_6hrly'] = data4['time_hrly'][::6]
    data5['time_6hrly'] = data5['time_hrly'][::6]
    data1[var + '_6hrly'] = data1[var + '_hrly'][::6]
    data2[var + '_6hrly'] = data2[var + '_hrly'][::6]
    data3[var + '_6hrly'] = data3[var + '_hrly'][::6]
    data4[var + '_6hrly'] = data4[var + '_hrly'][::6]
    data5[var + '_6hrly'] = data5[var + '_hrly'][::6]

    print(data2[var + '_6hrly'].shape)

    #### ---------------------------------------------------------------
    #### index to only look at altitudes <10km
    #### ---------------------------------------------------------------
    iTim = 0        ### initialised
    iObs = np.where(obs['sondes']['gpsaltitude'][:,iTim] <= 11000)
    iUM = np.where(data1['height'] <= 11000)
    iGLM = np.where(data5['height'] <= 11000)
    if ifs_flag == True:
        iIFS = np.where(data3['height'][iTim,:] <= 11000)
    else:
        iIFS = np.where(data3['height'] <= 11000)

    #### ---------------------------------------------------------------
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    if ifs_flag == True:
        data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems
            ####        further on
        data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### START INTERPOLATION
    #### ---------------------------------------------------------------
    if ifs_flag == True:
        print ('')
        print ('Defining IFS temperature profile as a function:')
        print ('using ifs.height[i,:] to define temperature profiles...')
        data3[var + '_hrly_UM'] = np.zeros([np.size(data3['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
        for iTim in range(0,np.size(data3['time_hrly'],0)):
            # print (iTim)
            iIFSind = np.where(data3['height_hrly'][iTim,:] <= 11000)
            if np.all(data3['height_hrly'][iTim,:] == 0.0):
                data3[var + '_hrly_UM'][iTim,:] = np.nan
            else:
                fnct_IFS = interp1d(np.squeeze(data3['height_hrly'][iTim,iIFSind]), np.squeeze(data3[var + '_hrly'][iTim,iIFSind]))
                data3[var + '_hrly_UM'][iTim,:] = fnct_IFS(data1['height'][iUM[0][3:]].data)
        print ('...')
        print ('IFS(UM Grid) function worked!')
        print (var + ' IFS data now on UM vertical grid')
        print ('*****')
        ### assign for easier indexing later
        data3[var + '_6hrly_UM'] = data3[var + '_hrly_UM'][::6,:]
        data3[var + '_6hrly'] = data3[var + '_hrly'][::6,:]
        data3['height_6hrly'] = data3['height_hrly'][::6,:]  ### need to explicitly save since height coord changes at each timedump
    else:
        data3['universal_height'] = data3['height'][iUM[0][3:]]
        data3['universal_height_UMindex'] = iUM[0][3:]
        data3['height_6hrly'] = np.zeros([np.size(data3['time_6hrly']), np.size(data3['universal_height_UMindex'])])
        for i in range(0, len(data3['time_6hrly'])): data3['height_6hrly'][i,:] = data3['universal_height'][:]
        data3[var + '_6hrly_UM'] = data3[var + '_hrly'][::6]


    #### INTERPOLATION TESTING:
    # print (data3['temp_hrly_UM'].shape),'radr_refl'
    # print (data3['time_hrly'][::6].shape)
    # print (data1['temp_hrly'][:,iUM[0][3:]].shape)
    # print (data1['time_hrly'][::6].shape)
    # for i in range(0, np.size(data3['temp_6hrly_UM'],0)):
    #     fig = plt.figure()
    #     plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], label = 'interpd')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), label = 'height indexed')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height'][0,iIFS]), label = 'height0')
    #     plt.title('IFS test ' + str(data3['time_6hrly'][i]))
    #     plt.legend()
    #     plt.savefig('../FIGS/regrid/IFS_test_doy' + str(data3['time_6hrly'][i]) + '.png')
    #     if i == 0:
    #         plt.show()
    #     else:
    #         plt.close()

    print ('')
    print ('Defining Sonde temperature profile as a function for the UM:')
    obs['sondes'][var + '_allSondes_UM'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(obs['sondes']['doy'],0)):
        # print 'iTim = ', str(iTim)
        fnct_Obs = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
        obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_Obs(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('Sonde(UM Grid) function worked!')
    print ('All ' + var + ' sonde data now on UM vertical grid.')
    print ('*****')
    #
    # print ('')
    # print ('Defining Sonde temperature profile as a function for the IFS:')
    # obs['sondes'][var + '_allSondes_IFS'] = np.zeros([np.size(obs['sondes']['doy'],0),len(data1['height'][0,iIFS])])
    # for iTim in range(0,np.size(obs['sondes']['doy'],0)):
    #     # print 'iTim = ', str(iTim)
    #     iIFS = np.where(data3['height'][iTim,:] <= 11000)
    #     fnct_ObsIFS = interp1d(np.squeeze(obs['sondes']['gpsaltitude'][iObs,iTim]), np.squeeze(obs['sondes'][varlist[0]][iObs,iTim]))
    #     obs['sondes'][var + '_allSondes_UM'][iTim,:] = fnct_ObsIFS(data3['height'][iTim,iIFS])
    # print ('...')
    # print ('Sonde(IFS Grid) function worked!')
    # print ('All ' + var + ' sonde data now on IFS_DATA vertical grid.')
    # print ('*****')

    print ('')
    print ('Defining GLM profile as a function:')
    data5[var + '_hrly_UM'] = np.zeros([np.size(data5['time_hrly'],0),len(data1['height'][iUM[0][3:]])])
    for iTim in range(0,np.size(data5['time_hrly'],0)):
        # print (iTim)
        fnct_GLM = interp1d(np.squeeze(data5['height'][iGLM]), np.squeeze(data5[var + '_hrly'][iTim,iGLM]))
        data5[var + '_hrly_UM'][iTim,:] = fnct_GLM(data1['height'][iUM[0][3:]].data)
    print ('...')
    print ('GLM(UM Grid) function worked!')
    print (var + ' GLM data now on UM(lam) vertical grid')
    print ('*****')
    ### assign for easier indexing later
    data5[var + '_6hrly_UM'] = data5[var + '_hrly_UM'][::6,:]
    # print (np.size(data5[var + '_6hrly_UM'],0))

    #### ---------------------------------------------------------------
    #### ONLY LOOK AT SONDES FROM THE DRIFT
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['sondes']['doy'] >= 225.9, obs['sondes']['doy'] <= 258.0))
    subset = np.where(np.logical_and(obs['sondes']['doy'] >= doy[0], obs['sondes']['doy'] <= doy[-1] + 0.05))
    # drift = np.where(np.logical_and(obs['sondes']['doy'] >= doy[0], obs['sondes']['doy'] <= doy[-1] + 0.05))

    # print (obs['sondes']['doy'][drift[0]])
    # print (obs['sondes']['doy'][drift[-1]])

    ### save in dict for ease
    obs['sondes']['doy_drift'] = obs['sondes']['doy'][drift]
    obs['sondes']['drift'] = drift
    obs['sondes'][var + '_driftSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][drift[0],:]

    # print (obs['sondes'][var + '_driftSondes_UM'].shape)
    # print (data2[var + '_6hrly'].shape)
    # print (data4[var + '_6hrly'].shape)

    ### save subset of drift in dict for ease
    obs['sondes']['doy_subset'] = obs['sondes']['doy'][subset]
    obs['sondes']['subset'] = subset
    obs['sondes'][var + '_subsetSondes_UM'] = obs['sondes'][var + '_allSondes_UM'][subset[0],:]

    # #### INTERPOLATION TESTING - IFS + SONDE + UM_RA2M:
    # print (obs['sondes']['doy_drift'].shape)
    # print (obs['sondes']['temp_allSondes_UM'][drift[0],:].shape)
    # if var == 'temp':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['temperature'][iObs,drift[0][i]]) + 273.15,np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['temp_driftSondes_UM'][i,:] + 273.15,data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'gold', label = 'ifs-Zindexed')
    #         plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = 'gold', label = 'ifs-interpd')
    #         plt.plot(data1['temp_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2['temp_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.plot(data4['temp_6hrly'][i,iUM[0][3:]], data4['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2t')
    #         plt.plot(np.squeeze(data5['temp_6hrly'][i,iGLM]), data5['height'][iGLM], '--', color = 'grey', label = 'um_glm')
    #         plt.plot(data5['temp_6hrly_UM'][i,:], data1['height'][iUM[0][3:]], color = 'grey', label = 'um_glm-interpd')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)) + '\n Model time = ' + str(data1['time_6hrly'][i]))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid_v3/REGRID_Ttest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'q':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['mr'][iObs,drift[0][i]]), np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes'][var + '_driftSondes_UM'][i,:], data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3[var + '_6hrly'][i,iIFS])*1e3,np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'gold', label = 'ifs-Zindexed')
    #         plt.plot(data3[var + '_6hrly_UM'][i,:]*1e3,data1['height'][iUM[0][3:]], color = 'gold', label = 'ifs-interpd')
    #         plt.plot(data1[var + '_6hrly'][i,iUM[0][3:]]*1e3, data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2[var + '_6hrly'][i,iUM[0][3:]]*1e3, data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.plot(data4[var + '_6hrly'][i,iUM[0][3:]]*1e3, data4['height'][iUM[0][3:]], color = 'steelblue', label = 'um_ra2t')
    #         plt.plot(np.squeeze(data5['q_6hrly'][i,iGLM])*1e3, data5['height'][iGLM], '--', color = 'grey', label = 'um_glm')
    #         plt.plot(data5['q_6hrly_UM'][i,:]*1e3, data1['height'][iUM[0][3:]], color = 'grey', label = 'um_glm-interpd')
    #         plt.title('REGRID test ' + str(np.round(obs['sondes']['doy_drift'][i],2)) + '\n Model time = ' + str(data1['time_6hrly'][i]))
    #         plt.legend()
    #         plt.savefig('../FIGS/regrid_v3/REGRID_Qtest_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()
    # elif var == 'thetaE':
    #     for i in range(0, np.size(obs['sondes']['doy_drift'])):
    #         plt.plot(np.squeeze(obs['sondes']['thetaE'][iObs,drift[0][i]]),np.squeeze(obs['sondes']['gpsaltitude'][iObs,drift[0][i]]), '--', color = 'k', label = 'sonde-original')
    #         plt.plot(obs['sondes']['thetaE_driftSondes_UM'][i,:],data1['height'][iUM[0][3:]], color = 'k', label = 'sonde-interpd')
    #         plt.plot(np.squeeze(data3['thetaE_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), '--', color = 'gold', label = 'ifs-Zindexed')
    #         plt.plot(data3['thetaE_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], color = 'gold', label = 'ifs-interpd')
    #         plt.plot(data1['thetaE_6hrly'][i,iUM[0][3:]], data1['height'][iUM[0][3:]], color = 'darkblue', label = 'um_ra2m')
    #         plt.plot(data2['thetaE_6hrly'][i,iUM[0][3:]], data2['height'][iUM[0][3:]], color = 'mediumseagreen', label = 'um_casim-100')
    #         plt.title('REGRID test DOY ' + str(np.round(obs['sondes']['doy_drift'][i],2)))
    #         plt.xlabel('$\Theta_{E}$ [K]')
    #         plt.ylabel('Z [m]')
    #         plt.ylim([0,3000])
    #         plt.xlim([260,320])
    #         plt.legend()
    #         plt.savefig('../FIGS/inversionIdent/REGRID_ThetaE_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
    #         if i == 0:
    #             plt.show()
    #         else:
    #             plt.close()

    #### ---------------------------------------------------------------
    #### make some dictionary assignments for use later
    #### ---------------------------------------------------------------
    data1['universal_height'] = data1['height'][iUM[0][3:]]
    data1['universal_height_UMindex'] = iUM[0][3:]

    # print (data2[var + '_6hrly'][:,data1['universal_height_UMindex']].shape)
    # print (data4[var + '_6hrly'][:,data1['universal_height_UMindex']].shape)

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    np.save('working_data1',data1)
    np.save('working_data2',data2)
    np.save('working_data3',data3)
    np.save('working_dataObs',obs['sondes'])
    # outfiles = write_reGrid(data1, data2, data3, obs, var)

    return data1, data2, data3, data4, data5, obs, drift

def write_reGrid(data1, data2, data3, obs, var):

    #################################################################
    ## Write regridded data to netCDF
    #################################################################

    from netCDF4 import num2date, date2num
    import time
    from datetime import datetime, timedelta

    print ('******')
    print ('')
    print ('Writing out data:')
    print ('')

    outfiles = ['REGRID-' + var + '-SONDES.nc','REGRID-' + var + '-UM_RA2M.nc','REGRID-' + var + '-UM_CASIM-100.nc','REGRID-' + var + '-ECMWF_IFS.nc']

    ###################################
    ## Open File
    ###################################
    nc0 = Dataset(outfiles[0], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc0.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc0.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc0.createDimension('time', np.size(obs['sondes']['doy_drift']))
    height = nc0.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    tempvar = np.zeros(np.size(obs['sondes']['doy_drift']))
    tempvar[:] = obs['sondes']['doy_drift'][:]
    times = nc0.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = tempvar[:]

    #### height
    height = nc0.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat0 = nc0.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat0.scale_factor = float(1)
    dat0.add_offset = float(0)
    if var == 'temp':
        dat0.units = 'K'
        dat0.long_name = 'temperature'
        dat0[:,:] = obs['sondes'][var + '_driftSondes_UM'][:,:] + 273.15
    elif var == 'q':
        dat0.units = 'g/kg'
        dat0.long_name = 'water vapour mixing ratio'
        dat0[:,:] = obs['sondes'][var + '_driftSondes_UM'][:,:]

    nc0.title = 'Radiosonde interpolated ' + var + ' data for the AO2018 drift period.'
    nc0.description = var + ' data up to 10 km, interpolated on to the UM vertical grid.'
    nc0.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc0.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc0.comment = 'Revision no. 0: Preliminary data.'
    nc0.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc0.close()


    ###################################
    ## Open File
    ###################################
    nc1 = Dataset(outfiles[1], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc1.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc1.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc1.createDimension('time', np.size(data1['time_6hrly']))
    height = nc1.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc1.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data1['time_6hrly'][:]

    #### height
    height = nc1.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat1 = nc1.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat1.scale_factor = float(1)
    dat1.add_offset = float(0)
    if var == 'temp':
        dat1.units = 'K'
        dat1.long_name = 'temperature'
        dat1[:,:] = data1[var + '_6hrly'][:,data1['universal_height_UMindex']]
    elif var == 'q':
        dat1.units = 'g/kg'
        dat1.long_name = 'water vapour mixing ratio'
        dat1[:,:] = data1[var + '_6hrly'][:,data1['universal_height_UMindex']]*1e3

    nc1.title = 'UM_RA2M ' + var + ' data for the AO2018 drift period.'
    nc1.description = var + ' data up to 10 km, referencing UM vertical grid.'
    nc1.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc1.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc1.comment = 'Revision no. 0: Preliminary data.'
    nc1.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc1.close()

    ###################################
    ## Open File
    ###################################
    nc2 = Dataset(outfiles[2], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc2.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc2.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc2.createDimension('time', np.size(data2['time_6hrly']))
    height = nc2.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc2.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data2['time_6hrly'][:]

    #### height
    height = nc2.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat2 = nc2.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat2.scale_factor = float(1)
    dat2.add_offset = float(0)
    if var == 'temp':
        dat2.units = 'K'
        dat2.long_name = 'temperature'
        dat2[:,:] = data2[var + '_6hrly'][:,data1['universal_height_UMindex']]
    elif var == 'q':
        dat2.units = 'g/kg'
        dat2.long_name = 'water vapour mixing ratio'
        dat2[:,:] = data2[var + '_6hrly'][:,data1['universal_height_UMindex']]*1e3

    nc2.title = 'UM_CASIM-100 ' + var + ' data for the AO2018 drift period.'
    nc2.description = var + ' data up to 10 km, referencing UM vertical grid.'
    nc2.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc2.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc2.comment = 'Revision no. 0: Preliminary data.'
    nc2.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc2.close()


    ###################################
    ## Open File
    ###################################
    nc3 = Dataset(outfiles[3], 'w', format ='NETCDF4_CLASSIC')
    print ('')
    print (nc3.file_format)
    print ('')

    ###################################
    ## Switch off automatic filling
    ###################################
    nc3.set_fill_off()

    ###################################
    ## Data dimensions
    ####################################
    times = nc3.createDimension('time', np.size(data3['time_6hrly']))
    height = nc3.createDimension('height', np.size(data1['universal_height']))

    ###################################
    ## Dimensions variables
    ###################################
    #### forecast_period
    times = nc3.createVariable('time', np.float64, ('time',), fill_value='-9999')
    times.scale_factor = float(1)
    times.add_offset = float(0)
    times.comment = 'DOY in AO2018 drift.'
    times.units = 'hours'
    times.long_name = 'time'
    times[:] = data3['time_6hrly'][:]

    #### height
    height = nc3.createVariable('height', np.float64, ('height',), fill_value='-9999')
    height.scale_factor = float(1)
    height.add_offset = float(0)
    height.comment = 'Height interpolated on to UM vertical grid (where appropriate)'
    height.units = 'm'
    height.long_name = 'height'
    height[:] = data1['universal_height'][:]

    ###################################
    ## Create DIAGNOSTICS
    ###################################
    dat3 = nc3.createVariable(var, np.float64, ('time','height',), fill_value='-9999')
    dat3.scale_factor = float(1)
    dat3.add_offset = float(0)
    if var == 'temp':
        dat3.units = 'K'
        dat3.long_name = 'temperature'
        dat3[:,:] = data3[var + '_6hrly_UM'][:,:]
    elif var == 'q':
        dat3.units = 'g/kg'
        dat3.long_name = 'water vapour mixing ratio'
        dat3[:,:] = data3[var + '_6hrly_UM'][:,:]*1e3

    nc3.title = 'ECMWF_IFS interpolated ' + var + ' data for the AO2018 drift period.'
    nc3.description = var + ' data up to 10 km, interpolated on to the UM vertical grid.'
    nc3.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ' by Gillian Young <G.Young1@leeds.ac.uk> using Python (netCDF4).'
    nc3.project = 'Arctic Ocean 2018 (AO2018) expedition.'
    nc3.comment = 'Revision no. 0: Preliminary data.'
    nc3.institution = 'University of Leeds.'

    ###################################
    ## Write out file
    ###################################
    nc3.close()

    return outfiles

def checkInvbaseBelow(invbaseID, thetaEDiff, dthresh):

    if np.nanmax(thetaEDiff) >= 0.0:
        if int(invbaseID)-1 > 0:        ### so we don't get a negative index
            if thetaEDiff[int(invbaseID)-1] > dthresh:
                invbaseID = int(invbaseID) - 1

    return invbaseID

def inversionIdent(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    print ('******')
    print ('')
    print ('Identifying stable layers from thetaE profiles:')
    print ('')

    #### change matlab time to doy
    obs['sondes']['doy'] = calcTime_Mat2DOY(np.squeeze(obs['sondes']['mday']))

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    #### set flagged values to nans
    data1['temperature'][data1['temperature'] == -9999] = np.nan
    data2['temperature'][data2['temperature'] == -9999] = np.nan
    data3['temperature'][data3['temperature'] <= 0] = np.nan
    data1['pressure'][data1['pressure'] == -9999] = np.nan
    data2['pressure'][data2['pressure'] == -9999] = np.nan
    data3['pressure'][data3['pressure'] <= 0] = np.nan
    data1['q'][data1['q'] == -9999] = np.nan
    data2['q'][data2['q'] == -9999] = np.nan
    data3['q'][data3['q'] <= 0] = np.nan

    #### ---------------------------------------------------------------
    #### calculate equivalent potential temperature
    #### ---------------------------------------------------------------
    data1['theta'], data1['thetaE'] = calcThetaE(data1['temperature'], data1['pressure'], data1['q'])
    data2['theta'], data2['thetaE'] = calcThetaE(data2['temperature'], data2['pressure'], data2['q'])
    data3['theta'], data3['thetaE'] = calcThetaE(data3['temperature'], data3['pressure'], data3['q'])

    obs['sondes']['theta'], obs['sondes']['thetaE'] = calcThetaE(np.transpose(obs['sondes']['temperature'])+273.15,
        np.transpose(obs['sondes']['pressure'])*1e2, np.transpose(obs['sondes']['mr'])/1e3)

    obs['sondes']['thetaE'] = np.transpose(obs['sondes']['thetaE'])         ### for consistency with original sonde dimensions

    #### ---------------------------------------------------------------
    #### save out working data for debugging
    #### ---------------------------------------------------------------
    # np.save('working_data1',data1)
    # np.save('working_data3',data3)
    # np.save('working_dataObs',obs['sondes'])

    # #### ---------------------------------------------------------------
    # #### Save out line profiles for reference
    # #### ---------------------------------------------------------------
    # for i in range(0, np.size(data3['temp_6hrly_UM'],0)):
    #     fig = plt.figure()
    #     plt.plot(data3['temp_6hrly_UM'][i,:],data1['height'][iUM[0][3:]], label = 'interpd')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height_6hrly'][i,iIFS]), label = 'height indexed')
    #     plt.plot(np.squeeze(data3['temp_6hrly'][i,iIFS]),np.squeeze(data3['height'][0,iIFS]), label = 'height0')
    #     plt.title('IFS test ' + str(data3['time_6hrly'][i]))
    #     plt.legend()
    #     plt.savefig('../FIGS/regrid/IFS_test_doy' + str(data3['time_6hrly'][i]) + '.png')
    #     if i == 0:
    #         plt.show()
    #     else:
    #         plt.close()

    #################################################################
    ## save data into temp variables to allow subsampling
    #################################################################
    bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    if ifs_flag == True:
        bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    else:
        bldepth3 = data3['bl_depth'][data3['hrly_flag']]

    #### ---------------------------------------------------------------
    #### re-grid sonde and IFS data to UM vertical grid <10km
    #### ---------------------------------------------------------------

    print ('...')
    print ('Re-gridding sonde and ifs data...')
    print ('')
    data1, data2, data3, obs, drift = reGrid_Sondes(data1, data2, data3, obs, doy, 'thetaE')
    print ('')
    print ('Done!')

    #### ---------------------------------------------------------------
    #### calculate dThetaE and dZ
    #### ---------------------------------------------------------------
    obs['sondes']['dThetaE'] = obs['sondes']['thetaE_driftSondes_UM'][:,1:] - obs['sondes']['thetaE_driftSondes_UM'][:,0:-1]
    data1['dThetaE'] = data1['thetaE_6hrly'][:,data1['universal_height_UMindex'][1:]].data - data1['thetaE_6hrly'][:,data1['universal_height_UMindex'][0:-1]].data
    data2['dThetaE'] = data2['thetaE_6hrly'][:,data1['universal_height_UMindex'][1:]].data - data2['thetaE_6hrly'][:,data1['universal_height_UMindex'][0:-1]].data
    data3['dThetaE'] = data3['thetaE_6hrly_UM'][:,1:] - data3['thetaE_6hrly_UM'][:,0:-1]
    data1['dZ'] = data1['universal_height'][1:27].data - data1['universal_height'][:26].data

    #### ---------------------------------------------------------------
    #### calculate dThetaE/dZ
    #### ---------------------------------------------------------------
    obs['sondes']['dThetaEdZ'] = obs['sondes']['dThetaE'][:,:26] / data1['dZ']
    data1['dThetaEdZ'] = data1['dThetaE'][:,:26] / data1['dZ']
    data2['dThetaEdZ'] = data2['dThetaE'][:,:26] / data1['dZ']
    data3['dThetaEdZ'] = data3['dThetaE'][:,:26] / data1['dZ']

    #### ---------------------------------------------------------------
    #### build bespoke thetaE arrays for casim-100 and ra2m on universal  (for ease)
    #### ---------------------------------------------------------------
    data1['thetaE_6hrly_UM'] = data1['thetaE_6hrly'][:,data1['universal_height_UMindex']].data
    data2['thetaE_6hrly_UM'] = data2['thetaE_6hrly'][:,data1['universal_height_UMindex']].data

    #### ---------------------------------------------------------------
    #### choose "inversion" gradient dthreshold (K)
    ####        dthresh = threshold for top of decoupled surface layer (smaller than sthresh)
    ####        sthresh = threshold for secondary inversion above main inversion layer (used for main and secondary inversion)
    #### ---------------------------------------------------------------
    dthresh = 1.0
    sthresh = 2.0

    #### ---------------------------------------------------------------
    #### save inversion positions
    #### ---------------------------------------------------------------
    #### initialise arrays
    lt3000 = np.where(data1['universal_height'] < 3000)
    obs['sondes']['thetaE_invbaseID'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['thetaE_invbaseID'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['thetaE_invbaseID'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['thetaE_invbaseID'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['thetaE_invbaseID'][:] = np.nan           ## fill with nans to account for missing files when populating
    obs['sondes']['dThetaEdZ_invbaseID'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['dThetaEdZ_invbaseID'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['dThetaEdZ_invbaseID'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['dThetaEdZ_invbaseID'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['dThetaEdZ_invbaseID'][:] = np.nan           ## fill with nans to account for missing files when populating
    obs['sondes']['thetaE_invbase'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['thetaE_invbase'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['thetaE_invbase'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['thetaE_invbase'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['thetaE_invbase'][:] = np.nan           ## fill with nans to account for missing files when populating
    obs['sondes']['dThetaEdZ_2ndinvID'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['dThetaEdZ_2ndinvID'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['dThetaEdZ_2ndinvID'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['dThetaEdZ_2ndinvID'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['dThetaEdZ_2ndinvID'][:] = np.nan           ## fill with nans to account for missing files when populating
    obs['sondes']['dThetaEdZ_decoupID'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['dThetaEdZ_decoupID'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['dThetaEdZ_decoupID'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['dThetaEdZ_decoupID'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['dThetaEdZ_decoupID'][:] = np.nan           ## fill with nans to account for missing files when populating
    obs['sondes']['thetaE_decoupleZ'] = np.zeros([np.size(obs['sondes']['dThetaE'],0)])
    data1['thetaE_decoupleZ'] = np.zeros([np.size(data1['dThetaE'],0)])
    data2['thetaE_decoupleZ'] = np.zeros([np.size(data2['dThetaE'],0)])
    data3['thetaE_decoupleZ'] = np.zeros([np.size(data3['dThetaE'],0)])
    data3['thetaE_decoupleZ'][:] = np.nan           ## fill with nans to account for missing files when populating

    #### find maximum dThetaE
    for i in range(0, np.size(obs['sondes']['doy_drift'])):

        #### ---------------------------------------------------------------
        #### find strongest inversion <3000m
        #### ---------------------------------------------------------------
        # obs['sondes']['thetaE_invbaseID'][i] = np.where(np.squeeze(obs['sondes']['dThetaE'][i,1:27]) ==
        #         np.squeeze(np.nanmax(obs['sondes']['dThetaE'][i,1:27])))[0][0]
        # data1['thetaE_invbaseID'][i] = np.where(np.squeeze(data1['dThetaE'][i,1:27]) ==
        #         np.squeeze(np.nanmax(data1['dThetaE'][i,1:27])))[0][0]
        # data2['thetaE_invbaseID'][i] = np.where(np.squeeze(data2['dThetaE'][i,1:27]) ==
        #         np.squeeze(np.nanmax(data2['dThetaE'][i,1:27])))[0][0]
        # if np.nanmax(data3['dThetaE'][i,lt3000]) >= 0.0:       ### ignore missing files (filled with nans)
        #     data3['thetaE_invbaseID'][i] = np.where(np.squeeze(data3['dThetaE'][i,1:27]) ==
        #         np.squeeze(np.nanmax(data3['dThetaE'][i,1:27])))[0][0]

        #### ---------------------------------------------------------------
        #### find maximum local gradient <3km
        #### ---------------------------------------------------------------
        obs['sondes']['dThetaEdZ_invbaseID'][i] = np.where(obs['sondes']['dThetaEdZ'][i,:] ==
            np.nanmax(obs['sondes']['dThetaEdZ'][i,:]))[0][0]
        data1['dThetaEdZ_invbaseID'][i] = np.where(data1['dThetaEdZ'][i,:] ==
                np.nanmax(data1['dThetaEdZ'][i,:]))[0][0]
        data2['dThetaEdZ_invbaseID'][i] = np.where(data2['dThetaEdZ'][i,:] ==
                np.nanmax(data2['dThetaEdZ'][i,:]))[0][0]
        if np.nanmax(data3['dThetaE'][i,lt3000]) >= 0.0:       ### ignore missing files (filled with nans)
            data3['dThetaEdZ_invbaseID'][i] = np.where(data3['dThetaEdZ'][i,:] ==
                    np.nanmax(data3['dThetaEdZ'][i,:]))[0][0]

        #### ---------------------------------------------------------------
        #### find second largest gradient
        #### ---------------------------------------------------------------
        obs['sondes']['dThetaEdZ_2ndinvID'][i] = np.where(obs['sondes']['dThetaEdZ'][i,:] ==
            np.sort(obs['sondes']['dThetaEdZ'][i,:])[::-1][1])[0][0]
        data1['dThetaEdZ_2ndinvID'][i] = np.where(data1['dThetaEdZ'][i,:] ==
            np.sort(data1['dThetaEdZ'][i,:])[::-1][1])[0][0]
        data2['dThetaEdZ_2ndinvID'][i] = np.where(data2['dThetaEdZ'][i,:] ==
            np.sort(data2['dThetaEdZ'][i,:])[::-1][1])[0][0]
        if np.nanmax(data3['dThetaE'][i,lt3000]) >= 0.0:       ### ignore missing files (filled with nans)
            data3['dThetaEdZ_2ndinvID'][i] = np.where(data3['dThetaEdZ'][i,:] ==
                np.sort(data3['dThetaEdZ'][i,:])[::-1][1])[0][0]

        #### ---------------------------------------------------------------
        ### if second largest dThetaEdZ is at invbaseID - 1, then reassign invbaseID -> 2ndinvID
        ###         in this case, 2ndinvID set to 3rd largest gradient
        #### ---------------------------------------------------------------
        if obs['sondes']['dThetaEdZ_invbaseID'][i]-1 == obs['sondes']['dThetaEdZ_2ndinvID'][i]:
                obs['sondes']['dThetaEdZ_invbaseID'][i] = obs['sondes']['dThetaEdZ_2ndinvID'][i]
                obs['sondes']['dThetaEdZ_2ndinvID'][i] = np.where(obs['sondes']['dThetaEdZ'][i,:] ==
                    np.sort(obs['sondes']['dThetaEdZ'][i,:])[::-1][2])[0][0]
        if data1['dThetaEdZ_invbaseID'][i]-1 == data1['dThetaEdZ_2ndinvID'][i]:
                data1['dThetaEdZ_invbaseID'][i] = data1['dThetaEdZ_2ndinvID'][i]
                data1['dThetaEdZ_2ndinvID'][i] = np.where(data1['dThetaEdZ'][i,:] ==
                    np.sort(data1['dThetaEdZ'][i,:])[::-1][2])[0][0]
        if data2['dThetaEdZ_invbaseID'][i]-1 == data2['dThetaEdZ_2ndinvID'][i]:
                data2['dThetaEdZ_invbaseID'][i] = data2['dThetaEdZ_2ndinvID'][i]
                data2['dThetaEdZ_2ndinvID'][i] = np.where(data2['dThetaEdZ'][i,:] ==
                    np.sort(data2['dThetaEdZ'][i,:])[::-1][2])[0][0]
        if data3['dThetaEdZ_invbaseID'][i]-1 == data3['dThetaEdZ_2ndinvID'][i]:
                data3['dThetaEdZ_invbaseID'][i] = data3['dThetaEdZ_2ndinvID'][i]
                data3['dThetaEdZ_2ndinvID'][i] = np.where(data3['dThetaEdZ'][i,:] ==
                    np.sort(data3['dThetaEdZ'][i,:])[::-1][2])[0][0]

        #### --------------------------------------------------------------------------------------------------------
        #### ---------------------------------------------------------------
        #### if 2nd inv below 1000m, set decoupID to 2ndinvID
        #### ---------------------------------------------------------------
        #### --------------------------------------------------------------------------------------------------------
        hgt_index = 1000.0
        hind = np.where(data1['universal_height'] <= hgt_index)
        if np.logical_and(obs['sondes']['dThetaEdZ_2ndinvID'][i] < hind[0][-1], obs['sondes']['dThetaEdZ_2ndinvID'][i] < obs['sondes']['dThetaEdZ_invbaseID'][i]):
            obs['sondes']['dThetaEdZ_decoupID'][i] = obs['sondes']['dThetaEdZ_2ndinvID'][i]
            obs['sondes']['dThetaEdZ_2ndinvID'][i] = np.where(obs['sondes']['dThetaEdZ'][i,:] ==
                np.sort(obs['sondes']['dThetaEdZ'][i,:])[::-1][3])[0][0]
        if np.logical_and(data1['dThetaEdZ_2ndinvID'][i] < hind[0][-1], data1['dThetaEdZ_2ndinvID'][i] < data1['dThetaEdZ_invbaseID'][i]):
            data1['dThetaEdZ_decoupID'][i] = data1['dThetaEdZ_2ndinvID'][i]
            data1['dThetaEdZ_2ndinvID'][i] = np.where(data1['dThetaEdZ'][i,:] ==
                np.sort(data1['dThetaEdZ'][i,:])[::-1][3])[0][0]
        if np.logical_and(data2['dThetaEdZ_2ndinvID'][i] < hind[0][-1], data2['dThetaEdZ_2ndinvID'][i] < data2['dThetaEdZ_invbaseID'][i]):
            data2['dThetaEdZ_decoupID'][i] = data2['dThetaEdZ_2ndinvID'][i]
            data2['dThetaEdZ_2ndinvID'][i] = np.where(data2['dThetaEdZ'][i,:] ==
                np.sort(data2['dThetaEdZ'][i,:])[::-1][3])[0][0]
        if np.logical_and(data3['dThetaEdZ_2ndinvID'][i] < hind[0][-1], data3['dThetaEdZ_2ndinvID'][i] < data3['dThetaEdZ_invbaseID'][i]):
            data3['dThetaEdZ_decoupID'][i] = data3['dThetaEdZ_2ndinvID'][i]
            data3['dThetaEdZ_2ndinvID'][i] = np.where(data3['dThetaEdZ'][i,:] ==
                np.sort(data3['dThetaEdZ'][i,:])[::-1][3])[0][0]

        #### ---------------------------------------------------------------
        ### if decoupID is at invbaseID - 1, then reassign invbaseID -> decoupID
        ###         in this case, decoupID -> 0.0
        #### ---------------------------------------------------------------
        if obs['sondes']['dThetaEdZ_invbaseID'][i]-1 == obs['sondes']['dThetaEdZ_decoupID'][i]:
                obs['sondes']['dThetaEdZ_invbaseID'][i] = obs['sondes']['dThetaEdZ_decoupID'][i]
                obs['sondes']['dThetaEdZ_decoupID'][i] = 0.0
        if data1['dThetaEdZ_invbaseID'][i]-1 == data1['dThetaEdZ_decoupID'][i]:
                data1['dThetaEdZ_invbaseID'][i] = data1['dThetaEdZ_decoupID'][i]
                data1['dThetaEdZ_decoupID'][i] = 0.0
        if data2['dThetaEdZ_invbaseID'][i]-1 == data2['dThetaEdZ_decoupID'][i]:
                data2['dThetaEdZ_invbaseID'][i] = data2['dThetaEdZ_decoupID'][i]
                data2['dThetaEdZ_decoupID'][i] = 0.0
        if data3['dThetaEdZ_invbaseID'][i]-1 == data3['dThetaEdZ_decoupID'][i]:
                data3['dThetaEdZ_invbaseID'][i] = data3['dThetaEdZ_decoupID'][i]
                data3['dThetaEdZ_decoupID'][i] = 0.0

        #### --------------------------------------------------------------------------------------------------------
        #### ---------------------------------------------------------------
        #### search for stable surface layers
        #### ---------------------------------------------------------------
        #### --------------------------------------------------------------------------------------------------------
        ### if invbase = 0, set to 2nd largest gradient instead
        if obs['sondes']['dThetaEdZ_invbaseID'][i] == 0.0:
            obs['sondes']['dThetaEdZ_invbaseID'][i] = obs['sondes']['dThetaEdZ_2ndinvID'][i]
            obs['sondes']['dThetaEdZ_decoupID'][i] = -1.0
        if data1['dThetaEdZ_invbaseID'][i] == 0.0:
            data1['dThetaEdZ_invbaseID'][i] = data1['dThetaEdZ_2ndinvID'][i]
            data1['dThetaEdZ_decoupID'][i] = -1.0
        if data2['dThetaEdZ_invbaseID'][i] == 0.0:
            data2['dThetaEdZ_invbaseID'][i] = data2['dThetaEdZ_2ndinvID'][i]
            data2['dThetaEdZ_decoupID'][i] = -1.0
        if data3['dThetaEdZ_invbaseID'][i] == 0.0:
            data3['dThetaEdZ_invbaseID'][i] = data3['dThetaEdZ_2ndinvID'][i]
            data3['dThetaEdZ_decoupID'][i] = -1.0

        ### if invbase = 1 or 2, means previous criteria has shown this is a stable surface layer
        if np.logical_or(obs['sondes']['dThetaEdZ_invbaseID'][i] == 1.0, obs['sondes']['dThetaEdZ_invbaseID'][i] == 2.0):
            ###         look first for largest difference, if !=0, then set as invbase
            ###                 only look below 2km
            obs['sondes']['dThetaEdZ_invbaseID'][i] = np.where(obs['sondes']['dThetaE'][i,:] == np.nanmax(obs['sondes']['dThetaE'][i,:21]))[0][0]
            ###         look instead for 2nd largest difference in ThetaE (largest will likely also be invbase=0)
            if np.logical_or(obs['sondes']['dThetaEdZ_invbaseID'][i] == 0.0, obs['sondes']['dThetaEdZ_invbaseID'][i] == 1.0):
                obs['sondes']['dThetaEdZ_invbaseID'][i] = np.where(obs['sondes']['dThetaE'][i,:] == np.sort(obs['sondes']['dThetaE'][i,:21])[::-1][1])[0][0]
        if np.logical_or(data1['dThetaEdZ_invbaseID'][i] == 1.0, data1['dThetaEdZ_invbaseID'][i] == 2.0):
            data1['dThetaEdZ_invbaseID'][i] = np.where(data1['dThetaE'][i,:] == np.nanmax(data1['dThetaE'][i,:21]))[0][0]
            if np.logical_or(data1['dThetaEdZ_invbaseID'][i] == 0.0, data1['dThetaEdZ_invbaseID'][i] == 1.0):
                data1['dThetaEdZ_invbaseID'][i] = np.where(data1['dThetaE'][i,:] == np.sort(data1['dThetaE'][i,:21])[::-1][1])[0][0]
        if np.logical_or(data2['dThetaEdZ_invbaseID'][i] == 1.0, data2['dThetaEdZ_invbaseID'][i] == 2.0):
            data2['dThetaEdZ_invbaseID'][i] = np.where(data2['dThetaE'][i,:] == np.nanmax(data2['dThetaE'][i,:21]))[0][0]
            if np.logical_or(data2['dThetaEdZ_invbaseID'][i] == 0.0, data2['dThetaEdZ_invbaseID'][i] == 1.0):
                data2['dThetaEdZ_invbaseID'][i] = np.where(data2['dThetaE'][i,:] == np.sort(data2['dThetaE'][i,:21])[::-1][1])[0][0]
        if np.logical_or(data3['dThetaEdZ_invbaseID'][i] == 1.0, data3['dThetaEdZ_invbaseID'][i] == 2.0):
            data3['dThetaEdZ_invbaseID'][i] = np.where(data3['dThetaE'][i,:] == np.nanmax(data3['dThetaE'][i,:21]))[0][0]
            if np.logical_or(data3['dThetaEdZ_invbaseID'][i] == 0.0, data3['dThetaEdZ_invbaseID'][i] == 1.0):
                data3['dThetaEdZ_invbaseID'][i] = np.where(data3['dThetaE'][i,:] == np.sort(data3['dThetaE'][i,:21])[::-1][1])[0][0]

        #### ---------------------------------------------------------------
        #### save timeseries of invbase and decoupleZ heights
        #### ---------------------------------------------------------------
        obs['sondes']['thetaE_invbase'][i] = data1['universal_height'][int(obs['sondes']['dThetaEdZ_invbaseID'][i])] ##data1['universal_height'][int(obs['sondes']['thetaE_invbaseID'][i])]
        data1['thetaE_invbase'][i] = data1['universal_height'][int(data1['dThetaEdZ_invbaseID'][i])] ## data1['universal_height'][int(data1['thetaE_invbaseID'][i])]
        data2['thetaE_invbase'][i] = data1['universal_height'][int(data2['dThetaEdZ_invbaseID'][i])] ## data1['universal_height'][int(data2['thetaE_invbaseID'][i])]
        if data3['dThetaEdZ_invbaseID'][i] >= 0.0: data3['thetaE_invbase'][i] = data1['universal_height'][int(data3['dThetaEdZ_invbaseID'][i])]
                ##data3['thetaE_invbase'][i] = data1['universal_height'][int(data3['thetaE_invbaseID'][i])]

        obs['sondes']['thetaE_decoupleZ'][i] = data1['universal_height'][int(obs['sondes']['dThetaEdZ_decoupID'][i])]
        data1['thetaE_decoupleZ'][i] = data1['universal_height'][int(data1['dThetaEdZ_decoupID'][i])]
        data2['thetaE_decoupleZ'][i] = data1['universal_height'][int(data2['dThetaEdZ_decoupID'][i])]
        if data3['dThetaEdZ_decoupID'][i] >= 0.0: data3['thetaE_decoupleZ'][i] = data1['universal_height'][int(data3['dThetaEdZ_decoupID'][i])]

    #### ---------------------------------------------------------------
    #### load working cloudnet data to mark cloud depth
    #### ---------------------------------------------------------------
    obs_data = np.load('../Cloudnet/obs_Cv_data.npy').item()
    um_data = np.load('../Cloudnet/um_Cv_data.npy').item()
    misc_data = np.load('../Cloudnet/casim_Cv_data.npy').item()
    ifs_data = np.load('../Cloudnet/ifs_Cv_data.npy').item()

    obs_cv = obs_data['Cv'][::6,:]
    um_ra2m_cv = um_data['model_Cv_filtered'][::6,:].data
    um_casim_cv = misc_data['model_Cv_filtered'][::6,:].data
    ecmwf_ifs_cv = ifs_data['model_snow_Cv_filtered'][::6,:].data

    #### ---------------------------------------------------------------
    #### Look at drift sondes only from Jutta's inversion code
    #### ---------------------------------------------------------------
    drift = np.where(np.logical_and(obs['inversions']['doy'] >= 225.9, obs['inversions']['doy'] <= 258.0))
    ### save in dict for ease
    obs['inversions']['doy_drift'] = obs['inversions']['doy'][drift]

    #### ---------------------------------------------------------------
    #### include UM BL type information
    #### ---------------------------------------------------------------
    #### list of BL types from UM documentation
    doc = ['(1) Stable BL', '(2) Sc over stable SL', '(3) Well-mixed BL', '(4) Unstable BL, dSc not o/Cu',
        '(5) dSc o/Cu', '(6) Cu-capped BL', '(7) Shear-dom unstable BL']

    bltype1 = data1['bl_type'][data1['hrly_flag']][::6].data
    bltype2 = data2['bl_type'][data2['hrly_flag']][::6].data

    # Type I: Stable boundary layer (with or without cloud)
    # Type II: Boundary layer with stratocumulus over a stable near-surface layer
    # Type III: Well mixed boundary layer
    # Type IV: Unstable boundary layer with a decoupled Sc layer not over cumulus
    # Type V: Boundary layer with a decoupled Sc layer layer over cumulus
    # Type VI: Cumulus-capped boundary layer
    # Type VII: Shear-dominated unstable layer

    ##################################################
    ##################################################
    #### create figure and axes instances
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=SMALL_SIZE)

    #### ---------------------------------------------------------------
    #### save quicklooks for reference
    #### ---------------------------------------------------------------
    i = 0
    for i in range(0, np.size(obs['sondes']['doy_drift'])):
        fig = plt.figure(figsize=(11,6))
        ax  = fig.add_axes([0.1,0.1,0.5,0.8])   # left, bottom, width, height

        ########### RADIOSONDES
        plt.plot(obs['sondes']['thetaE_driftSondes_UM'][i,:],data1['universal_height'], color = 'k',)# label = 'sonde-interpd')
        plt.plot(np.squeeze(obs['sondes']['thetaE_driftSondes_UM'][i,int(obs['sondes']['dThetaEdZ_invbaseID'][i])]),
            np.squeeze(data1['universal_height'][int(obs['sondes']['dThetaEdZ_invbaseID'][i])]),
            's', markersize = 8, color = 'darkgrey', markeredgecolor = 'k', label = 'sonde-interpd max d$\Theta_{E}$/dZ')
        # plt.plot(np.squeeze(obs['sondes']['thetaE_driftSondes_UM'][i,int(obs['sondes']['thetaE_invbaseID'][i])]),
        #     np.squeeze(data1['universal_height'][int(obs['sondes']['thetaE_invbaseID'][i])]),
        #     's', markersize = 8, color = 'k', label = 'sonde-interpd max d$\Theta_{E}$')
        # if obs['sondes']['dThetaEdZ_2ndinvID'][i] > 0.0:
        #     plt.plot(np.squeeze(obs['sondes']['thetaE_driftSondes_UM'][i,int(obs['sondes']['dThetaEdZ_2ndinvID'][i])]),
        #         np.squeeze(data1['universal_height'][int(obs['sondes']['dThetaEdZ_2ndinvID'][i])]),
        #         '+', color = 'k', label = 'sonde-interpd 2nd max d$\Theta_{E}$')
        # if obs['sondes']['dThetaEdZ_decoupID'][i] > 0.0:
        #     plt.plot(np.squeeze(obs['sondes']['thetaE_driftSondes_UM'][i,int(obs['sondes']['dThetaEdZ_decoupID'][i])]),
        #         np.squeeze(data1['universal_height'][int(obs['sondes']['dThetaEdZ_decoupID'][i])]),
        #         'o', color = 'darkgrey', markeredgecolor = 'k', label = 'sonde-interpd decoup surf d$\Theta_{E}$/dz')

        ############# ECMWF_IFS
        plt.plot(data3['thetaE_6hrly_UM'][i,:],data1['universal_height'], color = 'gold')#, label = 'ifs-interpd')
        if data3['dThetaEdZ_invbaseID'][i] > 0.0:     ### ignore nans (missing files)
            plt.plot(np.squeeze(data3['thetaE_6hrly_UM'][i,int(data3['dThetaEdZ_invbaseID'][i])]),
                np.squeeze(data1['universal_height'][int(data3['dThetaEdZ_invbaseID'][i])]),
                's', markersize = 8, color = 'gold', markeredgecolor = 'saddlebrown', label = 'ifs-interpd max d$\Theta_{E}$/dZ')
        # if data3['thetaE_invbaseID'][i] > 0.0:     ### ignore nans (missing files)
        #     plt.plot(np.squeeze(data3['thetaE_6hrly_UM'][i,int(data3['thetaE_invbaseID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data3['thetaE_invbaseID'][i])]),
        #         's', markersize = 8, color = 'gold', label = 'ifs-interpd max d$\Theta_{E}$')
        # if data3['dThetaEdZ_2ndinvID'][i] > 0.0:     ### ignore nans (missing files)
        #     plt.plot(np.squeeze(data3['thetaE_6hrly_UM'][i,int(data3['dThetaEdZ_2ndinvID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data3['dThetaEdZ_2ndinvID'][i])]),
        #         '+', color = 'gold', label = 'ifs-interpd 2nd max d$\Theta_{E}$')
        # if data3['dThetaEdZ_decoupID'][i] > 0.0:     ### ignore nans (missing files)
        #     plt.plot(np.squeeze(data3['thetaE_6hrly_UM'][i,int(data3['dThetaEdZ_decoupID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data3['dThetaEdZ_decoupID'][i])]),
        #         'o', color = 'gold', markeredgecolor = 'saddlebrown', label = 'ifs-interpd decoup surf d$\Theta_{E}$/dz')

        ############## UM_RA2M
        plt.plot(data1['thetaE_6hrly'][i,data1['universal_height_UMindex']], data1['universal_height'], color = 'darkblue')#, label = 'um_ra2m')
        plt.plot(np.squeeze(data1['thetaE_6hrly_UM'][i,int(data1['dThetaEdZ_invbaseID'][i])]),
            np.squeeze(data1['universal_height'][int(data1['dThetaEdZ_invbaseID'][i])]),
            's', markersize = 8, color = 'darkblue', markeredgecolor = 'midnightblue', label = 'um_ra2m max d$\Theta_{E}$/dZ')
        # plt.plot(np.squeeze(data1['thetaE_6hrly_UM'][i,int(data1['thetaE_invbaseID'][i])]),
        #     np.squeeze(data1['universal_height'][int(data1['thetaE_invbaseID'][i])]),
        #     's', markersize = 8, color = 'darkblue', label = 'um_ra2m max d$\Theta_{E}$')
        # if data1['dThetaEdZ_2ndinvID'][i] > 0.0:
        #     plt.plot(np.squeeze(data1['thetaE_6hrly_UM'][i,int(data1['dThetaEdZ_2ndinvID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data1['dThetaEdZ_2ndinvID'][i])]),
        #         '+', color = 'darkblue', label = 'um_ra2m 2nd max d$\Theta_{E}$')
        # if data1['dThetaEdZ_decoupID'][i] > 0.0:
        #     plt.plot(np.squeeze(data1['thetaE_6hrly_UM'][i,int(data1['dThetaEdZ_decoupID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data1['dThetaEdZ_decoupID'][i])]),
        #         'o', color = 'darkblue', markeredgecolor = 'midnightblue', label = 'um_ra2m decoup surf d$\Theta_{E}$/dz')

        ############## UM_CASIM-100
        plt.plot(data2['thetaE_6hrly'][i,data1['universal_height_UMindex']], data1['universal_height'], color = 'mediumseagreen')#, label = 'um_casim-100')
        plt.plot(np.squeeze(data2['thetaE_6hrly_UM'][i,int(data2['dThetaEdZ_invbaseID'][i])]),
            np.squeeze(data1['universal_height'][int(data2['dThetaEdZ_invbaseID'][i])]),
            's', markersize = 8, color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = 'um_casim-100 max d$\Theta_{E}$/dZ')
        # plt.plot(np.squeeze(data2['thetaE_6hrly_UM'][i,int(data2['thetaE_invbaseID'][i])]),
        #     np.squeeze(data1['universal_height'][int(data2['thetaE_invbaseID'][i])]),
        #     's', markersize = 8, color = 'mediumseagreen', label = 'um_casim-100 max d$\Theta_{E}$')
        # if data2['dThetaEdZ_2ndinvID'][i] > 0.0:
        #     plt.plot(np.squeeze(data2['thetaE_6hrly_UM'][i,int(data2['dThetaEdZ_2ndinvID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data2['dThetaEdZ_2ndinvID'][i])]),
        #         '+', color = 'mediumseagreen', label = 'um_casim-100 2nd max d$\Theta_{E}$')
        # if data2['dThetaEdZ_decoupID'][i] > 0.0:
        #     plt.plot(np.squeeze(data2['thetaE_6hrly_UM'][i,int(data2['dThetaEdZ_decoupID'][i])]),
        #         np.squeeze(data1['universal_height'][int(data2['dThetaEdZ_decoupID'][i])]),
        #         'o', color = 'mediumseagreen', markeredgecolor = 'darkslategrey', label = 'um_casim-100 decoup surf d$\Theta_{E}$/dz')

        ########## plot model diagnosed boundary layer depth as marker
        plt.plot([298,308],[np.squeeze(obs['inversions']['invbase'][drift][i]),np.squeeze(obs['inversions']['invbase'][drift][i])],
            'x--', markersize = 7, color = 'k', linewidth = 1, label = 'radiosondes (Vuellers et al., 2020)')
        plt.plot([298,308],[bldepth1[::6][i],bldepth1[::6][i]],
            'x--', markersize = 7, color = 'darkblue', linewidth = 1, label = 'um_ra2m bl_depth')
        plt.plot([298,308],[bldepth2[::6][i],bldepth2[::6][i]],
            'x--', markersize = 7, color = 'mediumseagreen', linewidth = 1, label = 'um_casim-100 bl_depth')
        plt.plot([298,308],[bldepth3[::6][i],bldepth3[::6][i]],
            'x--', markersize = 7, color = 'gold', linewidth = 1, label = 'ecmwf_ifs sfc_bl_height')

        ##### plot cloudnet cloud fraction as markers
        #####       need to find the correct cloudnet index for this model timestep
        timeind = 0.0       ### initialise
        timeind = np.where(data1['time_6hrly'][i] == um_data['time'][::6])
        if np.size(timeind) > 0:        ### if there is cloudnet data at this sonde timestep
            # tmp1 = np.zeros(len(np.where(um_ra2m_cv[timeind[0][0],:29] > 0)[0]))
            tmp1 = np.zeros(len(um_ra2m_cv[timeind[0][0],:29]))
            tmp1[:] = 310
            # plt.plot(tmp1,um_data['height'][i,np.where(um_ra2m_cv[timeind[0][0],:29] > 0)[0]],
            #     '.', markersize = 6, color = 'darkblue', label = 'um_ra2m Cv > 0 at t = ' + str(um_data['time'][::6][timeind[0][0]]))
            ax.scatter(tmp1,um_data['height'][i,:29], s = 7, c = um_ra2m_cv[timeind[0][0],:29],
                vmin = 0, vmax = 1, cmap = mpl_cm.Blues, label = 'Cv at t = ' + str(um_data['time'][::6][timeind[0][0]]))
            # tmp2 = np.zeros(len(np.where(um_casim_cv[timeind[0][0],:29] > 0)[0]))
            tmp2 = np.zeros(len(um_casim_cv[timeind[0][0],:29]))
            tmp2[:] = 312
            # plt.plot(tmp2,misc_data['height'][i,np.where(um_casim_cv[timeind[0][0],:29] > 0)[0]],
            #     '.', markersize = 6, color = 'mediumseagreen', label = 'um_casim-100 Cv > 0')
            ax.scatter(tmp2,misc_data['height'][i,:29], s = 7, c = um_casim_cv[timeind[0][0],:29],
                vmin = 0, vmax = 1, cmap = mpl_cm.Greens)#, label = 'um_casim-100 Cv > 0')
            # tmp3 = np.zeros(len(np.where(ecmwf_ifs_cv[timeind[0][0],:29] > 0)[0]))
            tmp3 = np.zeros(len(ecmwf_ifs_cv[timeind[0][0],:29]))
            tmp3[:] = 314
            # plt.plot(tmp3,ifs_data['height'][i,np.where(ecmwf_ifs_cv[timeind[0][0],:29] > 0)[0]],
            #     '.', markersize = 6, color = 'gold', label = 'ecmwf_ifs Cv > 0')
            ax.scatter(tmp3,ifs_data['height'][i,:29], s = 7, c = ecmwf_ifs_cv[timeind[0][0],:29],
                vmin = 0, vmax = 1, cmap = mpl_cm.Oranges)#, label = 'ecmwf_ifs Cv > 0')
            # tmp0 = np.zeros(len(np.where(obs_cv[timeind[0][0],:29] > 0)[0]))
            tmp0 = np.zeros(len(obs_cv[timeind[0][0],:29]))
            tmp0[:] = 316
            # plt.plot(tmp0,obs_data['height'][i,np.where(obs_cv[timeind[0][0],:29] > 0)[0]],
            #     '.', markersize = 6, color = 'k', label = 'obs Cv > 0')
            ax.scatter(tmp0,obs_data['height'][i,:29], s = 7, c = obs_cv[timeind[0][0],:29],
                vmin = 0, vmax = 1, cmap = mpl_cm.Greys)#, label = 'obs Cv > 0')

        plt.title('Inversion identification test DOY ' + str(np.round(obs['sondes']['doy_drift'][i],2)) + '\n Cloudnet cloud fractions indicated on RHS')
        plt.xlabel('$\Theta_{E}$ [K]')
        plt.ylabel('Z [m]')
        plt.ylim([0,3100])
        plt.xlim([260,320])
        plt.grid(axis = 'y')

        # plt.annotate('UM Boundary layer classifications: \n' +
        #     'UM_RA2M: ' + doc[int(data1['bl_type'][i])-1] + '\n' +
        #     'UM_CASIM-100: ' + doc[int(data2['bl_type'][i])-1],
        #     xy=(362,2550),xytext=(362,2550),
        #     fontsize=12)

        # print(i)
        ax.text(323, 2550, 'UM Boundary layer classifications: \n' +
            'UM_RA2M: ' + doc[int(bltype1[i])-1] + '\n' +
            'UM_CASIM-100: ' + doc[int(bltype2[i])-1],
            fontsize = 12)

        legend2 = ax.legend(bbox_to_anchor=(0.7, 0.0, 1., .102), loc=4, ncol=1)
        plt.savefig('../FIGS/inversionIdent/InvIdent_ThetaE_doy' + str(np.round(obs['sondes']['doy_drift'][i],1)) + '.png')
        if i == 0:
            plt.show()
        else:
            plt.close()

    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(13,7))
    ax  = fig.add_axes([0.1,0.15,0.55,0.7])   # left, bottom, width, height
    plt.plot(np.squeeze(obs['inversions']['doy_drift']),np.squeeze(obs['inversions']['invbase'][drift]),
        color = 'k', label = 'JV_Radiosondes')
    plt.plot(obs['sondes']['doy_drift'], obs['sondes']['thetaE_invbase'],
        'ks', label = 'radiosondes')
    plt.plot(data1['time_6hrly'], data1['thetaE_invbase'],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = 'um_ra2m')
    plt.plot(data2['time_6hrly'], data2['thetaE_invbase'],
        'v', color = 'mediumseagreen',  markeredgecolor = 'darkslategrey', label = 'um_casim-100')
    plt.plot(data3['time_6hrly'], data3['thetaE_invbase'],
        'd', color = 'gold',  markeredgecolor = 'saddlebrown', label = 'ecmwf_ifs')
    plt.xlabel('DOY')
    plt.ylabel('Inversion base height [m]')
    plt.legend()

    ax  = fig.add_axes([0.74,0.3,0.22,0.38])   # left, bottom, width, height
    plt.plot([0, 3000],[0, 3000], '--', color = 'lightgrey')
    plt.plot(obs['sondes']['thetaE_invbase'],data1['thetaE_invbase'],
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = 'um_ra2m')
    plt.plot(obs['sondes']['thetaE_invbase'], data2['thetaE_invbase'],
        'v', color = 'mediumseagreen',  markeredgecolor = 'darkslategrey', label = 'um_casim-100')
    plt.plot(obs['sondes']['thetaE_invbase'], data3['thetaE_invbase'],
        'd', color = 'gold',  markeredgecolor = 'saddlebrown', label = 'ecmwf_ifs')
    plt.xlabel('Obs')
    plt.ylabel('Models')
    plt.ylim([0, 3000])
    plt.xlim([0, 3000])
    plt.savefig('../FIGS/inversionIdent/InvbaseTS_dThetaEdZ.png')
    plt.show()


    ###---------------------------------------------------------------------------------------------
    ###         Save inversions to numpy dictionary
    ###---------------------------------------------------------------------------------------------

    data_obs = {}
    data_obs['invbase'] = obs['sondes']['thetaE_invbase'][:]
    data_obs['time'] = obs['sondes']['doy_drift'][:]
    data_obs['sfmlheight'] = obs['sondes']['thetaE_decoupleZ'][:]
    data_obs['height'] = data1['universal_height']
    np.save('obs_inversions_v3', data_obs)

    data_1 = {}
    data_1['invbase'] = data1['thetaE_invbase'][:]
    data_1['time'] = data1['time_6hrly'][:]
    data_1['sfmlheight'] = data1['thetaE_decoupleZ'][:]
    data_1['height'] = data1['universal_height']
    np.save('um_ra2m_inversions_v3', data_1)

    data_2 = {}
    data_2['invbase'] = data2['thetaE_invbase'][:]
    data_2['time'] = data2['time_6hrly'][:]
    data_2['sfmlheight'] = data2['thetaE_decoupleZ'][:]
    data_2['height'] = data1['universal_height']
    np.save('um_casim-100_inversions_v3', data_2)

    data_3 = {}
    data_3['invbase'] = data3['thetaE_invbase'][:]
    data_3['time'] = data3['time_6hrly'][:]
    data_3['sfmlheight'] = data3['thetaE_decoupleZ'][:]
    data_3['height'] = data1['universal_height']
    np.save('ecmwf_ifs_inversions_v3', data_3)

    return data1, data2, data3, obs

def radarRefl_Sandeep(data1, data2, data3, data4, obs, doy, label1, label2, label3, label4):

    '''
    Plots model radar reflectivity timeseries for selected MOCCHA time periods
    '''

    ### dates as strings
    d2 = '20180911'
    d3 = '20180913'
    d4 = '20180831'
    d5 = '20180823'
    d6 = '20180903'
    d7 = '20180905'

    ### calculate DOYs
    doy2 = calcTime_Date2DOY(d2)
    doy3 = calcTime_Date2DOY(d3)
    doy4 = calcTime_Date2DOY(d4)
    doy5 = calcTime_Date2DOY(d5)
    doy6 = calcTime_Date2DOY(d6)
    doy7 = calcTime_Date2DOY(d7)

    ### hour intervals for each period
    h2 = np.arange(0,15)
    h3 = np.arange(0,15)
    h4 = np.arange(0,5)
    h5 = np.arange(12,25)
    h6 = np.arange(19,23)
    h7 = np.arange(5,16)

    ### define periods as time of day (doy)
    p2 = doy2 + h2/24.0
    p3 = doy3 + h3/24.0
    p4 = doy4 + h4/24.0
    p5 = doy5 + h5/24.0
    p6 = doy6 + h6/24.0
    p7 = doy7 + h7/24.0

    ### find model data for each period
    i2 = np.where(np.logical_and(data1['time_hrly'] >= p2[0], data1['time_hrly'] <= p2[-1]))
    i3 = np.where(np.logical_and(data1['time_hrly'] >= p3[0], data1['time_hrly'] <= p3[-1]))
    i4 = np.where(np.logical_and(data1['time_hrly'] >= p4[0], data1['time_hrly'] <= p4[-1]))
    i5 = np.where(np.logical_and(data1['time_hrly'] >= p5[0], data1['time_hrly'] <= p5[-1]))
    i6 = np.where(np.logical_and(data1['time_hrly'] >= p6[0], data1['time_hrly'] <= p6[-1]))
    i7 = np.where(np.logical_and(data1['time_hrly'] >= p7[0], data1['time_hrly'] <= p7[-1]))

    #### remove flagged data
    data1['radr_refl'][data1['radr_refl'] == -9999.0] = np.nan
    data2['radr_refl'][data2['radr_refl'] == -9999.0] = np.nan
    data4['radr_refl'][data4['radr_refl'] == -9999.0] = np.nan

    data1['cloud_fraction'][data1['cloud_fraction'] == -9999.0] = np.nan
    data2['cloud_fraction'][data2['cloud_fraction'] == -9999.0] = np.nan
    data4['cloud_fraction'][data4['cloud_fraction'] == -9999.0] = np.nan


    ####### ----------------------------------------------------------------
    #######     FIGURE
    ####### ----------------------------------------------------------------
    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.figure(figsize=(14,10))
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.2)

    ## color limits
    cmax = 10.0
    cmin = -35.0
    h = 5e3

    # plt.subplot(3, 2, 1)
    # plt.pcolormesh(data2['time_hrly'][i2], data2['height'], np.transpose(data2['radr_refl'][i2[0],:]),
    #     vmin = cmin, vmax = cmax)
    # # plt.xlim([p2[0], p2[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P2) ' + d2 + ': ' + str(h2[0]) + ' - ' + str(h2[-1]) + 'h')
    # ax = plt.gca()
    # ax.set_xticks(p2)
    # ax.set_xticklabels(h2)
    #
    # plt.subplot(3, 2, 2)
    # ax = plt.gca()
    # plt.pcolormesh(data2['time_hrly'][i3], data2['height'], np.transpose(data2['cloud_fraction'][i3[0],:]),
    #     #vmin = cmin, vmax = cmax
    #     )
    # # plt.xlim([p3[0], p3[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P3) ' + d3 + ': ' + str(h3[0]) + ' - ' + str(h3[-1]) + 'h')
    # ax.set_xticks(p3)
    # ax.set_xticklabels(h3)

    ax = plt.gca()
    plt.plot(data1['time_hrly'][i3], data1['bl_type'][i3], label = 'um_ra2m')
    plt.plot(data2['time_hrly'][i3], data2['bl_type'][i3], label = 'um_casim-100')
    plt.plot(data4['time_hrly'][i3], data4['bl_type'][i3], label = 'um_ra2t')
        #vmin = cmin, vmax = cmax
        # )
    # plt.xlim([p3[0], p3[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    plt.legend()
    plt.title('(P3) ' + d3 + ': ' + str(h3[0]) + ' - ' + str(h3[-1]) + 'h')
    ax.set_xticks(p3)
    ax.set_xticklabels(h3)
    #



    # plt.subplot(3, 2, 3)
    # ax = plt.gca()
    # plt.pcolormesh(data2['time_hrly'][i4], data2['height'], np.transpose(data2['radr_refl'][i4[0],:]),
    #     vmin = cmin, vmax = cmax)
    # # plt.xlim([p4[0], p4[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P4) ' + d4 + ': ' + str(h4[0]) + ' - ' + str(h4[-1]) + 'h')
    # ax.set_xticks(p4)
    # ax.set_xticklabels(h4)
    #
    # plt.subplot(3, 2, 4)
    # ax = plt.gca()
    # plt.pcolormesh(data2['time_hrly'][i5], data2['height'], np.transpose(data2['radr_refl'][i5[0],:]),
    #     vmin = cmin, vmax = cmax)
    # # plt.xlim([p5[0], p5[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P5) ' + d5 + ': ' + str(h5[0]) + ' - ' + str(h5[-1]) + 'h')
    # ax.set_xticks(p5)
    # ax.set_xticklabels(h5)
    #
    # plt.subplot(3, 2, 5)
    # ax = plt.gca()
    # plt.pcolormesh(data2['time_hrly'][i6], data2['height'], np.transpose(data2['radr_refl'][i6[0],:]),
    #     vmin = cmin, vmax = cmax)
    # # plt.xlim([p6[0], p6[-1]])
    # plt.colorbar()
    # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P6) ' + d6 + ': ' + str(h6[0]) + ' - ' + str(h6[-1]) + 'h')
    # ax.set_xticks(p6)
    # ax.set_xticklabels(h6)
    #
    # plt.subplot(3, 2, 6)
    # ax = plt.gca()
    # plt.pcolormesh(data2['time_hrly'][i7], data2['height'], np.transpose(data2['radr_refl'][i7[0],:]),
    #     vmin = cmin, vmax = cmax)
    # # plt.xlim([p7[0], p7[-1]])
    # plt.colorbar()
    # # plt.ylabel('Z [m]')
    # plt.ylim([0,h])
    # plt.title('(P7) ' + d7 + ': ' + str(h7[0]) + ' - ' + str(h7[-1]) + 'h')
    # ax.set_xticks(p7)
    # ax.set_xticklabels(h7)

    # plt.savefig('../FIGS/TKE/PeriodComparison_UM_CASIM-100_RadrRefl.png', dpi = 600)
    plt.show()

def TKEDissRate_Sandeep(data1, data2, data3, data4, obs, doy, label1, label2, label3, label4):

    '''
    Plots model radar reflectivity timeseries for selected MOCCHA time periods
    '''

    ### dates as strings
    d2 = '20180911'
    d3 = '20180913'
    d4 = '20180831'
    d5 = '20180823'
    d6 = '20180903'
    d7 = '20180905'

    ### calculate DOYs
    doy2 = calcTime_Date2DOY(d2)
    doy3 = calcTime_Date2DOY(d3)
    doy4 = calcTime_Date2DOY(d4)
    doy5 = calcTime_Date2DOY(d5)
    doy6 = calcTime_Date2DOY(d6)
    doy7 = calcTime_Date2DOY(d7)

    ### hour intervals for each period
    h2 = np.arange(0,15)
    h3 = np.arange(0,15)
    h4 = np.arange(0,5)
    h5 = np.arange(12,25)
    h6 = np.arange(19,23)
    h7 = np.arange(5,16)

    ### define periods as time of day (doy)
    p2 = doy2 + h2/24.0
    p3 = doy3 + h3/24.0
    p4 = doy4 + h4/24.0
    p5 = doy5 + h5/24.0
    p6 = doy6 + h6/24.0
    p7 = doy7 + h7/24.0

    ### find model data for each period
    i2 = np.where(np.logical_and(data1['time_hrly'] >= p2[0], data1['time_hrly'] <= p2[-1]))
    i3 = np.where(np.logical_and(data1['time_hrly'] >= p3[0], data1['time_hrly'] <= p3[-1]))
    i4 = np.where(np.logical_and(data1['time_hrly'] >= p4[0], data1['time_hrly'] <= p4[-1]))
    i5 = np.where(np.logical_and(data1['time_hrly'] >= p5[0], data1['time_hrly'] <= p5[-1]))
    i6 = np.where(np.logical_and(data1['time_hrly'] >= p6[0], data1['time_hrly'] <= p6[-1]))
    i7 = np.where(np.logical_and(data1['time_hrly'] >= p7[0], data1['time_hrly'] <= p7[-1]))

    #### remove flagged data
    data2['mixing_length_for_momentum'][data2['mixing_length_for_momentum'] == -9999.0] = np.nan
    data2['tke'][data2['tke'] == -9999.0] = np.nan

    #### calculate Km
    data2['Km'] = data2['mixing_length_for_momentum'] * np.sqrt(data2['tke'])

    ## color limits
    cmax = 10.0
    cmin = 0.0
    h = 5e3

    ####### ----------------------------------------------------------------
    #######     FIGURE
    ####### ----------------------------------------------------------------
    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.figure(figsize=(14,10))
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.2)

    ax = plt.gca()
    plt.pcolormesh(data2['time_hrly'][i2], data2['height'], np.transpose(data2['radr_refl'][i2[0],:]),
        # vmin = cmin, vmax = cmax
        )
    plt.colorbar()
    plt.ylim([0,h])
    plt.title('(P2) ' + d2 + ': ' + str(h2[0]) + ' - ' + str(h2[-1]) + 'h')
    ax.set_xticks(p2)
    ax.set_xticklabels(h2)

    # plt.savefig('../FIGS/TKE/PeriodComparison_UM_CASIM-100_P1_Km.png', dpi = 600)
    plt.show()

def RMSE_analysis(data1, data2, data3, obs):

    from sklearn.metrics import mean_squared_error
    from math import sqrt

    ###---------------------------------------------------------------------------------------------
    ###         Calculate root mean square error -
    ###---------------------------------------------------------------------------------------------
    rms = sqrt(mean_squared_error(y_actual, y_predicted))

def check_Radiation(data1, data2, data3, data4, obs, doy, out_dir1):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY

        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('filtering bad/missing radiation data points in ship dataset from all data:')
    print ('')

    #################################################################
    ## sort out obs['obs_temp']ervations' timestamp
    #################################################################
    # 0: Tship / (1)                         (time: 2324)
    # 1: LWdice / (1)                        (time3: 1293)
    # 2: LWuice / (1)                        (time3: 1293)
    # 3: precip / (1)                        (time4: 2352)
    # 4: Tice / (1)                          (time1: 1296)
    # 5: SWdship / (1)                       (time2: 2348)
    # 6: LWdship / (1)                       (time2: 2348)
    # 7: SWdice / (1)                        (time3: 1293)
    # 8: SWuice / (1)                        (time3: 1293)

    datenums_radice = obs['obs_temp'].variables['time3'][:] ### radiation (ice station) on different timestep
    time_radice_all = calcTime_Mat2DOY(datenums_radice)

    datenums_radship = obs['obs_temp'].variables['time2'][:] ### radiation (ship) on different timestep
    time_radship_all = calcTime_Mat2DOY(datenums_radship)

    datenums_tice = obs['obs_temp'].variables['time1'][:] ### ice camp data on different timestep
    time_tice = calcTime_Mat2DOY(datenums_tice)

    # UM -> IFS comparisons:
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

    ##############################################################################
    #####   first, need to transfer to hourly data for like-for-like comparison
    ##############################################################################
    ###  ---- ICE DATA
    time_radice = time_radice_all[:-1:2]
    lwdice = obs['obs_temp'].variables['LWdice'][:]
    lwuice = obs['obs_temp'].variables['LWuice'][:]
    lwdmeanice = np.mean(lwdice[:-1].reshape(-1, 2), axis=1)
    lwumeanice = np.mean(lwuice[:-1].reshape(-1, 2), axis=1)
    swdice = obs['obs_temp'].variables['SWdice'][:]
    swuice = obs['obs_temp'].variables['SWuice'][:]
    swdmeanice = np.mean(swdice[:-1].reshape(-1, 2), axis=1)
    swumeanice = np.mean(swuice[:-1].reshape(-1, 2), axis=1)
    lwnetice = lwdice - lwuice
    lwnetmeanice = np.mean(lwnetice[:-1].reshape(-1, 2), axis=1) # lwdmeanice - lwumeanice
    swnetice = swdice - swuice
    swnetmeanice = np.mean(swnetice[:-1].reshape(-1, 2), axis=1) # swdmeanice - swumeanice
    #### CODE USED ORIGINIALLY IN plot_paperRadiation
    # time_radice = time_radice_all[:-1:2]
    # lwd = obs['obs_temp'].variables['LWdice'][:]
    # lwu = obs['obs_temp'].variables['LWuice'][:]
    # lwdmean = np.nanmean(lwd[:-1].reshape(-1, 2), axis=1)
    # lwumean = np.nanmean(lwu[:-1].reshape(-1, 2), axis=1)
    # swd = obs['obs_temp'].variables['SWdice'][:]
    # swu = obs['obs_temp'].variables['SWuice'][:]
    # swdmean = np.nanmean(swd[:-1].reshape(-1, 2), axis=1)
    # swumean = np.nanmean(swu[:-1].reshape(-1, 2), axis=1)
    # netLW = lwdmean - lwumean
    # netSW = swdmean - swumean

    ###  ---- SHIP DATA
    time_radship = time_radship_all[:-1:2]
    lwd = obs['obs_temp'].variables['LWdship'][:]
    lwdmean = np.mean(lwd[:].reshape(-1, 2), axis=1)
    swd = obs['obs_temp'].variables['SWdship'][:]
    swdmean = np.mean(swd[:].reshape(-1, 2), axis=1)
    netlw = obs['obs_temp'].variables['LWnetship'][:]
    netlwmean = np.mean(netlw[:].reshape(-1, 2), axis=1)
    netsw = obs['obs_temp'].variables['SWnetship'][:]
    netswmean = np.mean(netsw[:].reshape(-1, 2), axis=1)

    # plt.plot(time_radice, swdmeanice)
    # plt.plot(time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])], swdmean[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])])
    # plt.show()

    # print (swdmeanice.shape)
    # print (swdmean[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])].shape)

    # print (time_radice)
    # print (time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])])

    # plt.plot(time_radice[:-1], time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])])
    # plt.show()
    #
    # plt.plot(time_radice[:-1] / time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])])
    # plt.show()

    ### only look at timesteps where we have both ice and ship data
    time_radship_icetimes = time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])]

    # print (time_radice[59])
    # print (time_radship_icetimes[59])
    ### 233.13541666666666, is missing from time_radship_icetimes

    ##############################################################################
    ### build new times array for ship times and ship data including the missing timestamp
    ##############################################################################
    times_ship = np.zeros([np.size(time_radship) + 1])
    swd_ship = np.zeros([np.size(time_radship) + 1])
    swnet_ship = np.zeros([np.size(time_radship) + 1])
    lwd_ship = np.zeros([np.size(time_radship) + 1])
    lwnet_ship = np.zeros([np.size(time_radship) + 1])

    ### find the missing index
    badindex = np.where(time_radship == 233.09375)

    ### fill new arrays with original data, with nans in missing data index
    times_ship[:int(badindex[0] + 1)] = time_radship[:int(badindex[0] + 1)]
    times_ship[int(badindex[0]) + 1 ] = 233.13541666666666
    times_ship[int(badindex[0] + 2):] = time_radship[int(badindex[0] + 1):]

    swd_ship[:int(badindex[0] + 1)] = swdmean[:int(badindex[0] + 1)]
    swd_ship[int(badindex[0]) + 1 ] = np.nan
    swd_ship[int(badindex[0] + 2):] = swdmean[int(badindex[0] + 1):]

    lwd_ship[:int(badindex[0] + 1)] = lwdmean[:int(badindex[0] + 1)]
    lwd_ship[int(badindex[0]) + 1 ] = np.nan
    lwd_ship[int(badindex[0] + 2):] = lwdmean[int(badindex[0] + 1):]

    lwnet_ship[:int(badindex[0] + 1)] = netlwmean[:int(badindex[0] + 1)]
    lwnet_ship[int(badindex[0]) + 1 ] = np.nan
    lwnet_ship[int(badindex[0] + 2):] = netlwmean[int(badindex[0] + 1):]

    swnet_ship[:int(badindex[0] + 1)] = netswmean[:int(badindex[0] + 1)]
    swnet_ship[int(badindex[0]) + 1 ] = np.nan
    swnet_ship[int(badindex[0] + 2):] = netswmean[int(badindex[0] + 1):]

    # plt.plot(time_radice[:] / times_ship[np.logical_and(times_ship >= time_radice[0], times_ship <= time_radice[-1])])
    # plt.show()
    # print (times_ship)
    # print (np.nanmax(times_ship))

    ###------------------------------------------------------------------------------------------------
    ### reassign sea ice observations (for fixing)
    obs['fixed_radiation'] = {}
    obs['fixed_radiation']['time_ice'] = time_radice
    obs['fixed_radiation']['LWnet_ice'] = lwnetmeanice
    obs['fixed_radiation']['SWnet_ice'] = swnetmeanice
    obs['fixed_radiation']['LWd_ice'] = lwdmeanice
    obs['fixed_radiation']['SWd_ice'] = swdmeanice
    obs['fixed_radiation']['LWu_ice'] = lwumeanice
    obs['fixed_radiation']['SWu_ice'] = swumeanice

    ##############################################################################
    ### pick out period where we have sea ice radiation measurements and set missing ship points to nans in ice array
    ##############################################################################
    iceindex = np.where(np.logical_and(times_ship >= time_radice[0], times_ship <= time_radice[-1]))

    swd_badpoints = np.isnan(swd_ship[iceindex[0]])
    obs['fixed_radiation']['SWd_ice'][swd_badpoints] = np.nan

    swnet_badpoints = np.isnan(swnet_ship[iceindex[0]])
    obs['fixed_radiation']['SWnet_ice'][swnet_badpoints] = np.nan

    lwd_badpoints = np.isnan(lwd_ship[iceindex[0]])
    obs['fixed_radiation']['LWd_ice'][lwd_badpoints] = np.nan

    lwnet_badpoints = np.isnan(lwnet_ship[iceindex[0]])
    obs['fixed_radiation']['LWnet_ice'][lwnet_badpoints] = np.nan

    ### reassign ship observations (fixed)
    obs['fixed_radiation']['time_ship'] = times_ship
    obs['fixed_radiation']['LWnet_ship'] = lwnet_ship
    obs['fixed_radiation']['SWnet_ship'] = swnet_ship
    obs['fixed_radiation']['LWd_ship'] = lwd_ship
    obs['fixed_radiation']['SWd_ship'] = swd_ship

    ##############################################################################
    ### now do the same for the model data
    ##############################################################################
    ### assign fixed_radiation dictionaries
    data1['fixed_radiation'] = {}
    data2['fixed_radiation'] = {}
    data3['fixed_radiation'] = {}
    data4['fixed_radiation'] = {}

    if out_dir1[-10:-5] == 'RadPA':
        modelindex = np.where(np.logical_and(obs['fixed_radiation']['time_ship'] >= data2['time_hrly'][0],
                obs['fixed_radiation']['time_ship'] <= data2['time_hrly'][-3]))

        # print (obs['fixed_radiation']['time_ship'][modelindex[0][0]])
        # print (data1['time_hrly'][0])
        # print (data2['time_hrly'][0])
        # print (data3['time_hrly'][0])
        # print (data4['time_hrly'][0])
        #
        # print (obs['fixed_radiation']['time_ship'][modelindex[0]].shape)
        # print (data1['time_hrly'].shape)
        # print (data2['time_hrly'].shape)
        # print (data3['time_hrly'].shape)
        # print (data4['time_hrly'].shape)
        #
        # print (obs['fixed_radiation']['time_ship'][modelindex[0][-1]])
        # print (data1['time_hrly'][-2])
        # print (data2['time_hrly'][-3])
        # print (data3['time_hrly'][-3])
        # print (data4['time_hrly'][-3])
        #
        # plt.plot(obs['fixed_radiation']['time_ship'][modelindex[0]])
        # plt.plot(data1['time_hrly'][:-2])
        # plt.plot(data2['time_hrly'][:-3])
        # plt.plot(data3['time_hrly'][:-3]);plt.show()

        data1['fixed_radiation']['time'] = data1['time_hrly'][:-2]
        data2['fixed_radiation']['time'] = data2['time_hrly'][:-3]
        data3['fixed_radiation']['time'] = data3['time_hrly'][:-4]
        data4['fixed_radiation']['time'] = data4['time_hrly'][:-3]

        model_swd_badpoints = np.isnan(swd_ship[modelindex[0]])
        data1['fixed_radiation']['SWd'] = data1['surface_downwelling_SW_radiation'][data1['hrly_flag']][:-2]
        data1['fixed_radiation']['SWd'][model_swd_badpoints] = np.nan
        data2['fixed_radiation']['SWd'] = data2['surface_downwelling_SW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['SWd'][model_swd_badpoints] = np.nan
        data3['fixed_radiation']['SWd'] = data3['sfc_down_sw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['SWd'][model_swd_badpoints] = np.nan
        data4['fixed_radiation']['SWd'] = data4['surface_downwelling_SW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['SWd'][model_swd_badpoints] = np.nan

        model_swnet_badpoints = np.isnan(swnet_ship[modelindex[0]])
        data1['fixed_radiation']['SWnet'] = data1['surface_net_SW_radiation'][data1['hrly_flag']][:-2]
        data1['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data2['fixed_radiation']['SWnet'] = data2['surface_net_SW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data3['fixed_radiation']['SWnet'] = data3['sfc_net_sw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data4['fixed_radiation']['SWnet'] = data4['surface_net_SW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan

        model_lwd_badpoints = np.isnan(lwd_ship[modelindex[0]])
        data1['fixed_radiation']['LWd'] = data1['surface_downwelling_LW_radiation'][data1['hrly_flag']][:-2]
        data1['fixed_radiation']['LWd'][model_lwd_badpoints] = np.nan
        data2['fixed_radiation']['LWd'] = data2['surface_downwelling_LW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['LWd'][model_lwd_badpoints] = np.nan
        data3['fixed_radiation']['LWd'] = data3['sfc_down_lw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['LWd'][model_lwd_badpoints] = np.nan
        data4['fixed_radiation']['LWd'] = data4['surface_downwelling_LW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['LWd'][model_lwd_badpoints] = np.nan

        model_lwnet_badpoints = np.isnan(lwnet_ship[modelindex[0]])
        data1['fixed_radiation']['LWnet'] = data1['surface_net_LW_radiation'][data1['hrly_flag']][:-2]
        data1['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data2['fixed_radiation']['LWnet'] = data2['surface_net_LW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data3['fixed_radiation']['LWnet'] = data3['sfc_net_lw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data4['fixed_radiation']['LWnet'] = data4['surface_net_LW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan

    else:
        modelindex = np.where(np.logical_and(obs['fixed_radiation']['time_ship'] >= data1['time_hrly'][0],
                obs['fixed_radiation']['time_ship'] <= data1['time_hrly'][-3]))

        data1['fixed_radiation']['time'] = data1['time_hrly'][:-3]
        data2['fixed_radiation']['time'] = data2['time_hrly'][:-3]
        data3['fixed_radiation']['time'] = data3['time_hrly'][:-3]
        data4['fixed_radiation']['time'] = data4['time_hrly'][:-3]

        print(data2['fixed_radiation']['time'].shape)
        print(data4['fixed_radiation']['time'].shape)
        print(data3['fixed_radiation']['time'].shape)

        model_swnet_badpoints = np.isnan(swnet_ship[modelindex[0]])
        data1['fixed_radiation']['SWnet'] = data1['surface_net_SW_radiation'][data1['hrly_flag']][:-3]
        data1['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data2['fixed_radiation']['SWnet'] = data2['surface_net_SW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data3['fixed_radiation']['SWnet'] = data3['sfc_net_sw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan
        data4['fixed_radiation']['SWnet'] = data4['surface_net_SW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['SWnet'][model_swnet_badpoints] = np.nan

        model_lwnet_badpoints = np.isnan(lwnet_ship[modelindex[0]])
        data1['fixed_radiation']['LWnet'] = data1['surface_net_LW_radiation'][data1['hrly_flag']][:-3]
        data1['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data2['fixed_radiation']['LWnet'] = data2['surface_net_LW_radiation'][data2['hrly_flag']][:-3]
        data2['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data3['fixed_radiation']['LWnet'] = data3['sfc_net_lw'][data3['hrly_flag']][:-3]
        data3['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan
        data4['fixed_radiation']['LWnet'] = data4['surface_net_LW_radiation'][data4['hrly_flag']][:-3]
        data4['fixed_radiation']['LWnet'][model_lwnet_badpoints] = np.nan

    # plt.plot(time_radice, swdmeanice)
    # plt.plot(time_radship[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])], swdmean[np.logical_and(time_radship >= time_radice[0], time_radship <= time_radice[-1])])
    # plt.plot(obs['fixed_radiation']['time_ice'],obs['fixed_radiation']['SWd_ice'])
    # plt.plot(obs['fixed_radiation']['time_ship'],obs['fixed_radiation']['SWd_ship'])
    # plt.show()
    #
    # plt.plot(time_radice_all, swnetice)
    # plt.plot(time_radship_all, netsw)
    # plt.plot(time_radice, swnetmeanice)
    # plt.plot(time_radship, netswmean)
    # # plt.plot(obs['fixed_radiation']['time_ice'],obs['fixed_radiation']['SWnet_ice'])
    # # plt.plot(obs['fixed_radiation']['time_ship'],obs['fixed_radiation']['SWnet_ship'])
    # plt.show()

    return data1, data2, data3, data4, obs

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')


    ### PYTHON3 FIXES
    # import numpy as np
    # np = py3_FixNPLoad(np)

    ### CHOOSE PLATFORM (OPTIONS BELOW)
    platform = 'LAPTOP'

    ### only works on laptop for now

    ### JASMIN
    ### LAPTOP
    ### MONSOON
    ### DESKTOP

    if platform == 'JASMIN':
        platformflag = 'jasmin'
        ship_filename = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
        um_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        misc_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/processed_models/'
        obs_root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MOCCHA/ODEN/DATA/'
    if platform == 'LAPTOP':
        platformflag = 'laptop'
        um_root_dir = '/home/gillian/MOCCHA/UM/DATA/'
        misc_root_dir = '/home/gillian/MOCCHA/ECMWF/'
        obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        ship_filename = '~/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'
    if platform == 'MONSOON':
        root_dir = '~/cylc-run/u-bg610/share/cycle/20160401T0000Z/HighArctic/1p5km/RA2M_CON/um/'
    if platform == 'DESKTOP':
        root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        ship_filename = '/nfs/a96/MOCCHA/working/gillian/ship/2018_shipposition_1hour.txt'
        position_filename = 'AUX_DATA/POSITION_UNROTATED.csv'

    ### CHOSEN RUN
    if platform == 'LAPTOP':
        out_dir1 = '25_u-cc568_RA2M_CON/OUT_R0/'
        out_dir2 = '23_u-cc278_RA1M_CASIM/OUT_R0/'
        # out_dir3 = 'MET_DATA/'
        out_dir3 = 'OUT_25H/'
        out_dir4 = '24_u-cc324_RA2T_CON/OUT_R0_LAM/' #'14_u-bu570_RA1M_CASIM/OUT_R1/' # '12_u-br210_RA1M_CASIM/OUT_R1/' #
        out_dir5 = '7_u-bn068_RA2T_CON/OUT_R2_glm/'
    elif platform == 'JASMIN':
        out_dir1 = 'UM_RA2M/'
        out_dir2 = 'UM_CASIM-100/'
        # out_dir3 = 'MET_DATA/'
        out_dir3 = 'ECMWF_IFS/'

    ### IFS: OUT_25H/
    ### 4_u-bg610_RA2M_CON/OUT_R1(_RadPA_25h)/
    ### 5_u-bl661_RA1M_CASIM/OUT_R0/            # 100/cc accum mode aerosol
    ### 6_u-bm410_RA1M_CASIM/                   # 200/cc accum mode aerosol
    ### 7_u-bn068_RA2T_CON/OUT_R2(R3_lam)\(_RadPA_25h)/              # RA2T_CON nest + global 4D stash
    ### 7_u-bn068_RA2T_CON/OUT_R2_glm/              # RA2T_CON nest + global 4D stash
    ### 8_u-bp738_RA2M_CON/OUT_R0/              # ERAI
    ### 10_u-bq791_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Fletcher Nice param
    ### 11_u-bq798_RA1M_CASIM/OUT_R0/      # CASIM with 100/cc accum mode soluble aerosol w/Meyers Nice param
    ### 12_u-br210_RA1M_CASIM/OUT_R1/           # UKCA daily averaged aerosol profiles, identical suite = u-bm507
    ### 13_u-br409_RA1M_CASIM/OUT_R0/           # 100/cc accum mode aerosol; ARG + Cooper; passive aerosol processing
    ### 14_u-bu570_RA1M_CASIM/OUT_R1/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit
    ### 15_u-bu687_RA2M_CON/OUT_24h/           # Wilson and Ballard 1999 uphys; new RHcrit
    ### 16_u-bv926_RA2T_CON/              # RA2T_CON nest + global 4D stash + no subgrid mp production
    ### 17_u-bz429_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics
    ### 18_u-ca011_RA1M_CASIM/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; 15 min res 3D diagnostics; 1A BL Scheme
    ### 19_u-ca012_RA2T_CON/              # RA2T_CON nest + global 4D stash; includes diagnosed turbulent dissipation rate
    ### 20_u-ca362_RA1M_CASIM/OUT_R0/           # 100/cc accum mode aerosol; ARG + Cooper; new RHcrit; CICE sea ice scheme
    ### 23_u-cc278_RA1M_CASIM/OUT_R0/           # 100/cc accum mde aerosol; ARG + Cooper; new RHcrit; sea ice albedo options as GLM
    ### 24_u-cc324_RA2T_CON/OUT_R0_LAM/             # RA2T_CON nest + global 4D stash. sea ice albedo (GLM+LAM) and extra BL diags (LAM) included
    ### 25_u-cc568_RA2M_CON/OUT_R0/             # Wilson and Ballard 1999 uphys. sea ice albedo and extra BL diags

    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')

    # -------------------------------------------------------------
    # Load ship track
    # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Load in ship track file:')
    print ('')
    ship_data = readfile(ship_filename)
    columns = assignColumns(ship_data)

    # -------------------------------------------------------------
    # Load observations
    # -------------------------------------------------------------
    print ('Loading observations:')
            # -------------------------------------------------------------
            # Which file does what?
            # -------------------------------------------------------------
            #### ice station: net LW / net SW
                    #### obs['ice_station_fluxes']/mast_radiation_30min_v2.3.mat
            #### obs['foremast']:
                    #### obs['foremast']/ACAS_AO2018_obs['foremast']_30min_v2_0.nc
            #### 7th deck: temperature, surface temperature, RH, downwelling SW, downwelling LW
                    #### 7thDeck/ACAS_AO2018_WX_30min_v2_0.nc

    obs = {}

    if platform == 'LAPTOP':
        print ('Load temporary ice station data from Jutta...')
        obs['obs_temp'] = Dataset(obs_root_dir + 'MET_DATA/MetData_Gillian_V3_30minres.nc','r')
        print ('Load ice station flux data from Jutta...')
        obs['ice_station_fluxes'] = readMatlabStruct(obs_root_dir + 'ice_station/flux30qc_trhwxrel.mat')

        print ('Load HATPRO data used by Cloudnet...')
        dirname = '/home/gillian/MOCCHA/ODEN/DATA/hatpro/'
        dir = os.listdir(dirname)
        obs['hatpro'] = {}
        for file in dir:
            IWVtemp = readMatlabStruct(dirname + file)
            print (IWVtemp.keys())
            if file == '20180814_IWV_30s_V2.mat':       ### if it is the first file
                obs['hatpro']['IWV'] = np.squeeze(IWVtemp['IWV'])
                obs['hatpro']['mday'] = np.squeeze(IWVtemp['mday'])
                obs['hatpro']['LWP'] = np.squeeze(IWVtemp['IWV'])
                obs['hatpro']['rainflag'] = np.squeeze(IWVtemp['rainflag'])
            else:
                obs['hatpro']['IWV'] = np.append(np.squeeze(obs['hatpro']['IWV']),np.squeeze(IWVtemp['IWV']))
                obs['hatpro']['mday'] = np.append(np.squeeze(obs['hatpro']['mday']),np.squeeze(IWVtemp['mday']))
                obs['hatpro']['LWP'] = np.append(np.squeeze(obs['hatpro']['LWP']),np.squeeze(IWVtemp['LWP']))
                obs['hatpro']['rainflag'] = np.append(np.squeeze(obs['hatpro']['rainflag']),np.squeeze(IWVtemp['rainflag']))
        obs['hatpro']['doy'] = calcTime_Mat2DOY(obs['hatpro']['mday'])

        print ('Load albedo estimates from Michael...')
        obs['albedo'] = readMatlabStruct(obs_root_dir + 'MOCCHA_Albedo_estimates_Michael.mat')

    ### print ('Load ice station radiation data from Jutta...')
    ### obs['ice_station_radiation'] = readMatlabStruct(obs_root_dir + 'ice_station/mast_radiation_30min_v2.3.mat')

    print ('Load radiosonde data from Jutta...')
    obs['sondes'] = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print ('Load observations inversion height data from Jutta...')
    obs['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/InversionHeights_RSh05int_final_V03.mat')

    print ('Load foremast data from John...')
    obs['foremast'] = Dataset(obs_root_dir + 'foremast/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    print ('Load 7th deck weather station data from John...')
    obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print ('Load weather sensor data from John...')
    obs['pws'] = readMatlabStruct(obs_root_dir + '7thDeck/ACAS_AO2018_PWD_30min_v1_0.mat')

    print ('...')

    # # -------------------------------------------------------------
    # # Load cube
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin cube read in at ' + time.strftime("%c"))
    print (' ')

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------
    # tempnames = ['umnsaa_pa012_r0.nc','umnsaa_pb012_r0.nc','umnsaa_pc011_r0.nc','umnsaa_pd011_r0.nc','20180812_oden_metum.nc']
    Aug_names = ['20180813_oden_','20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_']

    Sep_names = ['20180901_oden_','20180902_oden_','20180903_oden_','20180904_oden_',
            '20180905_oden_','20180906_oden_','20180907_oden_','20180908_oden_',
            '20180909_oden_','20180910_oden_','20180911_oden_','20180912_oden_',
            '20180913_oden_','20180914_oden_']

    moccha_names = ['20180814_oden_','20180815_oden_','20180816_oden_',
            '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
            '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
            '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
            '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
            '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
            '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
            '20180910_oden_','20180911_oden_','20180912_oden_','20180913_oden_','20180914_oden_']

    Aug_missing_files = []

    Sep_missing_files = []

    moccha_missing_files = ['20180813_oden_','20180910_oden_']   ### cloud radar not working    #,'20180914_oden_'
    missing_files = [225,253]    # manually set missing files doy for now ## 230, , 257

    doy = np.arange(226,259)        ## set DOY for full drift figures (over which we have cloudnet data)
    # doy = np.arange(226,258)        ## exclude 2019014 for RadPA files
    # doy = np.arange(240,251)        ## set DOY for subset of drift figures (presentations)
    # doy = np.arange(240,248)        ## set DOY for UM_CASIM-100_CICE  (28th Aug to 4th Sep)
    # doy = np.arange(243,250)        ## set DOY for ERAI-GLM  (31st Aug to 5th Sep)
    # doy = np.arange(243,250)        ## set DOY for CASIM_Nice tests  (31st Aug to 5th Sep)
    # doy = np.arange(226,259)        ## set DOY for CASIM-AeroProf (14th Aug to 14th Sep)
    # doy = np.arange(226,259)        ## set DOY for CASIM-100_AP (1st Sep to 14th Sep)
    # doy = np.arange(244,250)        ## set DOY for UM_RA2T_noTurbMP (1st Sep to 5th Sep)
    # doy = np.arange(237,259)        ## set DOY for RA2M_newRHcrit (25th Aug to 14th Sep)

    # names = ['umnsaa_pa000','umnsaa_pc000.nc']       ### DEFAULT OUTPUT NAMES FOR TESTING

    ## Choose month:
    names = moccha_names
    # missing_files = moccha_missing_files
    month_flag = -1

    for i in range(0,len(names)):
        if out_dir1[-10:-5] == 'RadPA':
            filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum_a.nc'
            filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum_a.nc'
            if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
                print( '***IFS being compared***')
                ifs_flag = True
                filename_um3 = misc_root_dir + out_dir3 + names[i] + 'ecmwf.nc'
            else:
                print ('***IFS NOT being compared***')
                filename_um3 = um_root_dir + out_dir3 + names[i] + 'metum_a.nc'
                ifs_flag = False
            filename_um4 = um_root_dir + out_dir4 + names[i] + 'metum_a.nc'
            filename_um5 = um_root_dir + out_dir5 + names[i] + 'metum.nc'
        else:
            filename_um1 = um_root_dir + out_dir1 + names[i] + 'metum.nc'
            filename_um2 = um_root_dir + out_dir2 + names[i] + 'metum.nc'
            if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
                print( '***IFS being compared***')
                ifs_flag = True
                filename_um3 = misc_root_dir + out_dir3 + names[i] + 'ecmwf.nc'
            else:
                print ('***IFS NOT being compared***')
                filename_um3 = um_root_dir + out_dir3 + names[i] + 'metum.nc'
                ifs_flag = False
            filename_um4 = um_root_dir + out_dir4 + names[i] + 'metum.nc'
            filename_um5 = um_root_dir + out_dir5 + names[i] + 'metum.nc'
        print (filename_um1)
        print (filename_um2)
        print (filename_um3)
        print (filename_um4)
        print (filename_um5)
        print ('')

        #### LOAD DATA
        print( 'Loading first run diagnostics:')
        nc1 = Dataset(filename_um1,'r')
        print ('...')
        print( 'Loading second run diagnostics:')
        nc2 = Dataset(filename_um2,'r')
        print ('...')
        print ('Loading third run diagnostics:')
        nc3 = Dataset(filename_um3,'r')
        print ('...')
        print ('Loading fourth run diagnostics:')
        nc4 = Dataset(filename_um4,'r')
        print ('...')
        print ('Loading fifth run diagnostics:')
        nc5 = Dataset(filename_um5,'r')
        print ('...')
        # -------------------------------------------------------------
        # print 'i = ' + str(i)
        print ('')

        #### LOAD IN SPECIFIC DIAGNOSTICS
        # if out_dir == '4_u-bg610_RA2M_CON/OUT_R1/':
        if out_dir1[-10:-5] == 'RadPA':
            var_list1 = ['surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'toa_outgoing_longwave_flux','toa_incoming_shortwave_flux','toa_outgoing_shortwave_flux']
            var_list3 = ['sfc_net_sw','sfc_net_lw', 'sfc_down_lw', 'sfc_down_sw', 'sfc_albedo']
            var_list2 = var_list1
            var_list4 = ['surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'toa_outgoing_longwave_flux','toa_incoming_shortwave_flux','toa_outgoing_shortwave_flux','sfc_net_SW','sfc_net_LW','sfc_downwelling_SW',
                'sfc_downwelling_LW']
            var_list5 = ['cloud_fraction','qliq','qice','temperature','q']
        else:
            ### BASE UM RUNS (UM_RA2M/UM_RA2T)
            if out_dir1[:21] == '25_u-cc568_RA2M_CON':
                var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                    'sensible_heat_flux','latent_heat_flux','rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                    'cloud_fraction','radr_refl']#,'temp_1.5m']
            else:
                var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
                    'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                    'cloud_fraction','radr_refl']#,'temp_1.5m','latent_heat_flux']
            var_list4 = var_list1
            ### CASIM RUNS
            if np.logical_or(out_dir4[:21] == '12_u-br210_RA1M_CASIM',out_dir4[:21] == '14_u-bu570_RA1M_CASIM'):
                var_list2 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
                'air_temperature_at_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                'cloud_fraction','radr_refl','qnliq','qnice','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'toa_outgoing_longwave_flux','toa_incoming_shortwave_flux','toa_outgoing_shortwave_flux'] # , 'latent_heat_flux']
                var_list4 = var_list2
            else:
                var_list2 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
                'air_temperature_at_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','uwind','vwind','wwind',
                'cloud_fraction','radr_refl']#, 'tke']#,'latent_heat_flux']#,'tke'] #'qnliq','qnice','mixing_length_for_momentum','qice',
                #, 'latent_heat_flux']
            ### IFS DIAGS
            if ifs_flag: var_list3 = ['height','flx_height','temperature','sfc_net_sw','sfc_net_lw','sfc_down_lat_heat_flx','sfc_down_sens_heat_flx',
                'sfc_temp_2m','flx_ls_rain','flx_conv_rain','flx_ls_snow','flx_conv_snow','q','pressure','sfc_bl_height','uwind','vwind','wwind',
                'sfc_down_lw', 'sfc_down_sw', 'sfc_albedo']
            if not ifs_flag:
                if out_dir3[-4:-1] == 'glm':
                    var_list3 = ['cloud_fraction','qliq','qice']
                else:
                    var_list3 = var_list1
            ### GLM DIAGS
            var_list5 = ['cloud_fraction','qliq','qice','temperature','q']

        if i == 0:
            ## ------------------
            #### UM
            ## ------------------
            data1 = {}
            data2 = {}
            data3 = {}
            data4 = {}
            data5 = {}
            if month_flag == -1:
                time_um1 = doy[i] + (nc1.variables['forecast_time'][:]/24.0)
                time_um2 = doy[i] + (nc2.variables['forecast_time'][:]/24.0)
                if ifs_flag: time_um3 = doy[i] + (nc3.variables['time'][:]/24.0)
                if not ifs_flag: time_um3 = doy[i] + (nc3.variables['forecast_time'][:]/24.0)
                time_um4 = doy[i] + (nc4.variables['forecast_time'][:]/24.0)
                time_um5 = doy[i] + (nc5.variables['forecast_time'][:]/24.0)
            else:
                time_um1 = float(filename_um1[-16:-14]) + (nc1.variables['forecast_time'][:]/24.0)
                time_um2 = float(filename_um2[-16:-14]) + (nc2.variables['forecast_time'][:]/24.0)
                if ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['time'][:]/24.0)
                if not ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['forecast_time'][:]/24.0)
                time_um4 = float(filename_um4[-16:-14]) + (nc4.variables['forecast_time'][:]/24.0)
                time_um5 = float(filename_um5[-16:-14]) + (nc5.variables['forecast_time'][:]/24.0)

            ### define height arrays explicitly
            if out_dir1[-10:-5] != 'RadPA':
                data1['height'] = nc1.variables['height'][:]
                data2['height'] = nc2.variables['height'][:]
                if not ifs_flag: data3['height'] = nc3.variables['height'][:]
                data4['height'] = nc4.variables['height'][:]
                data5['height'] = nc5.variables['height'][:]

            print ('Starting on t=0 RA2M data:')
            for j in range(0,len(var_list1)):
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) >= 1:
                    data1[var_list1[j]] = nc1.variables[var_list1[j]][:]
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            print ('Starting on t=0 CASIM data:')
            for j in range(0,len(var_list2)):
                print (var_list2[j])
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
                    data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
            nc2.close()
            ## ------------------
            #### um3
            ## ------------------
            print ('Starting on t=0 IFS data:')
            for j in range(0,len(var_list3)):
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) >= 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
            nc3.close()
            ## ------------------
            #### um4
            ## ------------------
            print ('Starting on t=0 RA2T data, initialising:')
            for j in range(0,len(var_list4)):
                # print (var_list4[j])
                if var_list4[j] in nc4.variables:
                    if np.ndim(nc4.variables[var_list4[j]]) == 0:     # ignore horizontal_resolution
                        continue
                    elif np.ndim(nc4.variables[var_list4[j]]) >= 1:
                        data4[var_list4[j]] = nc4.variables[var_list4[j]][:]
            nc4.close()
            ## ------------------
            #### um5
            ## ------------------
            print ('Starting on t=0 Global data:')
            for j in range(0,len(var_list5)):
                if np.ndim(nc5.variables[var_list5[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc5.variables[var_list5[j]]) >= 1:
                    data5[var_list5[j]] = nc5.variables[var_list5[j]][:]
            nc5.close()
            print ('')
        else:
            if month_flag == -1:
                time_um1 = np.append(time_um1, doy[i] + (nc1.variables['forecast_time'][:]/24.0))
                time_um2 = np.append(time_um2, doy[i] + (nc2.variables['forecast_time'][:]/24.0))
                if ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['time'][:]/24.0))
                if not ifs_flag: time_um3 = np.append(time_um3, doy[i] + (nc3.variables['forecast_time'][:]/24.0))
                time_um4 = np.append(time_um4, doy[i] + (nc4.variables['forecast_time'][:]/24.0))
                time_um5 = np.append(time_um5, doy[i] + (nc5.variables['forecast_time'][:]/24.0))
            ## ------------------
            #### UM
            ## ------------------
            print ('Appending UM data:')
            for j in range(0,len(var_list1)):
                # print (var_list1[j])
                if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc1.variables[var_list1[j]]) == 1:
                    # data1[cube_um1[j].var_name] = cube_um1[j].data
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:])
                elif np.ndim(nc1.variables[var_list1[j]]) == 2:
                    data1[var_list1[j]] = np.append(data1[var_list1[j]],nc1.variables[var_list1[j]][:],0)
            # np.save('working_data1',data1)
            nc1.close()
            ## ------------------
            #### um2
            ## ------------------
            print ('Appending CASIM data:')
            for j in range(0,len(var_list2)):
                print (var_list2[j])
                if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc2.variables[var_list2[j]]) == 1:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]],nc2.variables[var_list2[j]][:])
                elif np.ndim(nc2.variables[var_list2[j]]) == 2:
                    data2[var_list2[j]] = np.append(data2[var_list2[j]],nc2.variables[var_list2[j]][:],0)
            nc2.close()
            ## ------------------
            #### um3 / ifs
            ## ------------------
            print ('Appending IFS data:')
            for j in range(0,len(var_list3)):
                # print (var_list3[j])
                if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc3.variables[var_list3[j]]) == 1:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:])
                elif np.ndim(nc3.variables[var_list3[j]]) == 2:
                    data3[var_list3[j]] = np.append(data3[var_list3[j]],nc3.variables[var_list3[j]][:],0)
            nc3.close()
            ## ------------------
            #### um4 (ra2t)
            ## ------------------
            print ('Appending UM data:')
            for j in range(0,len(var_list4)):
                if var_list4[j] in nc4.variables:
                    print (var_list4[j])
                    if np.ndim(nc4.variables[var_list4[j]]) == 0:     # ignore horizontal_resolution
                        print ('ndim=0')
                        continue
                    elif np.ndim(nc4.variables[var_list4[j]]) == 1:
                        print ('ndim=1')
                        data4[var_list4[j]] = np.append(data4[var_list4[j]],nc4.variables[var_list4[j]][:])
                    elif np.ndim(nc4.variables[var_list4[j]]) == 2:
                        print ('ndim=2')
                        data4[var_list4[j]] = np.append(data4[var_list4[j]],nc4.variables[var_list4[j]][:],0)
            # np.save('working_data1',data1)
            nc4.close()
            print ('')
            ## ------------------
            #### um5 (global)
            ## ------------------
            print ('Appending UM_GLM data:')
            for j in range(0,len(var_list5)):
                # print (var_list5[j])
                if np.ndim(nc5.variables[var_list5[j]]) == 0:     # ignore horizontal_resolution
                    continue
                elif np.ndim(nc5.variables[var_list5[j]]) == 1:
                    data5[var_list5[j]] = np.append(data5[var_list5[j]],nc5.variables[var_list5[j]][:])
                elif np.ndim(nc5.variables[var_list5[j]]) == 2:
                    data5[var_list5[j]] = np.append(data5[var_list5[j]],nc5.variables[var_list5[j]][:],0)
            # np.save('working_data1',data1)
            nc5.close()
            print ('')

    #################################################################
    ## save time to dictionary now we're not looping over all diags anymore
    #################################################################
    data1['time'] = time_um1
    data2['time'] = time_um2
    data3['time'] = time_um3
    data4['time'] = time_um4
    data5['time'] = time_um5

    ### stop double counting of 0000 and 2400 from model data
    temp1 = np.zeros([len(data1['time'])])
    temp2 = np.zeros([len(data2['time'])])
    temp3 = np.zeros([len(data3['time'])])
    temp4 = np.zeros([len(data4['time'])])
    # temp5 = np.zeros([len(data5['time'])])
    for i in range(0, len(temp1)-1):
        if data1['time'][i] == data1['time'][i+1]:
            continue
        else:
            temp1[i] = data1['time'][i]
    for i in range(0, len(temp2)-1):
        if data2['time'][i] == data2['time'][i+1]:
            continue
        else:
            temp2[i] = data2['time'][i]
    for i in range(0, len(temp3)-1):
        if data3['time'][i] == data3['time'][i+1]:
            continue
        else:
            temp3[i] = data3['time'][i]
    for i in range(0, len(temp4)-1):
        if data4['time'][i] == data4['time'][i+1]:
            continue
        else:
            temp4[i] = data4['time'][i]
    # for i in range(0, len(temp5)-1):
    #     if data5['time'][i] == data5['time'][i+1]:
    #         continue
    #     else:
    #         temp5[i] = data5['time'][i]
    ii1 = np.where(temp1 != 0.0)      ### picks out where data are non-zero
    ii2 = np.where(temp2 != 0.0)      ### picks out where data are non-zero
    ii3 = np.where(temp3 != 0.0)      ### picks out where data are non-zero
    ii4 = np.where(temp4 != 0.0)      ### picks out where data are non-zero
    # ii5 = np.where(temp5 != 0.0)      ### picks out where data are non-zero (shouldn't be any - not cloudnet-ed)

    ### can use temp for all model data since they are on the same (hourly) time binning
    data1['time_hrly'] = temp1[ii1]
    data2['time_hrly'] = temp2[ii2]
    data3['time_hrly'] = temp3[ii3]
    data4['time_hrly'] = temp4[ii4]
    data1['time_6hrly'] = data1['time_hrly'][::6]
    data2['time_6hrly'] = data2['time_hrly'][::6]
    data3['time_6hrly'] = data3['time_hrly'][::6]
    data4['time_6hrly'] = data4['time_hrly'][::6]
    data1['hrly_flag'] = ii1
    data2['hrly_flag'] = ii2
    data3['hrly_flag'] = ii3
    data4['hrly_flag'] = ii4
    # data5['hrly_flag'] = ii5

    ### hard code glm time since not cloudnet-ed
    data5['time_hrly'] = data5['time']
    data5['hrly_flag'] = np.arange(len(data5['time_hrly']))

    #### add override for data2 to allow 24h data to be used for testing purposes
    if out_dir2[-4:] == '24h/':
        data2['time_hrly'] = data2['time']
        data2['hrly_flag'] = np.arange(len(data2['time_hrly']))

    #### add override for data2 to allow 24h data to be used for testing purposes
    # if out_dir4[:-2] == '12_u-br210_RA1M_CASIM/OUT_R':
    #     data4['time_hrly'] = data4['time']
    #     data4['hrly_flag'] = np.arange(len(data4['time_hrly']))


    #################################################################
    ## load calculated model inversion heights
    #################################################################
    print ('Load calculated model inversion heights...')
    data1['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_RA2M_inversion_results.mat')
    data2['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_CASIM-100_inversion_results.mat')
    data3['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/ECMWF_IFS_inversion_results.mat')

    #################################################################
    ## filter radiation measurements for bad/missing values
    #################################################################
    data1, data2, data3, data4, obs = check_Radiation(data1, data2, data3, data4, obs, doy, out_dir1)

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '25_u-cc568': label1 = 'UM_RA2M'
    if out_dir1[:10] == '24_u-cc324': label1 = 'UM_RA2T_' + out_dir1[-4:-1]
    if out_dir1[:10] == '23_u-cc278': label1 = 'UM_CASIM-100_GA6alb'
    if out_dir1[:10] == '20_u-ca326': label1 = 'UM_CASIM-100_CICE'
    if out_dir1[:10] == '16_u-bv926': label1 = 'UM_RA2T_noTurbMP'
    if out_dir1[:10] == '15_u-bu687': label1 = 'UM_RA2M_newRHcrit'
    if out_dir1[:10] == '14_u-bu570': label1 = 'UM_CASIM-100'
    if out_dir1[:10] == '13_u-br409': label1 = 'UM_CASIM-100_AP'
    if out_dir1[:10] == '12_u-br210': label1 = 'UM_CASIM-AeroProf'
    if out_dir1[:10] == '11_u-bq798': label1 = 'UM_CASIM-100_Meyers'
    if out_dir1[:10] == '10_u-bq791': label1 = 'UM_CASIM-100_Fletcher'
    if out_dir1[:9] == '8_u-bp738': label1 = 'UM_RA2M-ERAI-GLM'
    if out_dir1[:9] == '7_u-bn068': label1 = 'UM_RA2T_' + out_dir1[-4:-1]
    if out_dir1[:9] == '6_u-bm410': label1 = 'UM_CASIM-200'
    if out_dir1[:9] == '5_u-bl661': label1 = 'UM_CASIM-100'
    if out_dir1[:9] == '4_u-bg610': label1 = 'UM_RA2M'
    if out_dir1 == 'UM_RA2M/': label1 = 'UM_RA2M'

    label2 = 'undefined_label'
    if out_dir2[:10] == '25_u-cc568': label2 = 'UM_RA2M'
    if out_dir2[:10] == '24_u-cc324': label2 = 'UM_RA2T_' + out_dir2[-4:-1]
    if out_dir2[:10] == '23_u-cc278': label2 = 'UM_CASIM-100_GA6alb'
    if out_dir2[:10] == '20_u-ca326': label2 = 'UM_CASIM-100_CICE'
    if out_dir2[:10] == '16_u-bv926': label2 = 'UM_RA2T_noTurbMP'
    if out_dir2[:10] == '15_u-bu687': label2 = 'UM_RA2M_newRHcrit'
    if out_dir2[:10] == '14_u-bu570': label2 = 'UM_CASIM-100'
    if out_dir2[:10] == '13_u-br409': label2 = 'UM_CASIM-100_AP'
    if out_dir2[:10] == '12_u-br210': label2 = 'UM_CASIM-AeroProf'
    if out_dir2[:10] == '11_u-bq798': label2 = 'UM_CASIM-100_Meyers'
    if out_dir2[:10] == '10_u-bq791': label2 = 'UM_CASIM-100_Fletcher'
    if out_dir2[:9] == '8_u-bp738': label2 = 'UM_RA2M-ERAI-GLM'
    if out_dir2[:9] == '7_u-bn068': label2 = 'UM_RA2T_' + out_dir2[-4:-1]
    if out_dir2[:9] == '6_u-bm410': label2 = 'UM_CASIM-200'
    if out_dir2[:9] == '5_u-bl661': label2 = 'UM_CASIM-100'
    if out_dir2[:9] == '4_u-bg610': label2 = 'UM_RA2M'
    if out_dir2 == 'UM_CASIM-100/': label2 = 'UM_CASIM-100'

    label3 = 'undefined_label'
    if out_dir3 == 'OUT_25H/': label3 = 'ECMWF_IFS'
    if out_dir3[:10] == '25_u-cc568': label3 = 'UM_RA2M'
    if out_dir3[:10] == '24_u-cc324': label3 = 'UM_RA2T_' + out_dir3[-4:-1]
    if out_dir3[:10] == '23_u-cc278': label3 = 'UM_CASIM-100_GA6alb'
    if out_dir3[:10] == '20_u-ca326': label3 = 'UM_CASIM-100_CICE'
    if out_dir3[:10] == '16_u-bv926': label3 = 'UM_RA2T_noTurbMP'
    if out_dir3[:10] == '15_u-bu687': label3 = 'UM_RA2M_newRHcrit'
    if out_dir3[:10] == '14_u-bu570': label3 = 'UM_CASIM-100'
    if out_dir3[:10] == '13_u-br409': label3 = 'UM_CASIM-100_AP'
    if out_dir3[:10] == '12_u-br210': label3 = 'UM_CASIM-AeroProf'
    if out_dir3[:10] == '11_u-bq798': label3 = 'UM_CASIM-100_Meyers'
    if out_dir3[:10] == '10_u-bq791': label3 = 'UM_CASIM-100_Fletcher'
    if out_dir3[:9] == '8_u-bp738': label3 = 'UM_RA2M-ERAI-GLM'
    if out_dir3[:9] == '7_u-bn068': label3 = 'UM_RA2T_' + out_dir3[-4:-1]
    if out_dir3[:9] == '6_u-bm410': label3 = 'UM_CASIM-200'
    if out_dir3[:9] == '5_u-bl661': label3 = 'UM_CASIM-100'
    if out_dir3[:9] == '4_u-bg610': label3 = 'UM_RA2M'
    if out_dir3 == 'ECMWF_IFS/': label3 = 'ECMWF_IFS'

    label4 = 'undefined_label'
    if out_dir4[:10] == '25_u-cc568': label4 = 'UM_RA2M'
    if out_dir4[:10] == '24_u-cc324': label4 = 'UM_RA2T_' + out_dir4[-4:-1]
    if out_dir4[:10] == '23_u-cc278': label4 = 'UM_CASIM-100_GA6alb'
    if out_dir4[:10] == '20_u-ca326': label4 = 'UM_CASIM-100_CICE'
    if out_dir4[:10] == '16_u-bv926': label4 = 'UM_RA2T_noTurbMP'
    if out_dir4[:10] == '15_u-bu687': label4 = 'UM_RA2M_newRHcrit'
    if out_dir4[:10] == '14_u-bu570': label4 = 'UM_CASIM-100'
    if out_dir4[:10] == '13_u-br409': label4 = 'UM_CASIM-100_AP'
    if out_dir4[:10] == '12_u-br210': label4 = 'UM_CASIM-AeroProf'
    if out_dir4[:10] == '11_u-bq798': label4 = 'UM_CASIM-100_Meyers'
    if out_dir4[:10] == '10_u-bq791': label4 = 'UM_CASIM-100_Fletcher'
    if out_dir4[:9] == '8_u-bp738': label4 = 'UM_ERAI-GLM'
    if out_dir4[:9] == '7_u-bn068': label4 = 'UM_RA2T_' + out_dir4[-4:-1]
    if out_dir4[:9] == '6_u-bm410': label4 = 'UM_CASIM-200'
    if out_dir4[:9] == '5_u-bl661': label4 = 'UM_CASIM-100'
    if out_dir4[:9] == '4_u-bg610': label4 = 'UM_RA2M'
    if out_dir4 == 'UM_RA2M/': label4 = 'UM_RA2M'

    label5 = 'undefined_label'
    if out_dir5[-4:-1] == 'glm': label5 = 'UM_GLM'

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_data4', data4)
    np.save('working_data5', data5)
    np.save('working_dataObs', obs['sondes'])

    # -------------------------------------------------------------
    # Plot combined timeseries as lineplot
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)

    # figure = plot_line_BLDepth(time_um1, time_um2, data1, data2, cube_um1, cube_um2, month_flag,
    #             missing_files, out_dir1, obs, doy)

    # figure = plot_line_RAD(data1, data2, data3, cube_um1, cube_um2, cube_um3,
    #     month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # CASIM plots
    # -------------------------------------------------------------
    # figure = plot_line_CASIM_NiceTest(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_CASIM_NdropTimeseries(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_CASIM_NiceTimeseries(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_CASIM_QliqTimeseries(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # Plot paper figures
    # -------------------------------------------------------------
    # figure = plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    figure = plot_paperRadiation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # figure = plot_Precipitation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # figure = plot_BLDepth(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_BLType(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_paperGLMAnalysis(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_paperRadiosondes(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_paperERAIProfiles(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_paperCASIMNiceProfiles(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_RadiosondesTemperature(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # figure = plot_RadiosondesQ(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # figure = period_Selection(data1, data2, data3, data4, data5, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4, label5)
    # figure = plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesTheta(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_line_RA2T(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_Cv_RA2T(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, out_dir4, obs, doy, label1, label2, label3, label4)
    # figure = plot_CWC_RA2T(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # figure = plot_line_subSect(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plotWinds(data1, data2, data3, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # Further analysis
    # -------------------------------------------------------------
    # data1, data2, data3, obs = inversionIdent(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # out = table_Radiation(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
    # out = table_Fluxes(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)
            ### need to use run #5 instead of run #14 for data2
    # out = plot_sfcAlbedo(data1, data2, data3, data4, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, label4)

    # -------------------------------------------------------------
    # -------------------------------------------------------------
    # -------------------------------------------------------------
    # Sandeep - MOCCHA TKE dissipation rates
    # -------------------------------------------------------------
    # out = radarRefl_Sandeep(data1, data2, data3, data4, obs, doy, label1, label2, label3, label4)
    # out = TKEDissRate_Sandeep(data1, data2, data3, data4, obs, doy, label1, label2, label3, label4)

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)
    np.save('working_data4', data4)
    np.save('working_data5', data5)
    np.save('working_dataObs', obs['sondes'])

    # -------------------------------------------------------------
    # save out radiosonde biases
    # -------------------------------------------------------------
    # np.save('../Cloudnet/UM_RA2M_SondeBiases', data1)
    # np.save('../Cloudnet/UM_CASIM-100_SondeBiases', data2)
    # np.save('../Cloudnet/ECMWF_IFS_SondeBiases', data3)
    # np.save('../Cloudnet/UM_RA2T_SondeBiases', data4)
    # np.save('../Cloudnet/UM_GLM_SondeBiases', data5)
    # np.save('../Cloudnet/SondeData_reGridded', obs['sondes'])

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')

    #### DIAGNOSTICS TO CHOOSE FROM:

    # UM -> um2 comparisons:
    # 1. snowfall_flux -> sfc_ls_snow
    # 2. rainfall_flux -> sfc_ls_rain
    # 3. sensible_heat_flux -> sfc_down_sens_heat_flx
    # 4. latent_heat_flux -> flx_turb_moist
    # 5. bl_depth -> sfc_bl_height
    # 6. sfc_pressure -> sfc_pressure
    # 7. sfc_temperature -> sfc_temperature
    # 8. surface_net_LW_radiation -> sfc_net_lw
    # 9. surface_net_SW_radiation -> sfc_net_sw

if __name__ == '__main__':

    main()
