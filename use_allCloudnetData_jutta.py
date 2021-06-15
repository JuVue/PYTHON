### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
### loads in the Cloudnet data and the raw model data for the cloudnet figures

from __future__ import print_function
import time
import datetime as dtime

import numpy as np
import pandas as pd
from netCDF4 import Dataset

# import diags_MOCCHA as diags
# import diags_varnames as varnames
#import cartopy.crs as ccrs
#import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.dates as mdates

import os
#import seaborn as sns
from scipy.interpolate import interp1d
from IPython import embed
#### import python functions
import sys
sys.path.insert(1, './py_functions/')
from time_functions import calcTime_Mat2DOY, date2datenum, datenum2date
from readMAT import readMatlabStruct
from physFuncts import calcThetaE, calcThetaVL
from manipFuncts import int2list
# from conversionFuncts import reGrid_Sondes

def readfile(filename):

    import pandas as pd

    # print '******'
    print( '')
    print( 'Reading .txt file with pandas')
    print( '')

    data = pd.read_csv(filename, sep = " ")
    values = data.values

    return data

#def assignColumns(data):

    columns = ['Year', 'Month', 'Day', 'Hour', 'Minutes', 'Seconds', 'Longitude', 'Latitude']

    return columns

#def iceDrift(data):

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

#def inIce(data):

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

#def trackShip(data, date):
    ###################################
    ## DEFINE METUM PERIOD (CLOUDNET COMPARISON)
    ###################################
    trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==14,data.values[:,1]==8),data.values[:,3]>=0))
    trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==25,data.values[:,1]==8),data.values[:,3]==1))
    # trackShip_start = np.where(np.logical_and(np.logical_and(data.values[:,2]==int(date[-2:]),data.values[:,1]==int(date[-4:-2])),data.values[:,3]>=0))
    # trackShip_end = np.where(np.logical_and(np.logical_and(data.values[:,2]==(int(date[-2:]) + 1),data.values[:,1]==int(date[-4:-2])),data.values[:,3]==1))
    trackShip_index = range(trackShip_start[0][0],trackShip_end[0][-1])

    print( '******')
    print( '')
    # print 'Mean lon/lat of ship track: (' + str(np.nanmedian(data.values[inIce_index,6])) + ', ' + str(np.nanmedian(data.values[inIce_index,7])) + ')'
    print( 'Lon/lat of start point: (' + str(data.values[trackShip_index[0],6]) + ', ' + str(data.values[trackShip_index[0],7]) + ')')
    print( 'Lon/lat of end point: (' + str(data.values[trackShip_index[-1],6]) + ', ' + str(data.values[trackShip_index[-1],7]) + ')')
    # print 'Start: ' + str(data.values[trackShip_start[0][0],0:4])
    # print 'End: ' + str(data.values[trackShip_end[0][-1],0:4])
    print( 'trackShip: ' + str(data.values[trackShip_index[0],0:4]) + ' - ' + str(data.values[trackShip_index[-1],0:4]))
    print( '')

    return trackShip_index

#def plot_CvProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs, obs_switch): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting Cv statistics for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	Plot mean profiles
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
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax = plt.gca()

    print (um_data.keys())

    #### set flagged um_data to nans
    # um_data['Cv'][um_data['Cv'] == -999] = np.nan
    # ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    # obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # obs_data['Cv_adv'][obs_data['Cv_adv'] == -999] = np.nan
    # # um_data['Cv'][um_data['Cv'] == 0] = np.nan
    # um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    # misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan
    # ra2t_data['model_Cv_filtered'][ra2t_data['model_Cv_filtered'] < 0.0] = np.nan

    if obs_switch == 'RADAR':
        # plt.plot(np.nanmean(obs_data['Cv'],0),obs_data['height'], 'k--', linewidth = 3, label = 'Obs-on-model-grid')
        # ax.fill_betweenx(obs_data['height'],np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
        #     np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
        # plt.plot(np.nanmean(obs_data['Cv_adv'],0),np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs')
        # plt.plot(np.nanmean(obs['cloudfractions']['cloudfraction_total'],0),np.squeeze(obs['cloudfractions']['height']),
        #     'k--', linewidth = 3, label = 'Obs-on-artificial-grid')
        ax.fill_betweenx(np.squeeze(obs['cloudfractions']['height']),np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) - np.nanstd(obs['cloudfractions']['cloudfraction_total'],0),
            np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) + np.nanstd(obs['cloudfractions']['cloudfraction_total'],0), color = 'lightgrey', alpha = 0.5)
    else:
        # plt.plot(np.nanmean(obs_data['Cv_adv'],0),np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs Cv_adv')
        plt.plot(np.nanmean(obs_data['Cv'],0),np.nanmean(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 5)
        ax.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
            np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
        # plt.plot(np.nanmean(obs['cloudfractions']['cloudfraction_total'],0),np.squeeze(obs['cloudfractions']['height']),
        #     '--', color = 'grey', linewidth = 3, label = 'Obs-on-artificial-grid')
        # ax.fill_betweenx(np.squeeze(obs['cloudfractions']['height']),np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) - np.nanstd(obs['cloudfractions']['cloudfraction_total'],0),
        #     np.nanmean(obs['cloudfractions']['cloudfraction_total'],0) + np.nanstd(obs['cloudfractions']['cloudfraction_total'],0), color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0), np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)

    ax.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0),
        np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0) - np.nanstd(ifs_data['model_snow_Cv_filtered'],0), np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0) + np.nanstd(ifs_data['model_snow_Cv_filtered'],0), np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)

    ax.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_Cv_filtered'],0) - np.nanstd(misc_data['model_Cv_filtered'],0),
        np.nanmean(misc_data['model_Cv_filtered'],0) + np.nanstd(misc_data['model_Cv_filtered'],0), color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0) - np.nanstd(misc_data['model_Cv_filtered'],0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0) + np.nanstd(misc_data['model_Cv_filtered'],0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(ra2t_data['model_Cv_filtered'],0) - np.nanstd(ra2t_data['model_Cv_filtered'],0),
        np.nanmean(ra2t_data['model_Cv_filtered'],0) + np.nanstd(ra2t_data['model_Cv_filtered'],0), color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmean(ra2t_data['model_Cv_filtered'],0) - np.nanstd(ra2t_data['model_Cv_filtered'],0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(ra2t_data['model_Cv_filtered'],0) + np.nanstd(ra2t_data['model_Cv_filtered'],0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0),
        np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0) - np.nanstd(um_data['model_Cv_filtered'],0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0) + np.nanstd(um_data['model_Cv_filtered'],0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(ifs_data['model_snow_Cv_filtered'],0),np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    plt.plot(np.nanmean(misc_data['model_Cv_filtered'],0),np.nanmean(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-AeroProf', zorder = 3)
    plt.plot(np.nanmean(ra2t_data['model_Cv_filtered'],0),np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(um_data['model_Cv_filtered'],0),np.nanmean(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    plt.xlabel('C$_{V}$')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(np.arange(0,9.01e3,0.5e3))
    ax.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    plt.xlim([0,1])
    plt.xticks(np.arange(0,1.01,0.1))
    ax.set_xticklabels([0,' ',0.2,' ',0.4,' ',0.6,' ',0.8,' ',1.0])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid-QF30_RA2M_IFS_CASIM-AeroProf_RA2T_Cv_226-257DOY_fixedRA2T_noOffsetLWP_wSetFlags.svg'
    plt.savefig(fileout)
    plt.show()

    print ('Z = ')
    print (np.nanmean(misc_data['height'],0))

    Zindex = np.where(np.nanmean(misc_data['height'],0) == 4.55500000e+03)
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_Cv_filtered'][:,Zindex[0]],0))
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_Cv_filtered'][:,Zindex[0]],0))
    print ('Obs = ')
    print (np.nanmean(obs_data['Cv'][:,Zindex[0]],0))

def plot_CvTimeseries(um_data, misc_data, ra2t_data, obs_data, dates, plots_out_dir, **args):

    if isinstance(args,dict):
        monc_data=args[list(args.keys())[0]]

    ylims=[0,5]
    yticks=np.arange(0,5e3,1e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting Cv timeseries for whole drift period:')
    print ('')


    SMALL_SIZE = 12
    MED_SIZE = 15
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(9.5,13))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    ### define axis instance
    ax = plt.gca()

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)

    plt.subplot(411)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['Cv']),
          np.arange(0,1.1,0.1),
          cmap = newcmp,
          zorder = 1)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['sfmlheight']), color = 'grey', linewidth = 1.0)

    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    #ax.xaxis.set_major_locator(mdates.DayLocator(interval=5))
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%d%m'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 35)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('C$_{V}$')
    # plt.colorbar()

    plt.subplot(412)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_Cv_filtered']),
         np.arange(0,1.1,0.1),
         cmap = newcmp,
         zorder = 1)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_CASIM', rotation = 270, labelpad = 17)
    ax2.set_yticks([])

    plt.subplot(413)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(ra2t_data['model_Cv_filtered']),
        np.arange(0,1.1,0.1),
        cmap = newcmp,
        zorder = 1)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data2['time_hrly'][::6], bldepth4[::6], 'k', linewidth = 1.0)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2T', rotation = 270, labelpad = 17)
    ax2.set_yticks([])


    plt.subplot(414)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_Cv_filtered']),
        np.arange(0,1.1,0.1),
        cmap = newcmp,
        zorder = 1)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.xlabel('Date')
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2M', rotation = 270, labelpad = 17)
    ax2.set_yticks([])

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')
    dstr=datenum2date(dates[1])
    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_CvTimeseries.png'
    plt.savefig(fileout)

def plot_LWCTimeseries(um_data, misc_data, ra2t_data, obs_data, plots_out_dir, dates, obs_switch): #, lon, lat):

    ylims=[0,5]
    yticks=np.arange(0,5e3,1e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting LWC timeseries for whole drift period:')
    print ('')

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(9.5,13))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    ### define axis instance
    ax = plt.gca()

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:1, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)


    ### set colour max as var
    cmax = 0.3

    plt.subplot(411)
    ax = plt.gca()
    if obs_switch == 'RADAR':
        img = plt.pcolor(obs_data['time'][::6], obs_data['height'][:394], np.transpose(obs_data['lwc'][::6,:394])*1e3,
            vmin = 0.0, vmax = cmax)
            #cmap = newcmp)
    else:
        img = plt.pcolor(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['lwc_adiabatic'])*1e3,
            cmap=newcmp,vmin = 0.0, vmax = cmax)
            #cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 35)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('LWC')

    plt.subplot(412)
    ax = plt.gca()
    plt.pcolor(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_lwc'])*1e3,
        cmap=newcmp,vmin = 0.0, vmax = cmax)
        # cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_CASIM', rotation = 270, labelpad = 17)
    ax2.set_yticks([])


    plt.subplot(413)
    ax = plt.gca()
    plt.pcolor(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(ra2t_data['model_lwc'])*1e3,
    cmap=newcmp,vmin = 0.0, vmax = cmax)
    # cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2T', rotation = 270, labelpad = 17)
    ax2.set_yticks([])


    plt.subplot(414)
    ax = plt.gca()
    plt.pcolor(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_lwc'])*1e3,
        cmap=newcmp,vmin = 0.0, vmax = cmax)
    # cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2M', rotation = 270, labelpad = 17)
    ax2.set_yticks([])

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')
    dstr=datenum2date(dates[1])

    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_LWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

def plot_IWCTimeseries(um_data,  misc_data, ra2t_data, obs_data, plots_out_dir, dates, obs_switch): #, lon, lat):

    ylims=[0,5]
    yticks=np.arange(0,5e3,1e3)
    ytlabels=yticks/1e3

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:1, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)

    print ('******')
    print ('')
    print ('Plotting IWC timeseries for whole drift period:')
    print ('')

    SMALL_SIZE = 12
    MED_SIZE = 12
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    fig = plt.figure(figsize=(9.5,13))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    ### define axis instance
    ax = plt.gca()

    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] == -999.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] == -999.0] = np.nan

    obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] <= 0.0] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] <= 0.0] = np.nan

    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    cmax = 0.05

    plt.subplot(411)
    ax = plt.gca()
    if obs_switch == 'RADAR':
        img = plt.pcolor(obs_data['time'], obs_data['height'], np.transpose(obs_data['iwc'])*1e3,
            cmap=newcmp,vmin = 0.0, vmax = cmax)
            #cmap = newcmp)
    else:
        img = plt.pcolor(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['iwc'])*1e3,
            cmap=newcmp,vmin = 0.0, vmax = cmax)
            #cmap = newcmp)

    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 35)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('IWC')


    plt.subplot(413)
    ax = plt.gca()
    plt.pcolor(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(ra2t_data['model_iwc_filtered'])*1e3,
        cmap=newcmp,vmin = 0.0, vmax = cmax)
        #cmap = newcmp)

    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2T', rotation = 270, labelpad = 35)
    ax2.set_yticks([])

    plt.subplot(412)
    ax = plt.gca()
    plt.pcolor(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_iwc_filtered'])*1e3,
        cmap=newcmp,vmin = 0.0, vmax = cmax)
        #cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_CASIM', rotation = 270, labelpad = 35)
    ax2.set_yticks([])

    plt.subplot(414)
    ax = plt.gca()
    plt.pcolor(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_iwc_filtered'])*1e3,
        cmap=newcmp,vmin = 0.0, vmax = cmax)
        #cmap = newcmp)

    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('UM_RA2M', rotation = 270, labelpad = 35)
    ax2.set_yticks([])

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')
    dstr=datenum2date(dates[1])

    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_IWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

def plot_TWCTimeseries(um_data, misc_data, ra2t_data, obs_data, plots_out_dir, dates, obs_switch,nanind,wcind):
    ylims=[0,5]
    yticks=np.arange(0,5e3,1e3)
    ytlabels=yticks/1e3

    from matplotlib.colors import LogNorm

    print ('******')
    print ('')
    print ('Plotting TWC timeseries for whole drift period:')
    print ('')

    ylims=[0,5]
    yticks=np.arange(0,5e3,1e3)
    ytlabels=yticks/1e3

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    # bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    # bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    # bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    # bldepth4 = data4['bl_depth'][data4['hrly_flag']]
    #
    twc0 = np.transpose(obs_data['twc'])*1e3

    cmax = 0.3

    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:10, :] = greyclr
    newcmp = ListedColormap(newcolors)

    SMALL_SIZE = 12
    MED_SIZE = 15
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)



    #####Plot Twc###############################################
    fig = plt.figure(figsize=(9.5,13))
    plt.subplots_adjust(top = 0.9, bottom = 0.06, right = 0.98, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(411)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), twc0,
        # np.arange(0,0.31,0.01),
        # locator=ticker.LogLocator(base = 10.0),
        levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
        cmap = newcmp)
        # )
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('Obs-' + obs_switch + 'grid')
    cbaxes = fig.add_axes([0.225, 0.95, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('TWC [g m$^{-3}$]')

    plt.subplot(412)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(misc_data['model_twc'])*1e3,
        # locator=ticker.LogLocator(base = 10.0),
        levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
        # np.arange(0,0.31,0.001),
        cmap = newcmp)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data2['time_hrly'][::6], bldepth2[::6], 'k', linewidth = 1.0)
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_CASIM')
    # plt.colorbar()

    plt.subplot(413)
    ax = plt.gca()
    # ax.set_facecolor('gainsboro')
    plt.contourf(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(ra2t_data['model_twc'])*1e3,
        # locator=ticker.LogLocator(base = 10.0),
        levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
        # np.arange(0,0.31,0.001),
        cmap = newcmp)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data4['time_hrly'][::6], bldepth4[::6], 'k', linewidth = 1.0)
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_RA2T')
    # plt.colorbar()

    plt.subplot(414)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_twc'])*1e3,
        # locator=ticker.LogLocator(base = 10.0),
        levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
        # np.arange(0,0.31,0.001),
        cmap = newcmp)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_RA2M')
    # plt.colorbar()

    plt.xlabel('Date')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    dstr=datenum2date(dates[1])
    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_TWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    ####            CREATE TWC MASK
    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------

    mask0 = np.zeros([np.size(obs_data['twc'],0), np.size(obs_data['twc'],1)])
    mask1 = np.zeros([np.size(um_data['model_twc'],0), np.size(um_data['model_twc'],1)])
    mask2 = np.zeros([np.size(misc_data['model_twc'],0), np.size(misc_data['model_twc'],1)])
    mask4 = np.zeros([np.size(ra2t_data['model_twc'],0), np.size(ra2t_data['model_twc'],1)])

    # print ('mask0 size = ' + str(mask0.shape))

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)


    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                elif obs_data['twc'][t,k] >= twc_thresh_um[k]:
                    mask0[t,k] = 1.0
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
            elif um_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask1[t,k] = 1.0
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
            elif misc_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask2[t,k] = 1.0
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
            elif ra2t_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask4[t,k] = 1.0


    mask0[nanind] = np.nan
    mask1[nanind] = np.nan
    mask2[nanind] = np.nan
    mask4[nanind] = np.nan

    mask0[wcind] = np.nan
    mask1[wcind] = np.nan
    mask2[wcind] = np.nan
    mask4[wcind] = np.nan



    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    fig = plt.figure(figsize=(9.5,13))
    plt.subplots_adjust(top = 0.9, bottom = 0.06, right = 0.98, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(411)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(mask0),
        np.arange(0,1.01,0.1),
        cmap = mpl_cm.viridis)
    #img = plt.pcolor(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(mask0),
    #            cmap = mpl_cm.viridis,vmin=0,vmax=0.1)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    ax = plt.gca()
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('Obs-' + obs_switch + 'grid')
    cbaxes = fig.add_axes([0.225, 0.95, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('Cloud mask')

    plt.subplot(412)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(mask2),
        np.arange(0,1.01,0.1),
        cmap = mpl_cm.viridis)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data2['time_hrly'][::6], bldepth2[::6], 'k', linewidth = 1.0)
    ax = plt.gca()
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_CASIM')
    # plt.colorbar()

    plt.subplot(413)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(mask4),
        np.arange(0,1.01,0.1),
        cmap = mpl_cm.viridis)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data4['time_hrly'][::6], bldepth4[::6], 'k', linewidth = 1.0)
    ax = plt.gca()
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_RA2T')
    # plt.colorbar()

    plt.subplot(414)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    plt.contourf(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(mask1),
        np.arange(0,1.01,0.1),
        cmap = mpl_cm.viridis)
    # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
    ax = plt.gca()
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('UM_RA2M')
    # plt.colorbar()
    plt.xlabel('Date')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_TWC-MASKTimeseries.png'
    plt.savefig(fileout)
    plt.close()

    # #####TCW MASK pcolor #############################################
    # fig = plt.figure(figsize=(9.5,13))
    # plt.subplots_adjust(top = 0.9, bottom = 0.06, right = 0.98, left = 0.08,
    #         hspace = 0.4, wspace = 0.2)
    #
    # plt.subplot(411)
    # ax = plt.gca()
    # img = plt.pcolor(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(mask0),
    #             cmap = mpl_cm.viridis)
    # ax = plt.gca()
    # nans = ax.get_ylim()
    # plt.ylabel('Z [km]')
    # plt.ylim(ylims)
    # plt.yticks(yticks)
    # ax.set_yticklabels(ytlabels)
    # plt.xlim([dates[0], dates[1]])
    # ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    # plt.title('Obs-' + obs_switch + 'grid')
    # cbaxes = fig.add_axes([0.225, 0.95, 0.6, 0.015])
    # cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    # plt.title('Cloud mask')
    #
    # plt.subplot(412)
    # ax = plt.gca()
    # plt.pcolor(misc_data['time'], np.squeeze(misc_data['height'][0,:]), np.transpose(mask2),
    #     cmap = mpl_cm.viridis)
    # ax = plt.gca()
    # nans = ax.get_ylim()
    # plt.ylabel('Z [km]')
    # plt.ylim(ylims)
    # plt.yticks(yticks)
    # ax.set_yticklabels(ytlabels)
    # plt.xlim([dates[0], dates[1]])
    # ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    # plt.title('UM_CASIM')
    # # plt.colorbar()
    #
    # plt.subplot(413)
    # ax = plt.gca()
    # # ax.set_facecolor('aliceblue')
    # plt.pcolor(ra2t_data['time'], np.squeeze(ra2t_data['height'][0,:]), np.transpose(mask4),
    #     cmap = mpl_cm.viridis)
    # # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # # plt.plot(data4['time_hrly'][::6], bldepth4[::6], 'k', linewidth = 1.0)
    # ax = plt.gca()
    # nans = ax.get_ylim()
    # plt.ylabel('Z [km]')
    # plt.ylim(ylims)
    # plt.yticks(yticks)
    # ax.set_yticklabels(ytlabels)
    # plt.xlim([dates[0], dates[1]])
    # ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    # plt.title('UM_RA2T')
    # # plt.colorbar()
    #
    # plt.subplot(414)
    # ax = plt.gca()
    # # ax.set_facecolor('aliceblue')
    # plt.pcolor(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(mask1),
    #     cmap = mpl_cm.viridis)
    # # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
    # # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
    # ax = plt.gca()
    # nans = ax.get_ylim()
    # plt.ylabel('Z [km]')
    # plt.ylim(ylims)
    # plt.yticks(yticks)
    # ax.set_yticklabels(ytlabels)
    # plt.xlim([dates[0], dates[1]])
    # ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    # ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    # plt.title('UM_RA2M')
    # # plt.colorbar()
    # plt.xlabel('Date')
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_TWC-MASKTimeseries_pcolor.png'
    # plt.savefig(fileout)
    # plt.close()
    #############################################################

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    plt.plot(np.nanmean(mask0,0),np.nanmean(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 5)
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(mask0,0) - np.nanstd(mask0,0),
        np.nanmean(mask0,0) + np.nanstd(mask0,0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(mask0,0) - np.nanstd(mask0,0), np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(mask0,0) + np.nanstd(mask0,0), np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(mask1,0) - np.nanstd(mask1,0),
        np.nanmean(mask1,0) + np.nanstd(mask1,0), color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(mask1,0) - np.nanstd(mask1,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(mask1,0) + np.nanstd(mask1,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(mask2,0) - np.nanstd(mask2,0),
        np.nanmean(mask2,0) + np.nanstd(mask2,0), color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(mask2,0) - np.nanstd(mask2,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(mask2,0) + np.nanstd(mask2,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(mask4,0) - np.nanstd(mask4,0),
        np.nanmean(mask4,0) + np.nanstd(mask4,0), color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(mask4,0) - np.nanstd(mask4,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(mask4,0) + np.nanstd(mask4,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(mask2,0),np.nanmean(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-AeroProf', zorder = 3)
    plt.plot(np.nanmean(mask4,0),np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(mask1,0),np.nanmean(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    # plt.plot(np.nanmedian(mask3,0),np.nanmedian(ifs_data['height'],0), '--', color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    # plt.plot(np.nanmedian(mask2,0),np.nanmedian(misc_data['height'],0), '--', color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    # plt.plot(np.nanmedian(mask4,0),np.nanmedian(ra2t_data['height'],0), '--', color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    # plt.plot(np.nanmedian(mask1,0),np.nanmedian(um_data['height'],0), '--', color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    plt.xlabel('TWC Cloud mask')
    plt.ylabel('Z [km]')
    plt.ylim([0,4000])
    plt.yticks(np.arange(0,9.01e3,0.5e3))
    ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    plt.xlim([0,1])
    plt.xticks(np.arange(0,1.01,0.1))
    ax1.set_xticklabels([0,' ',0.2,' ',0.4,' ',0.6,' ',0.8,' ',1.0])
    plt.legend()

    # plt.grid('on')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_RA2M_CASIM_RA2T_TWC-MASK_profiles.png'
    plt.savefig(fileout)

    # print ('Z1 = ')
    # print (np.round(np.nanmean(um_data['height'],0),-2))
    # print ('Z3 = ')
    # print (np.round(np.nanmean(ifs_data['height'],0),-2))
    #
    # Zindex1 = np.where(np.round(np.nanmean(um_data['height'],0),-2) < 600)#<= 2e3) #== 5.e+02)#
    # Zindex3 = np.where(np.round(np.nanmean(ifs_data['height'],0),-2) < 600)# <= 2e3) #== 5.e+02)#
    # print ('Zindex1 = ')
    # print (np.nanmean(um_data['height'][:,Zindex1[0]],0))
    # print ('Zindex3 = ')
    # print (np.nanmean(ifs_data['height'][:,Zindex3[0]],0))
    # print ('UM_RA2M = ')
    # print (np.nanmean(mask1[:,Zindex1[0]],0))
    # print ('UM_RA2T = ')
    # print (np.nanmean(mask4[:,Zindex1[0]],0))
    # print ('UM_CASIM-100 = ')
    # print (np.nanmean(mask2[:,Zindex1[0]],0))
    # print ('ECMWF_IFS = ')
    # print (np.nanmean(mask3[:,Zindex3[0]],0))
    # print ('Obs = ')
    # print (np.nanmean(mask0[:,Zindex1[0]],0))

    # print ('ECMWF_IFS/Obs = ')
    # print ((np.nanmean(mask3[:,Zindex3[0]],0) - np.nanmean(mask0[:,Zindex1[0]],0)) / np.nanmean(mask0[:,Zindex1[0]],0))


def plot_lwcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting LWC statistics for whole drift period:')
    print ('')

    # print (um_data.keys())

    ###----------------------------------------------------------------
    ###         1) LWC Method
    ###----------------------------------------------------------------
    # #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] == -999] = np.nan
    # obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    # obs_data['lwc'][obs_data['lwc'] >= 1e-3] = np.nan       ## exclude >1g/m3
    # obs_data['lwc'][obs_data['lwc'] < 1e-6] = np.nan       ## exclude <0.001g/m3
    # um_data['model_lwc'][um_data['model_lwc'] <= 0.0] = np.nan
    # um_data['model_lwc'][um_data['model_lwc'] < 1e-6] = np.nan
    # ifs_data['model_lwc'][ifs_data['model_lwc'] <= 0.0] = np.nan
    # ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 1e-6] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] <= 0] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 1e-6] = np.nan
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] <= 0] = np.nan
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 1e-6] = np.nan

    ###----------------------------------------------------------------
    ###         2) TWC Method
    ###----------------------------------------------------------------
    #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 5.0e-3] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] == -999] = 0.0
    # obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] == -999] = 0.0
    # um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                    obs_data['lwc'][t,k] = np.nan
                if obs_data['twc_ad'][t,k] < twc_thresh_um[k]:
                    obs_data['twc_ad'][t,k] = np.nan
                    obs_data['lwc_adiabatic'][t,k] = np.nan
                if obs_data['twc_ad_nolwp'][t,k] < twc_thresh_um[k]:
                    obs_data['twc_ad_nolwp'][t,k] = np.nan
                    obs_data['lwc_adiabatic_inc_nolwp'][t,k] = np.nan
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
                um_data['model_lwc'][t,k] = np.nan
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
                misc_data['model_lwc'][t,k] = np.nan
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
                ra2t_data['model_lwc'][t,k] = np.nan
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
                ifs_data['model_lwc'][t,k] = np.nan

    ###----------------------------------------------------------------
    ###         Plot figure - Mean profiles
    ###----------------------------------------------------------------

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()
    if obs_switch == 'RADAR':
        plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs_' + obs_switch + 'grid')
        ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
            np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.xlim([0,1.0])
    else:
        # #### SCALED LWC
        # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'grey', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-obsLWP', zorder = 5)
        # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
        #     np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.3)
        # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'grey', linewidth = 0.5)
        # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'grey', linewidth = 0.5)
        #### ADIABATIC LWC (where there are HATPRO LWP data available)
        plt.plot(np.nanmean(obs_data['lwc_adiabatic'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic', zorder = 5)
        # plt.plot(np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3,np.nanmedian(obs_data['height'],0), '--', color = 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic', zorder = 5)
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc_adiabatic'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
            np.nanmean(obs_data['lwc_adiabatic'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.xlim([0,0.2])
        plt.plot(np.nanmean(obs_data['lwc_adiabatic'],0)*1e3 - np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data['lwc_adiabatic'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        # #### ADIABATIC LWC (all times)
        # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic-incNoLWP', zorder = 5)
        # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
        #     np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.xlim([0,0.2])
        # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 - np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)
        # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)


    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3,
        np.nanmean(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_lwc'],0)*1e3 - np.nanstd(misc_data['model_lwc']*1e3,0),
        np.nanmean(misc_data['model_lwc'],0)*1e3 + np.nanstd(misc_data['model_lwc'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(misc_data['model_lwc'],0)*1e3 - np.nanstd(misc_data['model_lwc'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(misc_data['model_lwc'],0)*1e3 + np.nanstd(misc_data['model_lwc'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(ra2t_data['model_lwc'],0)*1e3 - np.nanstd(ra2t_data['model_lwc']*1e3,0),
        np.nanmean(ra2t_data['model_lwc'],0)*1e3 + np.nanstd(ra2t_data['model_lwc'],0)*1e3, color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmean(ra2t_data['model_lwc'],0)*1e3 - np.nanstd(ra2t_data['model_lwc'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(ra2t_data['model_lwc'],0)*1e3 + np.nanstd(ra2t_data['model_lwc'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc']*1e3,0),
        np.nanmean(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(ifs_data['model_lwc'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    plt.plot(np.nanmean(misc_data['model_lwc'],0)*1e3,np.nanmean(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    plt.plot(np.nanmean(ra2t_data['model_lwc'],0)*1e3,np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(um_data['model_lwc'],0)*1e3,np.nanmean(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    # plt.plot(np.nanmedian(ifs_data['model_lwc'],0)*1e3,np.nanmedian(ifs_data['height'],0), '--', color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    # plt.plot(np.nanmedian(misc_data['model_lwc'],0)*1e3,np.nanmedian(misc_data['height'],0), '--', color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    # plt.plot(np.nanmedian(ra2t_data['model_lwc'],0)*1e3,np.nanmedian(ra2t_data['height'],0), '--', color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    # plt.plot(np.nanmedian(um_data['model_lwc'],0)*1e3,np.nanmedian(um_data['height'],0), '--', color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    plt.xlabel('Liquid water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    # plt.ylim([0,5e3])
    # plt.yticks(np.arange(0,5.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
    plt.ylim([0,9000])
    plt.yticks(np.arange(0,9.01e3,0.5e3))
    ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    plt.xlim([0,0.2])
    plt.xticks(np.arange(0,0.21,0.025))
    ax1.set_xticklabels([0,' ',0.05,' ',0.1,' ',0.15,' ',0.2])
    plt.legend()

    # plt.grid('on')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid-V6_IFS_RA2M_CASIM-100_RA2T_LWC_MTThresholding-wLWCadiabatic-noOffsetLWP_226-257DOY_fixedRA2T_newColours_wSetFlags.png'
        # fileout = 'FIGS/Obs-' + obs_switch + 'grid-QF30_LWC_MTThresholding-wLWCadiabatic_noOffsetLWP_226-257DOY_blueNaNs_newColours.png'
    plt.savefig(fileout)
    plt.show()

    ###----------------------------------------------------------------
    ###         Plot figure - Median profiles
    ###----------------------------------------------------------------

    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.figure(figsize=(4.5,6))
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    #
    # ### define axis instance
    # ax1 = plt.gca()
    # if obs_switch == 'RADAR':
    #     plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs_' + obs_switch + 'grid')
    #     ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
    #         np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    #     plt.xlim([0,1.0])
    # else:
    #     # #### SCALED LWC
    #     # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'grey', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-obsLWP', zorder = 5)
    #     # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
    #     #     np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, color = 'lightgrey', alpha = 0.3)
    #     # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3, np.nanmean(obs_data['height'],0),
    #     #     '--', color = 'grey', linewidth = 0.5)
    #     # plt.plot(np.nanmean(obs_data['lwc'],0)*1e3 + np.nanstd(obs_data['lwc'],0)*1e3, np.nanmean(obs_data['height'],0),
    #     #     '--', color = 'grey', linewidth = 0.5)
    #     #### ADIABATIC LWC (where there are HATPRO LWP data available)
    #     # plt.plot(np.nanmean(obs_data['lwc_adiabatic'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic', zorder = 5)
    #     plt.plot(np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3,np.nanmedian(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic', zorder = 5)
    #     ax1.fill_betweenx(np.nanmedian(obs_data['height'],0),np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
    #         np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    #     # plt.xlim([0,0.2])
    #     plt.plot(np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3 - np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, np.nanmedian(obs_data['height'],0),
    #         '--', color = 'k', linewidth = 0.5)
    #     plt.plot(np.nanmedian(obs_data['lwc_adiabatic'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic'],0)*1e3, np.nanmedian(obs_data['height'],0),
    #         '--', color = 'k', linewidth = 0.5)
    #     # #### ADIABATIC LWC (all times)
    #     # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'purple', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adiabatic-incNoLWP', zorder = 5)
    #     # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 - np.nanstd(obs_data['lwc'],0)*1e3,
    #     #     np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    #     # plt.xlim([0,0.2])
    #     # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 - np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, np.nanmean(obs_data['height'],0),
    #     #     '--', color = 'k', linewidth = 0.5)
    #     # plt.plot(np.nanmean(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3 + np.nanstd(obs_data['lwc_adiabatic_inc_nolwp'],0)*1e3, np.nanmean(obs_data['height'],0),
    #     #     '--', color = 'k', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(ifs_data['height'],0),np.nanmedian(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3,
    #     np.nanmedian(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmedian(ifs_data['model_lwc'],0)*1e3 - np.nanstd(ifs_data['model_lwc'],0)*1e3, np.nanmedian(ifs_data['height'],0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmedian(ifs_data['model_lwc'],0)*1e3 + np.nanstd(ifs_data['model_lwc'],0)*1e3, np.nanmedian(ifs_data['height'],0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(misc_data['height'],0),np.nanmedian(misc_data['model_lwc'],0)*1e3 - np.nanstd(misc_data['model_lwc']*1e3,0),
    #     np.nanmedian(misc_data['model_lwc'],0)*1e3 + np.nanstd(misc_data['model_lwc'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(misc_data['model_lwc'],0)*1e3 - np.nanstd(misc_data['model_lwc'],0)*1e3, np.nanmedian(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(misc_data['model_lwc'],0)*1e3 + np.nanstd(misc_data['model_lwc'],0)*1e3, np.nanmedian(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(ra2t_data['height'],0),np.nanmedian(ra2t_data['model_lwc'],0)*1e3 - np.nanstd(ra2t_data['model_lwc']*1e3,0),
    #     np.nanmedian(ra2t_data['model_lwc'],0)*1e3 + np.nanstd(ra2t_data['model_lwc'],0)*1e3, color = 'lightblue', alpha = 0.3)
    # plt.plot(np.nanmedian(ra2t_data['model_lwc'],0)*1e3 - np.nanstd(ra2t_data['model_lwc'],0)*1e3, np.nanmedian(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(ra2t_data['model_lwc'],0)*1e3 + np.nanstd(ra2t_data['model_lwc'],0)*1e3, np.nanmedian(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(um_data['height'],0),np.nanmedian(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc']*1e3,0),
    #     np.nanmedian(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(um_data['model_lwc'],0)*1e3 - np.nanstd(um_data['model_lwc'],0)*1e3, np.nanmedian(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(um_data['model_lwc'],0)*1e3 + np.nanstd(um_data['model_lwc'],0)*1e3, np.nanmedian(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(ifs_data['model_lwc'],0)*1e3,np.nanmedian(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    # plt.plot(np.nanmedian(misc_data['model_lwc'],0)*1e3,np.nanmedian(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    # plt.plot(np.nanmedian(ra2t_data['model_lwc'],0)*1e3,np.nanmedian(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    # plt.plot(np.nanmedian(um_data['model_lwc'],0)*1e3,np.nanmedian(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)
    #
    # plt.xlabel('Liquid water content [g m$^{-3}$]')
    # plt.ylabel('Z [km]')
    # # plt.ylim([0,5e3])
    # # plt.yticks(np.arange(0,5.01e3,0.5e3))
    # # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
    # plt.ylim([0,9000])
    # plt.yticks(np.arange(0,9.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    # plt.xlim([0,0.2])
    # plt.xticks(np.arange(0,0.21,0.025))
    # ax1.set_xticklabels([0,' ',0.05,' ',0.1,' ',0.15,' ',0.2])
    # plt.legend()
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # if month_flag == -1:
    #     fileout = 'FIGS/Obs-' + obs_switch + 'grid-QF30_IFS_RA2M_CASIM-100_RA2T_LWC_MTThresholding-wLWCadiabatic-noOffsetLWP_226-257DOY_fixedRA2T_newColours_wSetFlags_Medians.svg'
    #     # fileout = 'FIGS/Obs-' + obs_switch + 'grid-QF30_LWC_MTThresholding-wLWCadiabatic_noOffsetLWP_226-257DOY_blueNaNs_newColours.png'
    # # plt.savefig(fileout)
    # plt.show()

    print ('Z1 = ')
    print (np.round(np.nanmean(um_data['height'],0),-2))
    print ('Z3 = ')
    print (np.round(np.nanmean(ifs_data['height'],0),-2))

    Zindex1 = np.where(np.round(np.nanmean(um_data['height'],0),-2) == 5.e+02)# == 4600)#<= 2e3) #
    Zindex3 = np.where(np.round(np.nanmean(ifs_data['height'],0),-2) == 5.e+02)#== 4400)# <= 2e3) #
    print ('Zindex1 = ')
    print (np.nanmean(um_data['height'][:,Zindex1[0]],0))
    print ('Zindex3 = ')
    print (np.nanmean(ifs_data['height'][:,Zindex3[0]],0))
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_lwc'][:,Zindex1[0]],0)*1e3)
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_lwc'][:,Zindex1[0]],0)*1e3)
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_lwc'][:,Zindex1[0]],0)*1e3)
    print ('ECMWF_IFS = ')
    print (np.nanmean(ifs_data['model_lwc'][:,Zindex3[0]],0)*1e3)
    print ('Obs = ')
    print (np.nanmean(obs_data['lwc_adiabatic'][:,Zindex1[0]],0)*1e3)

def plot_iwcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting IWC statistics for whole drift period:')
    print ('')

    # print (um_data.keys())

    ###----------------------------------------------------------------
    ###         1) IWC Method
    ###----------------------------------------------------------------
    # #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] == -999] = np.nan
    # obs_data['iwc'][obs_data['iwc'] < 1e-6] = np.nan
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 1e-6] = np.nan
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 1e-6] = np.nan
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 1e-6] = np.nan
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 1e-6] = np.nan

    ###----------------------------------------------------------------
    ###         2) TWC Method
    ###----------------------------------------------------------------
    #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 5.0e-3] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] == -999] = 0.0
    # # obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] == -999] = 0.0
    # um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                    obs_data['iwc'][t,k] = np.nan
                    obs_data['lwc_adiabatic'][t,k] = np.nan
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
                um_data['model_iwc_filtered'][t,k] = np.nan
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
                misc_data['model_iwc_filtered'][t,k] = np.nan
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
                ra2t_data['model_iwc_filtered'][t,k] = np.nan
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
                ifs_data['model_snow_iwc_filtered'][t,k] = np.nan

    ##################################################
    ##################################################
    #### 	Plot mean profiles
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
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    if obs_switch == 'RADAR':
        plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs-' + obs_switch + 'grid')
        ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
            np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    else:
        plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 5)
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
            np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3,
        np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, color = 'navajowhite', alpha = 0.4)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5, zorder = 3)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 - np.nanstd(misc_data['model_iwc_filtered']*1e3,0),
        np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 + np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 - np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered'],0)*1e3 + np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3 - np.nanstd(ra2t_data['model_iwc_filtered']*1e3,0),
        np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3 + np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, color = 'lightblue', alpha = 0.3)
    plt.plot(np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3 - np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3 + np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered']*1e3,0),
        np.nanmean(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered'],0)*1e3,np.nanmean(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    plt.plot(np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3,np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(um_data['model_iwc_filtered'],0)*1e3,np.nanmean(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    plt.xlabel('Ice water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(np.arange(0,9.01e3,0.5e3))
    ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    plt.xlim([0,0.06])
    plt.xticks(np.arange(0,0.061,0.0075))
    ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
    plt.legend()

    # plt.grid('on')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid-V6_UM_IFS_CASIM-100_IWC-MTThresholding-wLWCadiabatic-noOfsetLWP_226-257DOY_fixedRA2T_newColours_wSetFlags.png'
    plt.savefig(fileout)
    plt.show()

    ##################################################
    ##################################################
    #### 	Plot median profiles
    ##################################################
    ##################################################

    # SMALL_SIZE = 12
    # MED_SIZE = 14
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.rc('legend',fontsize=LARGE_SIZE)
    # plt.figure(figsize=(4.5,6))
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    #
    # ### define axis instance
    # ax1 = plt.gca()
    #
    # plt.plot(np.nanmedian(obs_data['iwc'],0)*1e3,np.nanmedian(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 5)
    # ax1.fill_betweenx(np.nanmedian(obs_data['height'],0),np.nanmedian(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
    #     np.nanmedian(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmedian(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3, np.nanmedian(obs_data['height'],0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmedian(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, np.nanmedian(obs_data['height'],0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(ifs_data['height'],0),np.nanmedian(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3,
    #     np.nanmedian(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, color = 'navajowhite', alpha = 0.4)
    # plt.plot(np.nanmedian(ifs_data['model_snow_iwc_filtered'],0)*1e3 - np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, np.nanmedian(ifs_data['height'],0),
    #     '--', color = 'gold', linewidth = 0.5, zorder = 3)
    # plt.plot(np.nanmedian(ifs_data['model_snow_iwc_filtered'],0)*1e3 + np.nanstd(ifs_data['model_snow_iwc_filtered'],0)*1e3, np.nanmedian(ifs_data['height'],0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(misc_data['height'],0),np.nanmedian(misc_data['model_iwc_filtered'],0)*1e3 - np.nanstd(misc_data['model_iwc_filtered']*1e3,0),
    #     np.nanmedian(misc_data['model_iwc_filtered'],0)*1e3 + np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmedian(misc_data['model_iwc_filtered'],0)*1e3 - np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, np.nanmedian(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmedian(misc_data['model_iwc_filtered'],0)*1e3 + np.nanstd(misc_data['model_iwc_filtered'],0)*1e3, np.nanmedian(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(ra2t_data['height'],0),np.nanmedian(ra2t_data['model_iwc_filtered'],0)*1e3 - np.nanstd(ra2t_data['model_iwc_filtered']*1e3,0),
    #     np.nanmedian(ra2t_data['model_iwc_filtered'],0)*1e3 + np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, color = 'lightblue', alpha = 0.3)
    # plt.plot(np.nanmedian(ra2t_data['model_iwc_filtered'],0)*1e3 - np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, np.nanmedian(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(ra2t_data['model_iwc_filtered'],0)*1e3 + np.nanstd(ra2t_data['model_iwc_filtered'],0)*1e3, np.nanmedian(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax1.fill_betweenx(np.nanmedian(um_data['height'],0),np.nanmedian(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered']*1e3,0),
    #     np.nanmedian(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmedian(um_data['model_iwc_filtered'],0)*1e3 - np.nanstd(um_data['model_iwc_filtered'],0)*1e3, np.nanmedian(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmedian(um_data['model_iwc_filtered'],0)*1e3 + np.nanstd(um_data['model_iwc_filtered'],0)*1e3, np.nanmedian(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmedian(ifs_data['model_snow_iwc_filtered'],0)*1e3,np.nanmedian(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    # plt.plot(np.nanmedian(misc_data['model_iwc_filtered'],0)*1e3,np.nanmedian(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    # plt.plot(np.nanmedian(ra2t_data['model_iwc_filtered'],0)*1e3,np.nanmedian(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    # plt.plot(np.nanmedian(um_data['model_iwc_filtered'],0)*1e3,np.nanmedian(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)
    #
    # plt.xlabel('Ice water content [g m$^{-3}$]')
    # plt.ylabel('Z [km]')
    # plt.ylim([0,9000])
    # plt.yticks(np.arange(0,9.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    # plt.xlim([0,0.06])
    # plt.xticks(np.arange(0,0.061,0.0075))
    # ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
    # plt.legend()
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # if month_flag == -1:
    #     fileout = 'FIGS/Obs-' + obs_switch + 'grid_UM_IFS_CASIM-100_IWC-MTThresholding-wLWCadiabatic-noOfsetLWP_226-257DOY_fixedRA2T_newColours_wSetFlags_Median.svg'
    # plt.savefig(fileout)
    # plt.show()

    # print ('Z1 = ')
    # print (np.round(np.nanmean(um_data['height'],0),-2))
    # print ('Z3 = ')
    # print (np.round(np.nanmean(ifs_data['height'],0),-2))

    Zindex1 = np.where(np.round(np.nanmean(um_data['height'],0),-2) <= 4000)#<= 2e3) # == 5.e+02)#
    Zindex3 = np.where(np.round(np.nanmean(ifs_data['height'],0),-2) <= 4000)# <= 2e3) #== 5.e+02)#
    print ('Zindex1 = ')
    print (np.nanmean(um_data['height'][:,Zindex1[0]],0))
    print ('Zindex3 = ')
    print (np.nanmean(ifs_data['height'][:,Zindex3[0]],0))
    # print ('UM_RA2M = ')
    # print (np.nanmean(um_data['model_iwc_filtered'][:,Zindex1[0]],0)*1e3)
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_iwc_filtered'],0)*1e3)
    # print ('UM_CASIM-100 = ')
    # print (np.nanmean(misc_data['model_iwc_filtered'][:,Zindex1[0]],0)*1e3)
    # print ('ECMWF_IFS = ')
    # print (np.nanmean(ifs_data['model_snow_iwc_filtered'][:,Zindex3[0]],0)*1e3)
    print ('Obs = ')
    print (np.nanmean(obs_data['iwc'],0)*1e3)

    print ('UM_RA2T - Obs = ')
    print ((np.nanmean(ra2t_data['model_iwc_filtered'][:,Zindex1[0]],0) - np.nanmean(obs_data['iwc'][:,Zindex1[0]],0))*1e3)

def plot_twcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    import matplotlib.colors as colors
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting TWC profiles for whole drift period:')
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
    plt.figure(figsize=(4.5,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    ax1 = plt.gca()

    print (um_data.keys())

    ###----------------------------------------------------------------
    ###         2) TWC Method
    ###----------------------------------------------------------------
    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 5.0e-3] = np.nan
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
    ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] == -999] = 0.0
    # obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] == -999] = 0.0
    um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
    ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                    obs_data['iwc'][t,k] = np.nan
                    obs_data['lwc_adiabatic'][t,k] = np.nan
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
                um_data['model_iwc_filtered'][t,k] = np.nan
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
                misc_data['model_iwc_filtered'][t,k] = np.nan
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
                ra2t_data['model_iwc_filtered'][t,k] = np.nan
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
                ifs_data['model_snow_iwc_filtered'][t,k] = np.nan

    # viridis = mpl_cm.get_cmap('viridis', 256)
    # newcolors = viridis(np.linspace(0, 1, 256))
    # greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    # newcolors[:20, :] = greyclr
    # newcmp = ListedColormap(newcolors)

    twc0 = np.transpose(obs_data['twc'])*1e3

    if obs_switch == 'RADAR':
        plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,obs_data['height'][:394], 'k--', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 3)
        ax1.fill_betweenx(obs_data['height'][:394],np.nanmean(obs_data['twc'],0)*1e3 - np.nanstd(obs_data['twc'],0)*1e3,
            np.nanmean(obs_data['twc'],0)*1e3 + np.nanstd(obs_data['twc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    else:
        plt.plot(np.nanmean(obs_data['twc'],0)*1e3,np.nanmean(obs_data['height'],0), 'k', linewidth = 3, label = 'Obs-' + obs_switch + 'grid', zorder = 5)
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['twc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
            np.nanmean(obs_data['twc'],0)*1e3 + np.nanstd(obs_data['twc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
        plt.plot(np.nanmean(obs_data['twc'],0)*1e3 - np.nanstd(obs_data['twc'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data['twc'],0)*1e3 + np.nanstd(obs_data['twc'],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['model_twc'],0)*1e3 - np.nanstd(ifs_data['model_twc'],0)*1e3,
        np.nanmean(ifs_data['model_twc'],0)*1e3 + np.nanstd(ifs_data['model_twc'],0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(ifs_data['model_twc'],0)*1e3 - np.nanstd(ifs_data['model_twc'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(ifs_data['model_twc'],0)*1e3 + np.nanstd(ifs_data['model_twc'],0)*1e3, np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(misc_data['height'],0),np.nanmean(misc_data['model_twc'],0)*1e3 - np.nanstd(misc_data['model_twc'],0)*1e3,
        np.nanmean(misc_data['model_twc'],0)*1e3 + np.nanstd(misc_data['model_twc'],0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(misc_data['model_twc'],0)*1e3 - np.nanstd(misc_data['model_twc'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(misc_data['model_twc'],0)*1e3 + np.nanstd(misc_data['model_twc'],0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(ra2t_data['model_twc'],0)*1e3 - np.nanstd(ra2t_data['model_twc'],0)*1e3,
        np.nanmean(ra2t_data['model_twc'],0)*1e3 + np.nanstd(ra2t_data['model_twc'],0)*1e3, color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(ra2t_data['model_twc'],0)*1e3 - np.nanstd(ra2t_data['model_twc'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(ra2t_data['model_twc'],0)*1e3 + np.nanstd(ra2t_data['model_twc'],0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax1.fill_betweenx(np.nanmean(um_data['height'],0),np.nanmean(um_data['model_twc'],0)*1e3 - np.nanstd(um_data['model_twc'],0)*1e3,
        np.nanmean(um_data['model_twc'],0)*1e3 + np.nanstd(um_data['model_twc'],0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(um_data['model_twc'],0)*1e3 - np.nanstd(um_data['model_twc'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(um_data['model_twc'],0)*1e3 + np.nanstd(um_data['model_twc'],0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(ifs_data['model_twc'],0)*1e3,np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'ECMWF_IFS', zorder = 4)
    plt.plot(np.nanmean(misc_data['model_twc'],0)*1e3,np.nanmean(misc_data['height'],0), color = 'mediumseagreen', linewidth = 3, label = 'UM_CASIM-100', zorder = 3)
    plt.plot(np.nanmean(ra2t_data['model_twc'],0)*1e3,np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(um_data['model_twc'],0)*1e3,np.nanmean(um_data['height'],0), color = 'darkblue', linewidth = 3, label = 'UM_RA2M', zorder = 1)

    plt.xlabel('Total water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(np.arange(0,9.01e3,0.5e3))
    ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5,' ',6,' ',7,' ',8,' ',9])
    plt.xlim([0,0.25])
    plt.xticks(np.arange(0,0.251,0.025))
    ax1.set_xticklabels([0,' ',0.05,' ',0.1,' ',0.15,' ',0.2,' ',0.25])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid-QF30_IFS_CASIM-100_RA2M_RA2T_TWC_226-257DOY_TWC-MTThresholding_newColours.svg'
    # plt.savefig(fileout)
    plt.show()


def plot_TWCTesting(um_data, ifs_data, misc_data, obs_data, data1, data2, data3, obs, month_flag, missing_files, doy):

    #### set flagged um_data to nans
    obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    um_data['model_iwc'][um_data['model_iwc'] < 0.0] = 0.0
    # ifs_data['model_iwc_filtered'][ifs_data['model_iwc_filtered'] <= 0.0] = np.nan
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 1.0] = 0.0
    misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
    um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = 0.0
    misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content - Cloudnet
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']

    ###----------------------------------------------------------------
    ###         Calculate total water content - Model
    ###----------------------------------------------------------------
    data1['qliq'][data1['qliq'] < 0.0] = 0.0
    data1['qice'][data1['qice'] < 0.0] = 0.0
    data1['qtot'] = data1['qliq'] + data1['qice']
    data2['qliq'][data2['qliq'] < 0.0] = 0.0
    data2['qice'][data2['qice'] < 0.0] = 0.0
    data2['qtot'] = data2['qliq'] + data2['qice']
    data3['ql'][data3['ql'] < 0.0] = 0.0
    data3['qi'][data3['qi'] < 0.0] = 0.0
    data3['qtot'] = data3['ql'] + data3['qi']


    ###----------------------------------------------------------------
    ###         New figure - timeseries
    ###----------------------------------------------------------------

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 1.0, left = 0.1,
            hspace = 0.4, wspace = 0.05)

    ### define axis instance
    ax = plt.gca()

    plt.subplot(211)
    plt.pcolormesh(um_data['time'], np.squeeze(um_data['height'][0,:]), np.transpose(um_data['model_twc'])*1e3,
        # )
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlim([226,257])
    plt.title('TWC [g/m3]; Cloudnet')
    plt.colorbar()

    plt.subplot(212)
    plt.pcolormesh(data1['time'], data1['height'], np.transpose(data1['qtot'])*1e3,
        # )
        vmin = 0.0, vmax = 0.5)
        #cmap = newcmp)
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.xlim([226,257])
    plt.title('TWC [g/kg]; Model')
    plt.colorbar()

    plt.show()

    ###----------------------------------------------------------------
    ###         New figure - profiles
    ###----------------------------------------------------------------
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.98, left = 0.1,
            hspace = 0.45, wspace = 0.25)
    ### define axis instance
    ax = plt.gca()

    plt.subplot(331)
    plt.plot(np.nanmean(um_data['model_lwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data1['qliq']*1e3,0), data1['height'], label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: LWC')

    plt.subplot(332)
    plt.plot(np.nanmean(um_data['model_iwc_filtered']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data1['qice']*1e3,0), data1['height'], label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: IWC')

    plt.subplot(333)
    plt.plot(np.nanmean(um_data['model_twc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data1['qtot']*1e3,0), data1['height'], label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_RA2M: TWC')

    plt.subplot(334)
    plt.plot(np.nanmean(misc_data['model_lwc']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data2['qliq']*1e3,0), data2['height'], label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: LWC')

    plt.subplot(335)
    plt.plot(np.nanmean(misc_data['model_iwc_filtered']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data2['qice']*1e3,0), data2['height'], label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: IWC')

    plt.subplot(336)
    plt.plot(np.nanmean(misc_data['model_twc']*1e3,0), np.squeeze(misc_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data2['qtot']*1e3,0), data2['height'], label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    # plt.xlim([226,257])
    plt.title('UM_CASIM-100: TWC')

    plt.subplot(337)
    plt.plot(np.nanmean(ifs_data['model_lwc']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data3['ql']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    # plt.xlim([226,257])
    plt.title('ECMWF_IFS: LWC')

    plt.subplot(338)
    plt.plot(np.nanmean(ifs_data['model_snow_iwc_filtered']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, filtered, g/m3')
    plt.plot(np.nanmean(data3['qi']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    # plt.plot(np.nanmean(um_data['model_iwc']*1e3,0), np.squeeze(um_data['height'][0,:]), label = 'Cloudnet, not filtered, g/m3')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    # plt.legend()
    plt.xlim([0.0,0.011])
    plt.title('ECMWF: IWC')

    plt.subplot(339)
    plt.plot(np.nanmean(ifs_data['model_twc']*1e3,0), np.squeeze(ifs_data['height'][0,:]), label = 'Cloudnet, g/m3')
    plt.plot(np.nanmean(data3['qtot']*1e3,0), np.nanmean(data3['height'],0), label = 'Model, g/kg')
    # plt.ylabel('Height [m]')
    plt.ylim([0,9000])
    plt.legend()
    plt.xlim([0.0,0.12])
    plt.title('ECMWF_IFS: TWC')

    plt.savefig('../FIGS/comparisons/CN-ModelComparison_LWC-IWC-TWC_profiles.png')
    plt.show()

def plot_LWP(um_data, ifs_data, misc_data, ra2t_data, obs_data, obs, month_flag, missing_files, um_out_dir, doy, obs_switch):#,lwp): #, lon, lat):

    ###################################
    ## PLOT TIMESERIES
    ###################################

    print ('******')
    print ('')
    print ('Plotting LWP timeseries for whole drift period:')
    print ('')

    ##################################################
    ##################################################
    #### 	SET UP FIGURE PROPERTIES
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
    fig = plt.figure(figsize=(12,3.5))
    plt.subplots_adjust(top = 0.9, bottom = 0.14, right = 0.96, left = 0.1,
            hspace = 0.4, wspace = 0.1)

    ### define axis instance
    # ax = plt.gca()

    print (um_data.keys())
    print (obs_data.keys())

    #### set flagged and bad data to nans
    um_data['model_lwp'][um_data['model_lwp'] <= 0] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] <= 0] = np.nan
    misc_data['model_lwp'][misc_data['model_lwp'] <= 0] = np.nan
    ra2t_data['model_lwp'][ra2t_data['model_lwp'] <= 0] = np.nan
    # um_data['model_lwp'][um_data['model_lwp'] >= 1000] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] >= 1.0] = np.nan
    # misc_data['model_lwp'][misc_data['model_lwp'] >= 1000] = np.nan
    if obs_switch == 'RADAR':
        obs_data['lwp'][obs_data['lwp'] < 0] = np.nan     ### index 0 is mean
        # obs_data['lwp'][obs_data['lwp'] > 1.2] = np.nan    ### >1.2 == >1200g/m2
    else:
        obs_data['lwp'][obs_data['lwp'][:,0] == -999.0, 0] = np.nan     ### index 0 is mean
        # lwp[lwp[:,0] == -999.0, 0] = np.nan
        # obs_data['lwp'][obs_data['lwp'][:,0] > 1.2, 0] = np.nan    ### >1.2 == >1200g/m2
        # obs_data['lwp'][obs_data['lwp'][:,1] < 0, 1] = np.nan     ### index 1 is min
        # obs_data['lwp'][obs_data['lwp'][:,1] > 1.2, 1] = np.nan    ###>1.2 == >1200g/m2
        # obs_data['lwp'][obs_data['lwp'][:,2] < 0, 2] = np.nan     ### index 2 is max
        # obs_data['lwp'][obs_data['lwp'][:,2] > 1.2, 2] = np.nan    ### >1.2 == >1200g/m2

    res = 3 ### hourly resolution to plot

    ax  = fig.add_axes([0.09,0.18,0.6,0.76])   # left, bottom, width, height
    plt.plot(obs['deck7th']['doy'][:],obs['deck7th']['lwp'][:], color = 'grey', label = 'Obs_HATPRO', zorder = 2)
    if obs_switch == 'RADAR':
        plt.plot(obs_data['time'][:],obs_data['lwp'][:]*1e3, color = 'purple', label = 'Obs_' + obs_switch + 'grid')
    else:
        plt.plot(obs_data['time'][:],obs_data['lwp'][:,0]*1e3, color = 'black', label = 'Obs_' + obs_switch + 'grid', zorder = 3)
        # plt.plot(obs_data['time'][:],lwp[:,0]*1e3, color = 'purple', label = 'Obs_' + obs_switch + 'grid2', zorder = 3)
        # plt.plot(obs_data['time'][:],obs_data['lwp'][:,1]*1e3, 'k--')
        # plt.plot(obs_data['time'][:],obs_data['lwp'][:,2]*1e3, 'k--')
        # ax.fill_between(obs_data['time'][:], obs_data['lwp'][:,1]*1e3, obs_data['lwp'][:,2]*1e3, color = 'grey', alpha = 0.2)
    plt.plot(ifs_data['time'][::res],ifs_data['model_lwp'][::res]*1e3,
        'v', color = 'gold', markeredgecolor = 'orange', label = 'ECMWF_IFS')
    plt.plot(misc_data['time'][::res],misc_data['model_lwp'][::res]*1e3,
        '<', color = 'mediumseagreen', markeredgecolor = 'darkgreen', label = 'UM_CASIM-AeroProf')
    plt.plot(ra2t_data['time'][::res],ra2t_data['model_lwp'][::res]*1e3,
        '>', color = 'steelblue', markeredgecolor = 'darkslategrey', label = 'UM_RA2T')
    plt.plot(um_data['time'][::res],um_data['model_lwp'][::res]*1e3,
        '^', color = 'darkblue', markeredgecolor = 'midnightblue', label = 'UM_RA2M')
    plt.ylabel('LWP [g m$^{2}$]')
    plt.ylim([0,800])
    # plt.xlim([245,250])
    # plt.xlabel('DOY')
    plt.xlim([doy[0],doy[-1]])
    plt.xticks([230,235,240,245,250,255])
    ax.set_xticklabels(['18/8','23/8','28/8','2/9','7/9','12/9'])
    plt.xlabel('Date')
    plt.legend(bbox_to_anchor=(0.0, 0.7, 1., .102), loc=3, ncol=3)
    ax = plt.gca()
    nans = ax.get_ylim()
    for file in missing_files:
        ax.fill_between(np.arange(file, file + 1, 1/24.0), nans[0], nans[-1],
            facecolor = 'white',
            # hatch = 'x',
            zorder = 2)

    ax  = fig.add_axes([0.76,0.27,0.22,0.6])   # left, bottom, width, height
    yEmax = 0.012
    plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    sns.distplot(um_data['model_lwp']*1e3, hist=False, color="darkblue", kde_kws={"shade": True})
    sns.distplot(ra2t_data['model_lwp']*1e3, hist=False, color="steelblue", kde_kws={"shade": True})
    sns.distplot(misc_data['model_lwp']*1e3, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    sns.distplot(ifs_data['model_lwp']*1e3, hist=False, color="gold", kde_kws={"shade": True})
    sns.distplot(obs_data['lwp'][:,0]*1e3, hist=False, color="black")
    sns.distplot(obs['deck7th']['lwp'][:], hist=False, color="grey")
    plt.xlim([-50,500])
    plt.ylim([0,yEmax])
    plt.xlabel('LWP [g m$^{-2}$]')

    # print (obs['deck7th']['doy'][1500:2000])
    # print (obs_data['time'])
    # print (ifs_data['time'])

    # print (obs['deck7th']['doy'][np.where(np.logical_and(obs['deck7th']['doy'] >= 250.0, obs['deck7th']['doy'] <= 257.0))])
    # print (obs_data['time'][np.where(np.logical_and(obs_data['time'] >= 234, obs_data['time'] <= 236))])
    # print (obs_data['lwp'][np.where(np.logical_and(obs_data['time'] >= 234, obs_data['time'] <= 236)),0]*1e3)

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs-' + obs_switch + 'grid-qf30_IFS_RA2M_CASIM-AeroProf_RA2T_LWP_226-257DOY_newColours_Date_noOffsetLWP-LWPbugfixed_fixedRA2T.svg'
    # plt.savefig(fileout)
    plt.show()

    print (np.nanmax(um_data['model_lwp']*1e3))
    print (np.nanmax(obs_data['lwp'][:,0]*1e3))

    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    ####    DEFINE PERIODS
    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    p3 = np.where(np.logical_and(obs_data['time'] >= doy[0], obs_data['time'] < 230.0))
    p4 = np.where(np.logical_and(obs_data['time'] >= 230.0, obs_data['time'] < 240.0))
    p5 = np.where(np.logical_and(obs_data['time'] >= 240.0, obs_data['time'] < 247.0))
    p6 = np.where(np.logical_and(obs_data['time'] >= 247.0, obs_data['time'] < 251.0))

    print ('')
    print ('*****')
    print ('P3 = ')
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_lwp'][p3]*1e3))
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_lwp'][p3]*1e3))
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_lwp'][p3]*1e3))
    print ('ECMWF_IFS = ')
    print (np.nanmean(ifs_data['model_lwp'][p3]*1e3))
    print ('Obs = ')
    print (np.nanmean(obs_data['lwp'][p3,0]*1e3))

    print ('*****')
    print ('P4 = ')
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_lwp'][p4]*1e3))
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_lwp'][p4]*1e3))
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_lwp'][p4]*1e3))
    print ('ECMWF_IFS = ')
    print (np.nanmean(ifs_data['model_lwp'][p4]*1e3))
    print ('Obs = ')
    print (np.nanmean(obs_data['lwp'][p4,0]*1e3))

    print ('*****')
    print ('P5 = ')
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_lwp'][p5]*1e3))
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_lwp'][p5]*1e3))
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_lwp'][p5]*1e3))
    print ('ECMWF_IFS = ')
    print (np.nanmean(ifs_data['model_lwp'][p5]*1e3))
    print ('Obs = ')
    print (np.nanmean(obs_data['lwp'][p5,0]*1e3))

    print ('*****')
    print ('P6 = ')
    print ('UM_RA2M = ')
    print (np.nanmean(um_data['model_lwp'][p6]*1e3))
    print ('UM_RA2T = ')
    print (np.nanmean(ra2t_data['model_lwp'][p6]*1e3))
    print ('UM_CASIM-100 = ')
    print (np.nanmean(misc_data['model_lwp'][p6]*1e3))
    print ('ECMWF_IFS = ')
    print (np.nanmean(ifs_data['model_lwp'][p6]*1e3))
    print ('Obs = ')
    print (np.nanmean(obs_data['lwp'][p6,0]*1e3))

def plot_BiVAR(um_data, ifs_data, misc_data, ra2t_data, obs_data, obs, month_flag, missing_files, cn_um_out_dir, doy, obs_switch, data1, data2, data3, data4, nanind, wcind):

    from sklearn.metrics import r2_score
    from sklearn import linear_model
    from scipy import stats

    '''
        Following Jonny's comment:
        The one thing I was thinking that would be good to see in addition is a plot that
        really nails down the link between the errors in cloud cover and the incoming
        radiation terms. I think plotting bivariate histograms of Cv and LWP vs SWdown/LWdown/Rnet
        could really highlight the importance of these missing clear sky periods for the energy
        balance, and hint at knock on issues for the sea ice evolution at longer timescales.
        Since this is a common issue it may have relevance for cliamte projections as well.
    '''

    ###################################
    ## PLOT TIMESERIES
    ###################################

    print ('******')
    print ('')
    print ('Plotting LWP timeseries for whole drift period:')
    print ('')

    #### set flagged and bad data to nans
    um_data['model_lwp'][um_data['model_lwp'] <= 0] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] <= 0] = np.nan
    misc_data['model_lwp'][misc_data['model_lwp'] <= 0] = np.nan
    ra2t_data['model_lwp'][ra2t_data['model_lwp'] <= 0] = np.nan
    ifs_data['model_lwp'][ifs_data['model_lwp'] >= 1.0] = np.nan
    if obs_switch == 'RADAR':
        obs_data['lwp'][obs_data['lwp'] < 0] = np.nan     ### index 0 is mean
    else:
        obs_data['lwp'][obs_data['lwp'][:,0] == -999.0, 0] = np.nan     ### index 0 is mean

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    ####            CREATE TWC MASK
    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------

    mask0 = np.zeros([np.size(obs_data['twc'],0), np.size(obs_data['twc'],1)])
    mask1 = np.zeros([np.size(um_data['model_twc'],0), np.size(um_data['model_twc'],1)])
    mask2 = np.zeros([np.size(misc_data['model_twc'],0), np.size(misc_data['model_twc'],1)])
    mask3 = np.zeros([np.size(ifs_data['model_twc'],0), np.size(ifs_data['model_twc'],1)])
    mask4 = np.zeros([np.size(ra2t_data['model_twc'],0), np.size(ra2t_data['model_twc'],1)])

    # print ('mask0 size = ' + str(mask0.shape))

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                elif obs_data['twc'][t,k] >= twc_thresh_um[k]:
                    mask0[t,k] = 1.0
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
            elif um_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask1[t,k] = 1.0
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
            elif misc_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask2[t,k] = 1.0
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
            elif ra2t_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask4[t,k] = 1.0
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
            elif ifs_data['model_twc'][t,k] >= twc_thresh_ifs[k]:
                mask3[t,k] = 1.0

    mask0[nanind] = np.nan
    mask1[nanind] = np.nan
    mask2[nanind] = np.nan
    mask3[nanind] = np.nan
    mask4[nanind] = np.nan

    mask0[wcind] = np.nan
    mask1[wcind] = np.nan
    mask2[wcind] = np.nan
    mask3[wcind] = np.nan
    mask4[wcind] = np.nan

    print (mask0.shape)
    print (obs_data['time'][0])
    print (obs_data['time'][-4])
    print (data1['fixed_radiation']['time'][0])
    print (data1['fixed_radiation']['time'][-1])
    print (obs_data['lwp'][:-3,0].shape)

    ### index drift period in ship radiation measurements
    drift_ship = np.where(np.logical_and(obs['fixed_radiation']['time_ship'] >= data1['fixed_radiation']['time'][0],
                    obs['fixed_radiation']['time_ship'] < data1['fixed_radiation']['time'][-1] + 0.05))
                                ### obs time array ever so slightly offset (+) from models

    print (drift_ship[0].shape)

    print (data1['fixed_radiation']['time'][-1])
    print (obs['fixed_radiation']['time_ship'][drift_ship[0][-1]])

    ### --------------------------------------
    ### compute linear regression
    ### --------------------------------------
    ### define dictionaries
    data1['biases'] = {}
    data2['biases'] = {}
    data3['biases'] = {}
    data4['biases'] = {}

    bias_list = ['SWd', 'LWd', 'SWnet', 'LWnet','Rnet']
    for bias in bias_list:
        data1['biases'][bias + '-LWP_regres'] = {}
        data1['biases']['LWP'] = um_data['model_lwp'][:-3]*1e3 - obs_data['lwp'][:-3,0]*1e3
        data1['biases'][bias] = data1['fixed_radiation'][bias] - obs['fixed_radiation'][bias + '_ship'][drift_ship[0]]
        mask = ~np.isnan(data1['biases']['LWP']) & ~np.isnan(data1['biases'][bias])
        data1['biases'][bias + '-LWP_regres']['slope'], data1['biases'][bias + '-LWP_regres']['intercept'], data1['biases'][bias + '-LWP_regres']['r_value'], data1['biases'][bias + '-LWP_regres']['p_value'], data1['biases'][bias + '-LWP_regres']['std_err'] = stats.linregress(data1['biases']['LWP'][mask],data1['biases'][bias][mask])
        data1['biases'][bias + '-LWP_regres']['line'] = data1['biases'][bias + '-LWP_regres']['slope'] * data1['biases']['LWP'] + data1['biases'][bias + '-LWP_regres']['intercept']
        print('r-' + bias + '-um_ra2m:', data1['biases'][bias + '-LWP_regres']['r_value'])
        print('r-squared-' + bias + '-um_ra2m:', data1['biases'][bias + '-LWP_regres']['r_value']**2)

        data2['biases'][bias + '-LWP_regres'] = {}
        data2['biases']['LWP'] = misc_data['model_lwp'][:-3]*1e3 - obs_data['lwp'][:-3,0]*1e3
        data2['biases'][bias] = data2['fixed_radiation'][bias] - obs['fixed_radiation'][bias + '_ship'][drift_ship[0]]
        mask = ~np.isnan(data2['biases']['LWP']) & ~np.isnan(data2['biases'][bias])
        data2['biases'][bias + '-LWP_regres']['slope'], data2['biases'][bias + '-LWP_regres']['intercept'], data2['biases'][bias + '-LWP_regres']['r_value'], data2['biases'][bias + '-LWP_regres']['p_value'], data2['biases'][bias + '-LWP_regres']['std_err'] = stats.linregress(data2['biases']['LWP'][mask],data2['biases'][bias][mask])
        data2['biases'][bias + '-LWP_regres']['line'] = data2['biases'][bias + '-LWP_regres']['slope'] * data2['biases']['LWP'] + data2['biases'][bias + '-LWP_regres']['intercept']
        print('r-' + bias + '-um_casim-100:', data2['biases'][bias + '-LWP_regres']['r_value'])
        print('r-squared-' + bias + '-um_casim-100:', data2['biases'][bias + '-LWP_regres']['r_value']**2)

        data3['biases'][bias + '-LWP_regres'] = {}
        data3['biases']['LWP'] = ifs_data['model_lwp'][:-3]*1e3 - obs_data['lwp'][:-3,0]*1e3
        data3['biases'][bias] = data3['fixed_radiation'][bias] - obs['fixed_radiation'][bias + '_ship'][drift_ship[0]]
        mask = ~np.isnan(data3['biases']['LWP']) & ~np.isnan(data3['biases'][bias])
        data3['biases'][bias + '-LWP_regres']['slope'], data3['biases'][bias + '-LWP_regres']['intercept'], data3['biases'][bias + '-LWP_regres']['r_value'], data3['biases'][bias + '-LWP_regres']['p_value'], data3['biases'][bias + '-LWP_regres']['std_err'] = stats.linregress(data3['biases']['LWP'][mask],data3['biases'][bias][mask])
        data3['biases'][bias + '-LWP_regres']['line'] = data3['biases'][bias + '-LWP_regres']['slope'] * data3['biases']['LWP'] + data3['biases'][bias + '-LWP_regres']['intercept']
        print('r-' + bias + '-ecmwf_ifs:', data3['biases'][bias + '-LWP_regres']['r_value'])
        print('r-squared-' + bias + '-ecmwf_ifs:', data3['biases'][bias + '-LWP_regres']['r_value']**2)

        data4['biases'][bias + '-LWP_regres'] = {}
        data4['biases']['LWP'] = ra2t_data['model_lwp'][:-3]*1e3 - obs_data['lwp'][:-3,0]*1e3
        data4['biases'][bias] = data4['fixed_radiation'][bias] - obs['fixed_radiation'][bias + '_ship'][drift_ship[0]]
        mask = ~np.isnan(data4['biases']['LWP']) & ~np.isnan(data4['biases'][bias])
        data4['biases'][bias + '-LWP_regres']['slope'], data4['biases'][bias + '-LWP_regres']['intercept'], data4['biases'][bias + '-LWP_regres']['r_value'], data4['biases'][bias + '-LWP_regres']['p_value'], data4['biases'][bias + '-LWP_regres']['std_err'] = stats.linregress(data4['biases']['LWP'][mask],data4['biases'][bias][mask])
        data4['biases'][bias + '-LWP_regres']['line'] = data4['biases'][bias + '-LWP_regres']['slope'] * data4['biases']['LWP'] + data4['biases'][bias + '-LWP_regres']['intercept']
        print('r-' + bias + '-um_ra2t:', data4['biases'][bias + '-LWP_regres']['r_value'])
        print('r-squared-' + bias + '-um_ra2t:', data4['biases'][bias + '-LWP_regres']['r_value']**2)

    ##################################################
    ##################################################
    #### plot figure
    ##################################################
    ##################################################
    # ### -------------------------------
    # ### Build simple radiation-LWP bias figures
    # ### -------------------------------
    # SMALL_SIZE = 12
    # MED_SIZE = 16
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=MED_SIZE)
    # plt.rc('axes',titlesize=MED_SIZE)
    # plt.rc('axes',labelsize=MED_SIZE)
    # plt.rc('xtick',labelsize=MED_SIZE)
    # plt.rc('ytick',labelsize=MED_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)

    # fig = plt.figure(figsize=(14,7.5))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    #         hspace = 0.3, wspace = 0.3)
    #
    #
    # plt.subplot(231)
    # plt.plot(obs_data['lwp'][:-3,0]*1e3, obs['fixed_radiation']['SWd_ship'][drift_ship[0]], 'ko', alpha = 0.5, markersize = 4, label = 'Obs')
    # plt.plot(um_data['model_lwp'][:-3]*1e3,data1['fixed_radiation']['SWd'], '.', alpha = 0.5, color = 'darkblue', label = 'UM_RA2M')
    # plt.plot(ra2t_data['model_lwp'][:-3]*1e3,data4['fixed_radiation']['SWd'], '.', alpha = 0.5, color = 'steelblue', label = 'UM_RA2T')
    # plt.plot(misc_data['model_lwp'][:-3]*1e3,data2['fixed_radiation']['SWd'], '.', alpha = 0.5, color = 'mediumseagreen', label = 'UM_CASIM-100')
    # plt.plot(ifs_data['model_lwp'][:-3]*1e3,data3['fixed_radiation']['SWd'], '.', alpha = 0.5, color = 'gold', label = 'ECMWF_IFS')
    # plt.legend()
    # plt.xlabel('LWP [g m$^{-2}$]')
    # plt.ylabel('SW$_{\downarrow}$ [W m$^{-2}$]')
    #
    # plt.subplot(232)
    # plt.plot(obs_data['lwp'][:-3,0]*1e3, obs['fixed_radiation']['LWd_ship'][drift_ship[0]], 'ko', alpha = 0.5, markersize = 4)
    # plt.plot(um_data['model_lwp'][:-3]*1e3,data1['fixed_radiation']['LWd'], '.', alpha = 0.5, color = 'darkblue')
    # plt.plot(misc_data['model_lwp'][:-3]*1e3,data2['fixed_radiation']['LWd'], '.', alpha = 0.5, color = 'mediumseagreen')
    # plt.plot(ra2t_data['model_lwp'][:-3]*1e3,data4['fixed_radiation']['LWd'], '.', alpha = 0.5, color = 'steelblue')
    # plt.plot(ifs_data['model_lwp'][:-3]*1e3,data3['fixed_radiation']['LWd'], '.', alpha = 0.5, color = 'gold')
    # plt.xlabel('LWP [g m$^{-2}$]')
    # plt.ylabel('LW$_{\downarrow}$ [W m$^{-2}$]')
    #
    # plt.subplot(233)
    # plt.plot(obs_data['lwp'][:-3,0]*1e3, obs['fixed_radiation']['Rnet_ship'][drift_ship[0]], 'ko', alpha = 0.5, markersize = 4)
    # plt.plot(um_data['model_lwp'][:-3]*1e3,data1['fixed_radiation']['Rnet'], '.', alpha = 0.5, color = 'darkblue')
    # plt.plot(misc_data['model_lwp'][:-3]*1e3,data2['fixed_radiation']['Rnet'], '.', alpha = 0.5, color = 'mediumseagreen')
    # plt.plot(ra2t_data['model_lwp'][:-3]*1e3,data4['fixed_radiation']['Rnet'], '.', alpha = 0.5, color = 'steelblue')
    # plt.plot(ifs_data['model_lwp'][:-3]*1e3,data3['fixed_radiation']['Rnet'], '.', alpha = 0.5, color = 'gold')
    # plt.xlabel('LWP [g m$^{-2}$]')
    # plt.ylabel('R$_{net}$ [W m$^{-2}$]')
    #
    # plt.subplot(234)
    # alf = 0.3
    # xmin = -250
    # xmax = 400
    # y1min = -130
    # y1max = 75
    # plt.plot(data1['biases']['LWP'], data1['biases']['SWd'], '.', alpha = alf, color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['SWd'], '.', alpha = alf, color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['SWd'], '.', alpha = alf, color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['SWd'], '.', alpha = alf, color = 'steelblue')
    # plt.plot(data1['biases']['LWP'], data1['biases']['SWd-LWP_regres']['line'], color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['SWd-LWP_regres']['line'], color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['SWd-LWP_regres']['line'], color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['SWd-LWP_regres']['line'], color = 'steelblue')
    # # plt.grid('on')
    # plt.ylim([y1min,y1max])
    # plt.xlim([xmin,xmax])
    # plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.plot([0,0],[y1min,y1max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    # plt.ylabel('SW$_{\downarrow}$ [mod-obs]')
    #
    # plt.subplot(235)
    # y2min = -80
    # y2max = 90
    # plt.plot(data1['biases']['LWP'], data1['biases']['LWd'], '.', alpha = alf, color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['LWd'], '.', alpha = alf, color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['LWd'], '.', alpha = alf, color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['LWd'], '.', alpha = alf, color = 'steelblue')
    # plt.plot(data1['biases']['LWP'], data1['biases']['LWd-LWP_regres']['line'], color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['LWd-LWP_regres']['line'], color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['LWd-LWP_regres']['line'], color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['LWd-LWP_regres']['line'], color = 'steelblue')
    # # plt.grid('on')
    # plt.ylim([y2min,y2max])
    # plt.xlim([xmin,xmax])
    # plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    # plt.ylabel('LW$_{\downarrow}$ [mod-obs]')
    #
    # plt.subplot(236)
    # y3min = -80
    # y3max = 90
    # plt.plot(data1['biases']['LWP'], data1['biases']['Rnet'], '.', alpha = alf, color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['Rnet'], '.', alpha = alf, color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['Rnet'], '.', alpha = alf, color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['Rnet'], '.', alpha = alf, color = 'steelblue')
    # plt.plot(data1['biases']['LWP'], data1['biases']['Rnet-LWP_regres']['line'], color = 'darkblue')
    # plt.plot(data2['biases']['LWP'], data2['biases']['Rnet-LWP_regres']['line'], color = 'mediumseagreen')
    # plt.plot(data3['biases']['LWP'], data3['biases']['Rnet-LWP_regres']['line'], color = 'gold')
    # plt.plot(data4['biases']['LWP'], data4['biases']['Rnet-LWP_regres']['line'], color = 'steelblue')
    # # plt.grid('on')
    # # plt.ylim([y2min,y2max])
    # # plt.xlim([xmin,xmax])
    # plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    # plt.ylabel('R$_{net}$ [mod-obs]')
    #
    # plt.savefig('../FIGS/comparisons/Radiation-LWP_Correlations.png', dpi = 300)
    # plt.show()

    ### -------------------------------
    ### Build simple radiation-LWP bias figures
    ### -------------------------------

    ## create index for below 3km
    lt3km_um = np.where(um_data['height'][0,:] <= 1000)
    lt3km_ifs = np.where(ifs_data['height'][0,:] <= 1000)

    tempvar1 = np.nanmean(mask1[:,lt3km_um[0]],0) - np.nanmean(mask0[:,lt3km_um[0]],0)
    dZ_um = um_data['height'][0,lt3km_um[0][1:]] - um_data['height'][0,lt3km_um[0][:-1]]

    # print (tempvar1.shape)

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    fig = plt.figure(figsize=(14,14))
    plt.subplots_adjust(top = 0.95, bottom = 0.08, right = 0.95, left = 0.1,
            hspace = 0.35, wspace = 0.35)

    alf = 0.3
    xmin = -250
    xmax = 400
    y1min = -130
    y1max = 75
    y2min = -80
    y2max = 90
    y3min = -80
    y3max = 90

    vminn = 0.0
    vmaxx = 0.5

    cloud_fractions = [obs_data['Cv'][:-3,:], um_data['model_Cv_filtered'][:-3,:], misc_data['model_Cv_filtered'][:-3,:], ifs_data['model_snow_Cv_filtered'][:-3,:], ra2t_data['model_Cv_filtered'][:-3,:]]
    cloud_fractions_lt3km = [obs_data['Cv'][:-3,lt3km_um[0]], um_data['model_Cv_filtered'][:-3,lt3km_um[0]], misc_data['model_Cv_filtered'][:-3,lt3km_um[0]], ifs_data['model_snow_Cv_filtered'][:-3,lt3km_um[0]], ra2t_data['model_Cv_filtered'][:-3,lt3km_um[0]]]
    cloud_masks = [mask0[:-3,:], mask1[:-3,:], mask2[:-3,:], mask3[:-3,:], mask4[:-3,:]]
    cloud_masks_lt3km = [mask0[:-3,lt3km_um[0]], mask1[:-3,lt3km_um[0]], mask2[:-3,lt3km_um[0]], mask3[:-3,lt3km_ifs[0]], mask4[:-3,lt3km_um[0]]]
    lwc = [obs_data['lwc_adiabatic'][:-3,:]*1e3, um_data['model_lwc'][:-3,:]*1e3, misc_data['model_lwc'][:-3,:]*1e3, ifs_data['model_lwc'][:-3,:]*1e3, ra2t_data['model_lwc'][:-3,:]*1e3]
    iwc = [obs_data['iwc'][:-3,:]*1e3, um_data['model_iwc'][:-3,:]*1e3, misc_data['model_iwc'][:-3,:]*1e3, ifs_data['model_iwc_filtered'][:-3,:]*1e3, ra2t_data['model_iwc'][:-3,:]*1e3]

    iwv = [obs['hatpro']['IWV'][drift_ship[0]], 0, 0, 0, 0]

    print (obs['hatpro']['IWV'][drift_ship[0]].shape)
    print (data1['biases']['LWP'].shape)
    print (data1['biases']['SWd'].shape)
    print (obs_data['iwc'][:-3,:].shape)

    lwc_biases = [obs_data['lwc_adiabatic'][:-3,:]*1e3 - obs_data['lwc_adiabatic'][:-3,:]*1e3, um_data['model_lwc'][:-3,:]*1e3 - obs_data['lwc_adiabatic'][:-3,:]*1e3,
        misc_data['model_lwc'][:-3,:]*1e3 - obs_data['lwc_adiabatic'][:-3,:]*1e3, ifs_data['model_lwc'][:-3,:]*1e3,
        ra2t_data['model_lwc'][:-3,:]*1e3 - obs_data['lwc_adiabatic'][:-3,:]*1e3]

    var = lwc_biases
    # arg = [var[0], var[1], var[2], var[3], var[4]]    ## 1D timeseries data
    arg = [np.nanmean(var[0],1), np.nanmean(var[1],1), np.nanmean(var[2],1), np.nanmean(var[3],1), np.nanmean(var[4],1)]    ## means of cloudnet profile data

    print (np.nanmean(var[1],1).shape)

    ####------          SWd
    plt.subplot(431)
    plt.scatter(data1['biases']['LWP'], data1['biases']['SWd'], c = arg[1], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data1['biases']['LWP'], data1['biases']['SWd-LWP_regres']['line'], color = 'darkblue', zorder = 3)
    plt.ylim([y1min,y1max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y1min,y1max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_RA2M', fontsize = 12)
    plt.ylabel('SW$_{\downarrow}$ [mod-obs]')

    plt.subplot(434)
    plt.scatter(data2['biases']['LWP'], data2['biases']['SWd'], c = arg[2], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data2['biases']['LWP'], data2['biases']['SWd-LWP_regres']['line'], color = 'mediumseagreen', zorder = 3)
    plt.ylim([y1min,y1max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y1min,y1max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_CASIM-100', fontsize = 12)
    plt.ylabel('SW$_{\downarrow}$ [mod-obs]')

    plt.subplot(437)
    plt.scatter(data3['biases']['LWP'], data3['biases']['SWd'], c = arg[3], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data3['biases']['LWP'], data3['biases']['SWd-LWP_regres']['line'], color = 'gold', zorder = 3)
    plt.ylim([y1min,y1max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y1min,y1max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('ECMWF_IFS', fontsize = 12)
    plt.ylabel('SW$_{\downarrow}$ [mod-obs]')

    plt.subplot(4,3,10)
    plt.scatter(data4['biases']['LWP'], data4['biases']['SWd'], c = arg[4], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data4['biases']['LWP'], data4['biases']['SWd-LWP_regres']['line'], color = 'steelblue', zorder = 3)
    plt.ylim([y1min,y1max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y1min,y1max],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.title('UM_RA2T', fontsize = 12)
    plt.xlabel('LWP [mod-obs]')
    plt.ylabel('SW$_{\downarrow}$ [mod-obs]')

    ####------          LWd
    plt.subplot(432)
    plt.scatter(data1['biases']['LWP'], data1['biases']['LWd'], c = arg[1], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data1['biases']['LWP'], data1['biases']['LWd-LWP_regres']['line'], color = 'darkblue', zorder = 3)
    plt.ylim([y2min,y2max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_RA2M', fontsize = 12)
    plt.ylabel('LW$_{\downarrow}$ [mod-obs]')

    plt.subplot(435)
    plt.scatter(data2['biases']['LWP'], data2['biases']['LWd'], c = arg[2], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data2['biases']['LWP'], data2['biases']['LWd-LWP_regres']['line'], color = 'mediumseagreen', zorder = 3)
    plt.ylim([y2min,y2max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_CASIM-100', fontsize = 12)
    plt.ylabel('LW$_{\downarrow}$ [mod-obs]')

    plt.subplot(438)
    plt.scatter(data3['biases']['LWP'], data3['biases']['LWd'], c = arg[3], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data3['biases']['LWP'], data3['biases']['LWd-LWP_regres']['line'], color = 'gold', zorder = 3)
    plt.ylim([y2min,y2max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('ECMWF_IFS', fontsize = 12)
    plt.ylabel('LW$_{\downarrow}$ [mod-obs]')

    plt.subplot(4,3,11)
    plt.scatter(data4['biases']['LWP'], data4['biases']['LWd'], c = arg[4], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data4['biases']['LWP'], data4['biases']['LWd-LWP_regres']['line'], color = 'steelblue', zorder = 3)
    plt.ylim([y2min,y2max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y2min,y2max],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.title('UM_RA2T', fontsize = 12)
    plt.xlabel('LWP [mod-obs]')
    plt.ylabel('LW$_{\downarrow}$ [mod-obs]')

    ####------          Rnet
    plt.subplot(433)
    plt.scatter(data1['biases']['LWP'], data1['biases']['Rnet'], c = arg[1], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data1['biases']['LWP'], data1['biases']['Rnet-LWP_regres']['line'], color = 'darkblue', zorder = 3)
    plt.ylim([y3min,y3max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y3min,y3max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_RA2M', fontsize = 12)
    plt.ylabel('R$_{net}$ [mod-obs]')

    plt.subplot(436)
    plt.scatter(data2['biases']['LWP'], data2['biases']['Rnet'], c = arg[2], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data2['biases']['LWP'], data2['biases']['Rnet-LWP_regres']['line'], color = 'mediumseagreen', zorder = 3)
    plt.ylim([y3min,y3max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y3min,y3max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('UM_CASIM-100', fontsize = 12)
    plt.ylabel('R$_{net}$ [mod-obs]')

    plt.subplot(439)
    plt.scatter(data3['biases']['LWP'], data3['biases']['Rnet'], c = arg[3], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data3['biases']['LWP'], data3['biases']['Rnet-LWP_regres']['line'], color = 'gold', zorder = 3)
    plt.ylim([y3min,y3max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y3min,y3max],'--', color='lightgrey', linewidth = 1, zorder=0)
    # plt.xlabel('LWP [mod-obs]')
    plt.title('ECMWF_IFS', fontsize = 12)
    plt.ylabel('R$_{net}$ [mod-obs]')

    plt.subplot(4,3,12)
    plt.scatter(data4['biases']['LWP'], data4['biases']['Rnet'], c = arg[4], s = 6,
        # vmin = vminn, vmax = vmaxx,
        )
    plt.plot(data4['biases']['LWP'], data4['biases']['Rnet-LWP_regres']['line'], color = 'steelblue', zorder = 3)
    plt.ylim([y3min,y3max])
    plt.xlim([xmin,xmax])
    plt.plot([xmin,xmax],[0,0],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.plot([0,0],[y3min,y3max],'--', color='lightgrey', linewidth = 1, zorder=0)
    plt.title('UM_RA2T', fontsize = 12)
    plt.xlabel('LWP [mod-obs]')
    plt.ylabel('R$_{net}$ [mod-obs]')

    # plt.savefig('../FIGS/comparisons/Radiation-LWP_Correlations_cTWCMask-lt3km.png', dpi = 300)
    plt.show()

def plot_ObsGridComparison(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy): #, lon, lat):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Comparing observed Cv based on grid results interpolated on to:')
    print ('')

    ##################################################
    ##################################################
    #### 	CARTOPY
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
    plt.rc('legend',fontsize=LARGE_SIZE)
    plt.figure(figsize=(4,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.94, left = 0.15,
            hspace = 0.22, wspace = 0.15)

    # print um_data.keys()

    #### set flagged um_data to nans
    ifs_data['Cv'][ifs_data['Cv'] == -999] = np.nan
    um_data['Cv'][um_data['Cv'] == -999] = np.nan
    obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    ra2t_data['Cv'][ra2t_data['Cv'] == -999] = np.nan

    ax2 = plt.gca()
    ax2.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0),
        np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(obs_data['Cv'],0),np.nanmean(obs_data['height'],0), 'k--', linewidth = 3, label = 'Obs', zorder = 3)
    plt.plot(np.nanmean(obs_data['Cv'],0) - np.nanstd(obs_data['Cv'],0), np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(obs_data['Cv'],0) + np.nanstd(obs_data['Cv'],0), np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(ifs_data['height'],0),np.nanmean(ifs_data['Cv'],0) - np.nanstd(ifs_data['Cv'],0),
        np.nanmean(ifs_data['Cv'],0) + np.nanstd(ifs_data['Cv'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(ifs_data['Cv'],0),np.nanmean(ifs_data['height'],0), color = 'gold', linewidth = 3, label = 'Obs_IFS')
    plt.plot(np.nanmean(ifs_data['Cv'],0) - np.nanstd(ifs_data['Cv'],0), np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(ifs_data['Cv'],0) + np.nanstd(ifs_data['Cv'],0), np.nanmean(ifs_data['height'],0),
        '--', color = 'gold', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(ra2t_data['height'],0),np.nanmean(ra2t_data['Cv'],0) - np.nanstd(ra2t_data['Cv'],0),
        np.nanmean(ra2t_data['Cv'],0) + np.nanstd(ra2t_data['Cv'],0), color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(ra2t_data['Cv'],0),np.nanmean(ra2t_data['height'],0), color = 'steelblue', linewidth = 3, label = 'Obs_UM')
    plt.plot(np.nanmean(ra2t_data['Cv'],0) - np.nanstd(ra2t_data['Cv'],0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(ra2t_data['Cv'],0) + np.nanstd(ra2t_data['Cv'],0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.xlabel('Cloud Fraction')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax2.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,1])
    plt.legend()

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    if month_flag == -1:
        fileout = 'FIGS/Obs_GridComparison_Cv_226-257DOY_wMissingFiles_newColours.svg'
    plt.savefig(fileout)
    plt.show()

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
    if np.logical_or(out_dir3 == 'OUT_25H/',out_dir3 == 'ECMWF_IFS/'):
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
    # plt.savefig(fileout, dpi=300)
    plt.show()

def plot_scaledBLCv_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Scaling cloudnet data by identified thetaE inversions:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    # #################################################################
    # ## save data into temp variables to allow subsampling
    # #################################################################
    # bldepth1 = data1['bl_depth'][data1['hrly_flag']]
    # bldepth2 = data2['bl_depth'][data2['hrly_flag']]
    # if ifs_flag == True:
    #     bldepth3 = data3['sfc_bl_height'][data3['hrly_flag']]
    # else:
    #     bldepth3 = data3['bl_depth'][data3['hrly_flag']]

    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

    # #### ---------------------------------------------------------------
    # #### ONLY LOOK AT SONDES FROM THE DRIFT
    # #### ---------------------------------------------------------------
    # drift = np.where(np.logical_and(obs['inversions']['thetaE']['time'] >= 225.9, obs['inversions']['thetaE']['time'] <= 258.0))

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = obs['inversions']['thetaE']['invbase']
    obsmlh = obs['inversions']['thetaE']['sfmlheight']

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    temp_time = obs['inversions']['thetaE']['time']         #### temporary time array for indexing
    temp_time2 = np.zeros(len(temp_time))
    temp_time2[:] = np.nan
    temp_time2[:nanindices[0]] = temp_time[:nanindices[0]]
    temp_time2[nanindices[3]+1:nanindices[4]] = temp_time[nanindices[3]+1:nanindices[4]]
    temp_time2[nanindices[7]+1:nanindices[8]] = temp_time[nanindices[7]+1:nanindices[8]]
    temp_inv = np.zeros(len(obsinv))
    temp_inv[:] = np.nan
    temp_inv[:nanindices[0]] = obsinv[:nanindices[0]]
    temp_inv[nanindices[3]+1:nanindices[4]] = obsinv[nanindices[3]+1:nanindices[4]]
    temp_inv[nanindices[7]+1:nanindices[8]] = obsinv[nanindices[7]+1:nanindices[8]]
    temp_sfml = np.zeros(len(obsmlh))
    temp_sfml[:] = np.nan
    temp_sfml[:nanindices[0]] = obsmlh[:nanindices[0]]
    temp_sfml[nanindices[3]+1:nanindices[4]] = obsmlh[nanindices[3]+1:nanindices[4]]
    temp_sfml[nanindices[7]+1:nanindices[8]] = obsmlh[nanindices[7]+1:nanindices[8]]

    ### reset time array with new one accounting for missing FILES
    obs['inversions']['thetaE']['time'] = temp_time2
    obsinv = temp_inv
    obsmlh = temp_sfml

    ### save non-nan values to dictionary
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['thetaE']['time'][~np.isnan(obs['inversions']['thetaE']['time'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data['Cv_6hrly'] = obs_data['Cv'][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); obssfml[:] = np.nan

    ### look for altitudes < invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    plt.figure()
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),obs['inversions']['thetaE']['invbase'])
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaledCv'] = {}
    obs['inversions']['scaledCv']['binned'] = {}
    obs['inversions']['scaledCv']['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['mean'][:] = np.nan
    obs['inversions']['scaledCv']['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['blCv'] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['blCv'][:] = np.nan

    ### fill arrays with cloudnet data below each invbase
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        ###         +1 includes invbase in array
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['blCv'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)] = obs_data['Cv_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['blCv'][i,tempvar]
            if np.size(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaledCv']['mean'][i,k] = np.nanmean(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaledCv']['stdev'][i,k] = np.nanstd(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    ####        data from thetaE algorithm is on radiosonde (6 hourly) timesteps already
    #### ---------------------------------------------------------------

    ### need to make sure times where we don't have an inversion (invbase == nan) in the IFS doesn't mess up the algorithm
    ###     set these cases to zero for now
    data3['inversions']['invbase'][np.isnan(data3['inversions']['invbase'])] = 0.0

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    ####        time
    tim1 = np.zeros(len(data1['inversions']['time']))
    tim1[:] = np.nan
    tim1[:nanindices[0]] = data1['inversions']['time'][:nanindices[0]]
    tim1[nanindices[3]+1:nanindices[4]] = data1['inversions']['time'][nanindices[3]+1:nanindices[4]]
    tim1[nanindices[7]+1:nanindices[8]] = data1['inversions']['time'][nanindices[7]+1:nanindices[8]]
    tim2 = tim1
    tim3 = tim1
    ####        inversions
    inv1 = np.zeros(len(data1['inversions']['invbase']))
    inv1[:] = np.nan
    inv1[:nanindices[0]] = data1['inversions']['invbase'][:nanindices[0]]
    inv1[nanindices[3]+1:nanindices[4]] = data1['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv1[nanindices[7]+1:nanindices[8]] = data1['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv2 = np.zeros(len(data2['inversions']['invbase']))
    inv2[:] = np.nan
    inv2[:nanindices[0]] = data2['inversions']['invbase'][:nanindices[0]]
    inv2[nanindices[3]+1:nanindices[4]] = data2['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv2[nanindices[7]+1:nanindices[8]] = data2['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv3 = np.zeros(len(data3['inversions']['invbase']))
    inv3[:] = np.nan
    inv3[:nanindices[0]] = data3['inversions']['invbase'][:nanindices[0]]
    inv3[nanindices[3]+1:nanindices[4]] = data3['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv3[nanindices[7]+1:nanindices[8]] = data3['inversions']['invbase'][nanindices[7]+1:nanindices[8]]

    ### save non-nan values
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    tim1 = tim1[~np.isnan(tim1)]
    tim2 = tim2[~np.isnan(tim2)]
    tim3 = tim3[~np.isnan(tim3)]
    inv1 = inv1[~np.isnan(inv1)]
    inv2 = inv2[~np.isnan(inv2)]
    inv3 = inv3[~np.isnan(inv3)]

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    # #### ---------------------------------------------------------------
    # #### remove flagged IFS heights
    # #### ---------------------------------------------------------------
    # data3['height'][data3['height'] == -9999] = 0.0
    #         #### set all heights to zero if flagged. setting to nan caused problems further on
    # data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    # #### ---------------------------------------------------------------
    # #### Define meteorological periods from Jutta's paper
    # #### ---------------------------------------------------------------
    #
    # ## Meteorological period definitions from Jutta's paper:
    # ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    # ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    # ##     radar measurements were not possible during heavy ice breaking because of excessive
    # ##     vibration, cloud characteristics and fog heights are not available during period 2.
    # ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    # ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    # ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    # ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    # ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    # ##     out to the ice edge. "
    # #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...
    #
    # ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    # p3 = np.where(um_data['time'] < 230.0)
    # p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    # p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    # p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    # p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    # p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time'].data[::6]
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time'].data[::6]
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time'].data[::6]

    #### define empty arrays of nans to fill with scaled data
    data1['scaledCv'] = {}
    data1['scaledCv']['binned'] = {}
    data1['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaledCv']['mean'][:] = np.nan
    data1['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaledCv']['stdev'][:] = np.nan
    data1['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(um_data['height'],1)]); data1['blCv'][:] = np.nan

    data2['scaledCv'] = {}
    data2['scaledCv']['binned'] = {}
    data2['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaledCv']['mean'][:] = np.nan
    data2['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaledCv']['stdev'][:] = np.nan
    data2['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(misc_data['height'],1)]); data2['blCv'][:] = np.nan

    data3['scaledCv'] = {}
    data3['scaledCv']['binned'] = {}
    data3['scaledCv']['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaledCv']['mean'][:] = np.nan
    data3['scaledCv']['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaledCv']['stdev'][:] = np.nan
    data3['blCv'] = np.zeros([np.size(data1['scaledTime']),np.size(ifs_data['height'],1)]); data3['blCv'][:] = np.nan

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data1['scaledCv']['inversion_Tindex'][:] = np.nan
    data1['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data1['scaledCv']['inversionForCloudnet'][:] = np.nan
    data2['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data2['scaledCv']['inversion_Tindex'][:] = np.nan
    data2['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data2['scaledCv']['inversionForCloudnet'][:] = np.nan
    data3['scaledCv']['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data3['scaledCv']['inversion_Tindex'][:] = np.nan
    data3['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data3['scaledCv']['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.size(np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))[0][0]
        if np.size(np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))[0][0]
        if np.size(np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaledCv']['inversion_Tindex'][i] = np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data1['scaledCv']['inversionForCloudnet'][int(data1['scaledCv']['inversion_Tindex'][i])] = inv1[i]
        if data2['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data2['scaledCv']['inversionForCloudnet'][int(data2['scaledCv']['inversion_Tindex'][i])] = inv2[i]
        if data3['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data3['scaledCv']['inversionForCloudnet'][int(data3['scaledCv']['inversion_Tindex'][i])] = inv3[i]
    #
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    #
    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(data1['scaledTime'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(data2['scaledTime'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(data3['scaledTime'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(data1['scaledTime'])):        ### all can go in this loop, data1['scaledTime'] == 6-hourly data

        ### main inversion base assignments
        ###         (1) find where UM height array matches the invbase index for index i
        ###         (2) find where UM height array matches the invbase index for index i
        ###         (3) find where IFS height array is less than or equal to the UM-gridded invbase index for index i
        if np.size(np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(misc_data['height'][i,:].data == data2['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(misc_data['height'][i,:].data == data2['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(ifs_data['height'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            temp = ifs_data['height'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    plt.figure()
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(data1['scaledTime'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1,inv1)
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(data2['scaledTime'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2, inv2)
    plt.ylim([0,3e3])
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(data3['scaledTime'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3, inv3)
    plt.ylim([0,3e3])
    plt.show()

    ### set 6 hourly cloudnet Cv arrays as tempvars
    ra2m_Cv = um_data['model_Cv_filtered'][::6,:]
    casim_Cv = misc_data['model_Cv_filtered'][::6,:]
    ifs_Cv = ifs_data['model_snow_Cv_filtered'][::6,:]

    ### find all Cv data below identified inversion
    for i in range(0,np.size(data1['scaledTime'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaledCv']['binned']['t' + str(i)] = {}
        data2['scaledCv']['binned']['t' + str(i)] = {}
        data3['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['blCv'][i,:int(data1['inversions']['invbase_kIndex'][i]+1)] = ra2m_Cv[i,:int(data1['inversions']['invbase_kIndex'][i]+1)]
        data2['blCv'][i,:int(data2['inversions']['invbase_kIndex'][i]+1)] = casim_Cv[i,:int(data2['inversions']['invbase_kIndex'][i]+1)]
        data3['blCv'][i,:int(data3['inversions']['invbase_kIndex'][i]+1)] = ifs_Cv[i,:int(data3['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            data1['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data1['blCv'][i,tempvar1]
            if np.size(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaledCv']['mean'][i,k] = np.nanmean(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data1['scaledCv']['stdev'][i,k] = np.nanstd(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data2['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data2['blCv'][i,tempvar2]
            if np.size(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaledCv']['mean'][i,k] = np.nanmean(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data2['scaledCv']['stdev'][i,k] = np.nanstd(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data3['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data3['blCv'][i,tempvar3]
            if np.size(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaledCv']['mean'][i,k] = np.nanmean(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data3['scaledCv']['stdev'][i,k] = np.nanstd(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    ### save working data for debug
    np.save('working_data1', data1)

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.show()

    ### obs
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['blCv'])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),np.squeeze(obs['inversions']['thetaE']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_ra2m
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],um_data['height'][0,:].data,np.transpose(data1['blCv'])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['time'],data1['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_casim-100
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],misc_data['height'][0,:].data,np.transpose(data2['blCv'])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['time'],data2['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### ecmwf_ifs
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],ifs_data['height'][0,:].data,np.transpose(data3['blCv'])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['time'],data3['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### profiles
    ###         loop through and set all zeros to nans
    obsmean = obs['inversions']['scaledCv']['mean']
    ra2mmean = data1['scaledCv']['mean']
    casimmean = data2['scaledCv']['mean']
    ifsmean = data3['scaledCv']['mean']
    # for i in range(0, len(data1['scaledTime'])):
    #     obsmean[i,obs['inversions']['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     ra2mmean[i,data1['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     casimmean[i,data2['scaledCv']['mean'][i,:] == 0.0] = np.nan
    #     ifsmean[i,data3['scaledCv']['mean'][i,:] == 0.0] = np.nan
    plt.plot(np.nanmean(obsmean,0),obs['inversions']['scaledZ'], '--', color = 'k', linewidth = 2, label = 'Obs')
    plt.plot(np.nanmean(ra2mmean,0),data1['scaledZ'], '^-', color = 'darkblue', linewidth = 2, label = label1)
    plt.plot(np.nanmean(casimmean,0),data2['scaledZ'], 'v-', color = 'mediumseagreen', linewidth = 2, label = label2)
    plt.plot(np.nanmean(ifsmean,0),data3['scaledZ'], 'd-', color = 'gold', linewidth = 2, label = label3)
    plt.legend()
    plt.show()

def plot_scaledBLCv_JVInv(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3):

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


    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
    um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan

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

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = np.squeeze(obs['inversions']['invbase'][drift])
    obsmlh = np.squeeze(obs['inversions']['sfmlheight'][drift])

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])
    obs['inversions']['doy_drift'][nanindices] = np.nan
    obsinv[nanindices] = np.nan
    obsmlh[nanindices] = np.nan

    ### save non-nan values to dictionary
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['doy_drift'][~np.isnan(obs['inversions']['doy_drift'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data['Cv_6hrly'] = obs_data['Cv'][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); obssfml[:] = np.nan

    ### look for altitudes < invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    plt.figure()
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']))
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaledCv'] = {}
    obs['inversions']['scaledCv']['binned'] = {}
    obs['inversions']['scaledCv']['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['mean'][:] = np.nan
    obs['inversions']['scaledCv']['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaledCv']['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['blCv'] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['blCv'][:] = np.nan

    ###
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['blCv'][i,:int(obs['inversions']['invbase_kIndex'][i])] = obs_data['Cv_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i])]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['blCv'][i,tempvar]
            if np.size(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaledCv']['mean'][i,k] = np.nanmean(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaledCv']['stdev'][i,k] = np.nanstd(obs['inversions']['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    #### ---------------------------------------------------------------
    #### ---------------------------------------------------------------
    #### ONLY LOOK AT DATA FROM THE DRIFT
    #### ---------------------------------------------------------------
    driftmod = np.where(np.logical_and(data1['inversions']['doy'] >= 226.0, data1['inversions']['doy'] <= 258.0))
    data1['inversions']['doy_drift'] = data1['inversions']['doy'][driftmod[0][1:-2]]
    data1['inversions']['invbase_drift'] = data1['inversions']['invbase'][driftmod[0][1:-2]]
    data2['inversions']['doy_drift'] = data2['inversions']['doy'][driftmod[0][1:-2]]
    data2['inversions']['invbase_drift'] = data2['inversions']['invbase'][driftmod[0][1:-2]]
    data3['inversions']['doy_drift'] = data3['inversions']['doy'][driftmod[0][1:-2]]
    data3['inversions']['invbase_drift'] = data3['inversions']['invbase'][driftmod[0][1:-2]]

    #### make inversion tempvars to allow for easy subsampling
    tim1 = np.squeeze(data1['inversions']['doy_drift'][data1['hrly_flag']])
    tim2 = np.squeeze(data2['inversions']['doy_drift'][data2['hrly_flag']])
    tim3 = np.squeeze(data3['inversions']['doy_drift'][data3['hrly_flag']])
    inv1 = np.squeeze(data1['inversions']['invbase_drift'][data1['hrly_flag'],0])
    inv2 = np.squeeze(data2['inversions']['invbase_drift'][data2['hrly_flag'],0])
    inv3 = np.squeeze(data3['inversions']['invbase_drift'][data3['hrly_flag'],0])

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    #### ---------------------------------------------------------------
    #### remove flagged IFS heights
    #### ---------------------------------------------------------------
    data3['height'][data3['height'] == -9999] = 0.0
            #### set all heights to zero if flagged. setting to nan caused problems further on
    data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    #### ---------------------------------------------------------------
    #### Define meteorological periods from Jutta's paper
    #### ---------------------------------------------------------------

    ## Meteorological period definitions from Jutta's paper:
    ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    ##     radar measurements were not possible during heavy ice breaking because of excessive
    ##     vibration, cloud characteristics and fog heights are not available during period 2.
    ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    ##     out to the ice edge. "
    #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...

    ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    p3 = np.where(um_data['time'] < 230.0)
    p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    #### define empty arrays of nans to fill with scaled data
    data1['scaledCv'] = {}
    data1['scaledCv']['binned'] = {}
    data1['scaledCv']['mean'] = np.zeros([np.size(um_data['height'],0),len(Zpts)]); data1['scaledCv']['mean'][:] = np.nan
    data1['scaledCv']['stdev'] = np.zeros([np.size(um_data['height'],0),len(Zpts)]); data1['scaledCv']['stdev'][:] = np.nan
    data1['blCv'] = np.zeros([np.size(um_data['height'],0),np.size(um_data['height'],1)]); data1['blCv'][:] = np.nan

    data2['scaledCv'] = {}
    data2['scaledCv']['binned'] = {}
    data2['scaledCv']['mean'] = np.zeros([np.size(misc_data['height'],0),len(Zpts)]); data2['scaledCv']['mean'][:] = np.nan
    data2['scaledCv']['stdev'] = np.zeros([np.size(misc_data['height'],0),len(Zpts)]); data2['scaledCv']['stdev'][:] = np.nan
    data2['blCv'] = np.zeros([np.size(misc_data['height'],0),np.size(misc_data['height'],1)]); data2['blCv'][:] = np.nan

    data3['scaledCv'] = {}
    data3['scaledCv']['binned'] = {}
    data3['scaledCv']['mean'] = np.zeros([np.size(ifs_data['height'],0),len(Zpts)]); data3['scaledCv']['mean'][:] = np.nan
    data3['scaledCv']['stdev'] = np.zeros([np.size(ifs_data['height'],0),len(Zpts)]); data3['scaledCv']['stdev'][:] = np.nan
    data3['blCv'] = np.zeros([np.size(ifs_data['height'],0),np.size(ifs_data['height'],1)]); data3['blCv'][:] = np.nan

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time']
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time']
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time']

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim1)); data1['scaledCv']['inversion_Tindex'][:] = np.nan
    data1['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim1)); data1['scaledCv']['inversionForCloudnet'][:] = np.nan
    data2['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim2)); data2['scaledCv']['inversion_Tindex'][:] = np.nan
    data2['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim2)); data2['scaledCv']['inversionForCloudnet'][:] = np.nan
    data3['scaledCv']['inversion_Tindex'] = np.zeros(np.size(tim3)); data3['scaledCv']['inversion_Tindex'][:] = np.nan
    data3['scaledCv']['inversionForCloudnet'] = np.zeros(np.size(tim3)); data3['scaledCv']['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.floor(tim1[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv1[i] = np.nan
            continue
        elif np.size(np.where(np.round(um_data['time'].data,3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaledCv']['inversion_Tindex'][i] = np.where(np.round(um_data['time'].data,3) == np.round(tim1[i],3))[0][0]
        if np.floor(tim2[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv2[i] = np.nan
            continue
        elif np.size(np.where(np.round(misc_data['time'].data,3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaledCv']['inversion_Tindex'][i] = np.where(np.round(misc_data['time'].data,3) == np.round(tim2[i],3))[0][0]
        if np.floor(tim3[i]) in missing_files:
            ### don't assign a cloudnet index for missing files
            ### instead, set inv1 to nan at these times
            inv3[i] = np.nan
            continue
        elif np.size(np.where(np.round(ifs_data['time'].data,3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaledCv']['inversion_Tindex'][i] = np.where(np.round(ifs_data['time'].data,3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data1['scaledCv']['inversionForCloudnet'][int(data1['scaledCv']['inversion_Tindex'][i])] = inv1[i]
        if data2['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data2['scaledCv']['inversionForCloudnet'][int(data2['scaledCv']['inversion_Tindex'][i])] = inv2[i]
        if data3['scaledCv']['inversion_Tindex'][i] >= 0.0:
            data3['scaledCv']['inversionForCloudnet'][int(data3['scaledCv']['inversion_Tindex'][i])] = inv3[i]

    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)

    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(um_data['time'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(misc_data['time'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(ifs_data['time'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(um_data['time'])):        ### all can go in this loop, um_data['time'] == hourly data

        ### main inversion base assignments
        if np.size(np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(data2['height'][1:].data == data2['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(data2['height'][1:].data == data2['scaledCv']['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(data3['height_hrly'][i].data <= data3['scaledCv']['inversionForCloudnet'][i])) > 0.0:
            temp = data3['height_hrly'][i,:].data <= data3['scaledCv']['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    plt.figure()
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(um_data['time'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1[:-22],inv1[:-22])
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(misc_data['time'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2[:-22], inv2[:-22])
    plt.ylim([0,3e3])
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(ifs_data['time'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3[:-22], inv3[:-22])
    plt.ylim([0,3e3])
    plt.show()

    # print (zind3)
    #### re-check inversion algorithm success rate to make sure no !0 values got dropped
    zzind1 = np.where(data1['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind1rate = np.size(zzind1) / np.float(np.size(inv1)) * 100.0
    zzind2 = np.where(data2['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind2rate = np.size(zzind2) / np.float(np.size(inv2)) * 100.0
    zzind3 = np.where(data3['inversions']['invbase_kIndex'] >= 0.0)  ## non-nan values
    zind3rate = np.size(zzind3) / np.float(np.size(inv3)) * 100.0

    ### try i = 0 first to see if it works
    ### this will go into a loop once tested
    # i = 110
    for i in range(0,np.size(um_data['time'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaledCv']['binned']['t' + str(i)] = {}
        data2['scaledCv']['binned']['t' + str(i)] = {}
        data3['scaledCv']['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['blCv'][i,:int(data1['inversions']['invbase_kIndex'][i])] = um_data['model_Cv_filtered'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        data2['blCv'][i,:int(data2['inversions']['invbase_kIndex'][i])] = misc_data['model_Cv_filtered'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        data3['blCv'][i,:int(data3['inversions']['invbase_kIndex'][i])] = ifs_data['model_snow_Cv_filtered'][i,:int(data3['inversions']['invbase_kIndex'][i])]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            data1['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data1['blCv'][i,tempvar1]
            if np.size(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaledCv']['mean'][i,k] = np.nanmean(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data1['scaledCv']['stdev'][i,k] = np.nanstd(data1['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data2['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data2['blCv'][i,tempvar2]
            if np.size(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaledCv']['mean'][i,k] = np.nanmean(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data2['scaledCv']['stdev'][i,k] = np.nanstd(data2['scaledCv']['binned']['t' + str(i)][Zpts[k]])

            data3['scaledCv']['binned']['t' + str(i)][Zpts[k]] = data3['blCv'][i,tempvar3]
            if np.size(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaledCv']['mean'][i,k] = np.nanmean(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])
            data3['scaledCv']['stdev'][i,k] = np.nanstd(data3['scaledCv']['binned']['t' + str(i)][Zpts[k]])

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean']), vmin = 0, vmax = 1)
    plt.show()

    ### obs
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['blCv'])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['doy_drift']),np.squeeze(obs['inversions']['invbase'][drift]),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_ra2m
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(um_data['time'].data,um_data['height'][0,:].data,np.transpose(data1['blCv'])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['doy_drift'],np.squeeze(data1['inversions']['invbase_drift']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### um_casim-100
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(misc_data['time'].data,misc_data['height'][0,:].data,np.transpose(data2['blCv'])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['doy'],np.squeeze(data2['inversions']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### ecmwf_ifs
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(ifs_data['time'].data,ifs_data['height'][0,:].data,np.transpose(data3['blCv'])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['doy'],np.squeeze(data3['inversions']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaledCv']['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.show()

    ### profiles
    plt.plot(np.nanmean(obs['inversions']['scaledCv']['mean'],0),obs['inversions']['scaledZ'], color = 'k', linewidth = 2, label = 'Obs')
    plt.plot(np.nanmean(data1['scaledCv']['mean'][::6,:],0),data1['scaledZ'], color = 'darkblue', linewidth = 2, label = label1)
    plt.plot(np.nanmean(data2['scaledCv']['mean'][::6,:],0),data2['scaledZ'], color = 'mediumseagreen', linewidth = 2, label = label2)
    plt.plot(np.nanmean(data3['scaledCv']['mean'][::6,:],0),data3['scaledZ'], color = 'gold', linewidth = 2, label = label3)
    plt.legend()
    plt.show()

def plot_scaledBL_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, var):

    ###################################
    ## PLOT TIMESERIES OF BL DEPTH
    ###################################

    print ('******')
    print ('')
    print ('Scaling cloudnet data by identified thetaE inversions:')
    print ('')

    # UM -> IFS comparisons:
    # 5. bl_depth -> sfc_bl_height

    ### set diagnostic naming flags for if IFS being used
    if np.logical_or(out_dir3 == 'OUT_25H/', out_dir3 == 'ECMWF_IFS/'):
        ifs_flag = True
    else:
        ifs_flag = False

    #################################################################
    ### interpolate obs cloudnet data for continuous array
    #################################################################
    obs_data = interpCloudnet(obs_data, month_flag, missing_files, doy)

    #### ---------------------------------------------------------------
    #### prepare cloudnet data
    #### ---------------------------------------------------------------

    #### set flagged data to nans
    if var == 'Cv':
        obs_data['Cv'][obs_data['Cv'] < 0.0] = np.nan
        um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
        ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
        misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan
    elif var == 'lwc':
        obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
        # obs_data['lwc'][obs_data['lwc'] == 0] = np.nan
        um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
        ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
        ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
        misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
        #### change units to g/m3
        obs_data['lwc'] = obs_data['lwc'] * 1e3
        um_data['model_lwc'] = um_data['model_lwc'] * 1e3
        misc_data['model_lwc'] = misc_data['model_lwc'] * 1e3
        ifs_data['model_lwc'] = ifs_data['model_lwc'] * 1e3
    elif var == 'iwc':
        obs_data['iwc'][obs_data['iwc'] == -999] = 0.0
        # obs_data['iwc'][obs_data['iwc'] == 0] = np.nan
        um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
        ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
        # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
        misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
        #### change units to g/m3
        obs_data['iwc'] = obs_data['iwc'] * 1e3
        um_data['model_iwc_filtered'] = um_data['model_iwc_filtered'] * 1e3
        misc_data['model_iwc_filtered'] = misc_data['model_iwc_filtered'] * 1e3
        ifs_data['model_snow_iwc_filtered'] = ifs_data['model_snow_iwc_filtered'] * 1e3

    # #### ---------------------------------------------------------------
    # #### ONLY LOOK AT SONDES FROM THE DRIFT
    # #### ---------------------------------------------------------------
    # drift = np.where(np.logical_and(obs['inversions']['thetaE']['time'] >= 225.9, obs['inversions']['thetaE']['time'] <= 258.0))

    #### ------------------------------------------------------------------------------
    #### load inversions data from RADIOSONDES (i.e. 6 hourly data)
    #### ------------------------------------------------------------------------------
    obsinv = obs['inversions']['thetaE']['invbase']
    obsmlh = obs['inversions']['thetaE']['sfmlheight']

    #### ------------------------------------------------------------------------------
    #### need to identify what cloudnet indices correspond to radiosondes
    #### ------------------------------------------------------------------------------
    missing_files = [225, 230, 253, 257]    # manually set missing files doy for now

    #### remove DOY 230, 253, 257 manually for now
    nanindices = np.array([16,17,18,19,108,109,110,111,124,125,126,127])

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    temp_time = obs['inversions']['thetaE']['time']         #### temporary time array for indexing
    temp_time2 = np.zeros(len(temp_time))
    temp_time2[:] = np.nan
    temp_time2[:nanindices[0]] = temp_time[:nanindices[0]]
    temp_time2[nanindices[3]+1:nanindices[4]] = temp_time[nanindices[3]+1:nanindices[4]]
    temp_time2[nanindices[7]+1:nanindices[8]] = temp_time[nanindices[7]+1:nanindices[8]]
    temp_inv = np.zeros(len(obsinv))
    temp_inv[:] = np.nan
    temp_inv[:nanindices[0]] = obsinv[:nanindices[0]]
    temp_inv[nanindices[3]+1:nanindices[4]] = obsinv[nanindices[3]+1:nanindices[4]]
    temp_inv[nanindices[7]+1:nanindices[8]] = obsinv[nanindices[7]+1:nanindices[8]]
    temp_sfml = np.zeros(len(obsmlh))
    temp_sfml[:] = np.nan
    temp_sfml[:nanindices[0]] = obsmlh[:nanindices[0]]
    temp_sfml[nanindices[3]+1:nanindices[4]] = obsmlh[nanindices[3]+1:nanindices[4]]
    temp_sfml[nanindices[7]+1:nanindices[8]] = obsmlh[nanindices[7]+1:nanindices[8]]

    ### reset time array with new one accounting for missing FILES
    obs['inversions']['thetaE']['time'] = temp_time2
    obsinv = temp_inv
    obsmlh = temp_sfml

    ### save non-nan values to dictionary
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    obs['inversions']['TimesForCloudnet'] = obs['inversions']['thetaE']['time'][~np.isnan(obs['inversions']['thetaE']['time'])]
    obs['inversions']['InvBasesForCloudnet'] = obsinv[~np.isnan(obsinv)]
    obs['inversions']['sfmlForCloudnet'] = obsmlh[~np.isnan(obsmlh)]

    #### ------------------------------------------------------------------------------
    #### fill obs array with Cloudnet height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    ### define scaledZ array to sort data in to
    ###     will act as mid point of vertical "boxes" of width 0.1
    binres = 0.1
    Zpts = np.arange(0.0 + binres/2.0, 1.0 + binres/2.0, binres)

    ### use 6hourly cloudnet data to compare radiosonde inversion heights to
    obs_data['height_6hrly'] = obs_data['height'][::6,:]
    obs_data[var + '_6hrly'] = obs_data[var][::6,:]
    obs_data['time_6hrly'] = obs_data['time'][::6]      ### 6 hourly cloudnet data

    ### initialise array to hold height indices, set all to nan before filling
    obsind = np.zeros(np.size(obs['inversions']['InvBasesForCloudnet'])); #obsind[:] = np.nan
    obssfml = np.zeros(np.size(obs['inversions']['sfmlForCloudnet'])); #obssfml[:] = np.nan

    ### look for altitudes <= invbase in obs cloudnet data
    for i in range(0, np.size(obs['inversions']['TimesForCloudnet'])):        ### time loop (from radiosondes)
        #### check if there are any height levels below the inversion
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])) > 0.0:
            ### if there is, set obsind to the last level before the inversion
            obsind[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['InvBasesForCloudnet'][i])[0][-1]
        #### check if there are any height levels below the sfmlheight
        if np.size(np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])) > 0.0:
            ### if there is, set obssfml to the last level before the sfmlheight
            obssfml[i] = np.where(obs_data['height_6hrly'][i,:] <= obs['inversions']['sfmlForCloudnet'][i])[0][-1]

    fig = plt.figure(figsize=(8,6))
    plt.title('temp fig: radiosonde invbase w/pulled cloudnet inv height')
    for i in range(0, np.size(obsind)): plt.plot(obs_data['time_6hrly'][i],obs_data['height_6hrly'][i,int(obsind[i])],'o')
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),obs['inversions']['thetaE']['invbase'])
    plt.xlabel('DOY')
    plt.ylabel('Z [m]')
    # plt.savefig('FIGS/' + var + '_obs_inversionDetection_timeseries.png')
    plt.show()

    ### save inversion base index into dictionary
    obs['inversions']['invbase_kIndex'] = obsind
    obs['inversions']['sfmlheight_kIndex'] = obssfml

    ### initialise scaled arrays in dictionary
    obs['inversions']['scaled' + var] = {}
    obs['inversions']['scaled' + var]['binned'] = {}
    obs['inversions']['scaled' + var]['mean'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaled' + var]['mean'][:] = np.nan
    obs['inversions']['scaled' + var]['stdev'] = np.zeros([np.size(obs_data['time_6hrly']),len(Zpts)]); obs['inversions']['scaled' + var]['stdev'][:] = np.nan
    obs['inversions']['scaledZ'] = Zpts
    obs['inversions']['scaledTime'] = obs_data['time_6hrly']
    obs['inversions']['bl' + var] = np.zeros([np.size(obs_data['height_6hrly'],0),np.size(obs_data['height_6hrly'],1)]); obs['inversions']['bl' + var][:] = np.nan

    ### fill arrays with cloudnet data below each invbase
    for i in range(0,np.size(obs['inversions']['TimesForCloudnet'])):     ## loop over radiosonde time
        print(str(i) + 'th timestep (radiosonde):')

        ### create new dictionary entry for i-th timestep
        obs['inversions']['scaled' + var]['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        ###         +1 includes invbase in array
        if obs['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts = obs_data['height_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts = hgts / obs_data['height_6hrly'][i,int(obs['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        obs['inversions']['bl' + var][i,:int(obs['inversions']['invbase_kIndex'][i]+1)] = obs_data[var + '_6hrly'][i,:int(obs['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar = np.where(np.logical_and(scaled_hgts >= Zpts[k] - binres/2.0, scaled_hgts < Zpts[k] + binres/2.0))
            obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = obs['inversions']['bl' + var][i,tempvar]
            if np.size(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                obs['inversions']['scaled' + var]['mean'][i,k] = np.nanmean(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            obs['inversions']['scaled' + var]['stdev'][i,k] = np.nanstd(obs['inversions']['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

    #### ---------------------------------------------------------------
    #### prepare model inversion data
    ####        data from thetaE algorithm is on radiosonde (6 hourly) timesteps already
    #### ---------------------------------------------------------------

    ### need to make sure times where we don't have an inversion (invbase == nan) in the IFS doesn't mess up the algorithm
    ###     set these cases to zero for now
    data3['inversions']['invbase'][np.isnan(data3['inversions']['invbase'])] = 0.0

    ### need to build new arrays manually, isn't allowing indexing + ==nan for some reason...
    ####        time
    tim1 = np.zeros(len(data1['inversions']['time']))
    tim1[:] = np.nan
    tim1[:nanindices[0]] = data1['inversions']['time'][:nanindices[0]]
    tim1[nanindices[3]+1:nanindices[4]] = data1['inversions']['time'][nanindices[3]+1:nanindices[4]]
    tim1[nanindices[7]+1:nanindices[8]] = data1['inversions']['time'][nanindices[7]+1:nanindices[8]]
    tim2 = tim1
    tim3 = tim1
    ####        inversions
    inv1 = np.zeros(len(data1['inversions']['invbase']))
    inv1[:] = np.nan
    inv1[:nanindices[0]] = data1['inversions']['invbase'][:nanindices[0]]
    inv1[nanindices[3]+1:nanindices[4]] = data1['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv1[nanindices[7]+1:nanindices[8]] = data1['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv2 = np.zeros(len(data2['inversions']['invbase']))
    inv2[:] = np.nan
    inv2[:nanindices[0]] = data2['inversions']['invbase'][:nanindices[0]]
    inv2[nanindices[3]+1:nanindices[4]] = data2['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv2[nanindices[7]+1:nanindices[8]] = data2['inversions']['invbase'][nanindices[7]+1:nanindices[8]]
    inv3 = np.zeros(len(data3['inversions']['invbase']))
    inv3[:] = np.nan
    inv3[:nanindices[0]] = data3['inversions']['invbase'][:nanindices[0]]
    inv3[nanindices[3]+1:nanindices[4]] = data3['inversions']['invbase'][nanindices[3]+1:nanindices[4]]
    inv3[nanindices[7]+1:nanindices[8]] = data3['inversions']['invbase'][nanindices[7]+1:nanindices[8]]

    ### save non-nan values
    ###         these arrays will be used for picking out inversions in the cloudnet files
    ###             (and so miss out dates where we're missing cloudnet data)
    tim1 = tim1[~np.isnan(tim1)]
    tim2 = tim2[~np.isnan(tim2)]
    tim3 = tim3[~np.isnan(tim3)]
    inv1 = inv1[~np.isnan(inv1)]
    inv2 = inv2[~np.isnan(inv2)]
    inv3 = inv3[~np.isnan(inv3)]

    #### calculate inversion algorithm success rate
    ind1 = np.where(inv1 >= 0.0)  ## non-nan values
    data1['inversions']['successRate'] = np.size(ind1[0]) / np.float(np.size(inv1)) * 100.0
    ind2 = np.where(inv2 >= 0.0)  ## non-nan values
    data2['inversions']['successRate'] = np.size(ind2[0]) / np.float(np.size(inv2)) * 100.0
    ind3 = np.where(inv3 >= 0.0)  ## non-nan values
    data3['inversions']['successRate'] = np.size(ind3[0]) / np.float(np.size(inv3)) * 100.0

    print (label1 + ' inversion algorithm success rate = ' + str(data1['inversions']['successRate']))
    print (label2 + ' inversion algorithm success rate = ' + str(data2['inversions']['successRate']))
    print (label3 + ' inversion algorithm success rate = ' + str(data3['inversions']['successRate']))
    print ('****')

    # #### ---------------------------------------------------------------
    # #### remove flagged IFS heights
    # #### ---------------------------------------------------------------
    # data3['height'][data3['height'] == -9999] = 0.0
    #         #### set all heights to zero if flagged. setting to nan caused problems further on
    # data3['height_hrly'] = np.squeeze(data3['height'][data3['hrly_flag'],:])  ### need to explicitly save since height coord changes at each timedump

    # #### ---------------------------------------------------------------
    # #### Define meteorological periods from Jutta's paper
    # #### ---------------------------------------------------------------
    #
    # ## Meteorological period definitions from Jutta's paper:
    # ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    # ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    # ##     radar measurements were not possible during heavy ice breaking because of excessive
    # ##     vibration, cloud characteristics and fog heights are not available during period 2.
    # ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    # ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    # ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    # ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    # ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    # ##     out to the ice edge. "
    # #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...
    #
    # ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    # p3 = np.where(um_data['time'] < 230.0)
    # p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    # p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    # p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    # p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    # p8 = np.where(um_data['time'] >= 255.5)

    #### ---------------------------------------------------------------
    #### Use extracted height indices to probe cloudnet data
    #### ---------------------------------------------------------------

    ### save new height and cloudnet time array into dictionary (latter to account for missing files)
    data1['scaledZ'] = Zpts
    data1['scaledTime'] = um_data['time'].data[::6]
    data2['scaledZ'] = Zpts
    data2['scaledTime'] = misc_data['time'].data[::6]
    data3['scaledZ'] = Zpts
    data3['scaledTime'] = ifs_data['time'].data[::6]

    #### define empty arrays of nans to fill with scaled data
    data1['scaled' + var] = {}
    data1['scaled' + var]['binned'] = {}
    data1['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaled' + var]['mean'][:] = np.nan
    data1['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data1['scaled' + var]['stdev'][:] = np.nan
    data1['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(um_data['height'],1)]); data1['bl' + var][:] = np.nan

    data2['scaled' + var] = {}
    data2['scaled' + var]['binned'] = {}
    data2['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaled' + var]['mean'][:] = np.nan
    data2['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data2['scaled' + var]['stdev'][:] = np.nan
    data2['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(misc_data['height'],1)]); data2['bl' + var][:] = np.nan

    data3['scaled' + var] = {}
    data3['scaled' + var]['binned'] = {}
    data3['scaled' + var]['mean'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaled' + var]['mean'][:] = np.nan
    data3['scaled' + var]['stdev'] = np.zeros([np.size(data1['scaledTime']),len(Zpts)]); data3['scaled' + var]['stdev'][:] = np.nan
    data3['bl' + var] = np.zeros([np.size(data1['scaledTime']),np.size(ifs_data['height'],1)]); data3['bl' + var][:] = np.nan

    ### ------------------------------------------------------------------------------------------
    ### find cloudnet timesteps which match the inversion timesteps
    ### ------------------------------------------------------------------------------------------
    data1['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data1['scaled' + var]['inversion_Tindex'][:] = np.nan
    data1['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data1['scaled' + var]['inversionForCloudnet'][:] = np.nan
    data2['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data2['scaled' + var]['inversion_Tindex'][:] = np.nan
    data2['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data2['scaled' + var]['inversionForCloudnet'][:] = np.nan
    data3['scaled' + var]['inversion_Tindex'] = np.zeros(np.size(data1['scaledTime'])); data3['scaled' + var]['inversion_Tindex'][:] = np.nan
    data3['scaled' + var]['inversionForCloudnet'] = np.zeros(np.size(data1['scaledTime'])); data3['scaled' + var]['inversionForCloudnet'][:] = np.nan

    for i in range(0, len(tim1)):
        ## find the cloudnet time INDEX which matches the inversion timestep
        if np.size(np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data1['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data1['scaledTime'],3) == np.round(tim1[i],3))[0][0]
        if np.size(np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data2['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data2['scaledTime'],3) == np.round(tim2[i],3))[0][0]
        if np.size(np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))) > 0:
            ### gives INDEX of CLOUDNET DATA corresponding to that timestep
            data3['scaled' + var]['inversion_Tindex'][i] = np.where(np.round(data3['scaledTime'],3) == np.round(tim3[i],3))[0][0]


        ### use inversion_Tindices to define new inversion height array on cloudnet timesteps for looping over
        ###         if inversion_Tindex is not NaN, use to index inv into new array (cninv)
        if data1['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data1['scaled' + var]['inversionForCloudnet'][int(data1['scaled' + var]['inversion_Tindex'][i])] = inv1[i]
        if data2['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data2['scaled' + var]['inversionForCloudnet'][int(data2['scaled' + var]['inversion_Tindex'][i])] = inv2[i]
        if data3['scaled' + var]['inversion_Tindex'][i] >= 0.0:
            data3['scaled' + var]['inversionForCloudnet'][int(data3['scaled' + var]['inversion_Tindex'][i])] = inv3[i]
    #
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    #
    #### ---------------------------------------------------------------
    #### Look at data below main inversion base only - model data
    #### ---------------------------------------------------------------
    #### create empty arrays to hold height index
    zind1 = np.zeros(np.size(data1['scaledTime'])); zind1[:] = np.nan
    zind2 = np.zeros(np.size(data2['scaledTime'])); zind2[:] = np.nan
    zind3 = np.zeros(np.size(data3['scaledTime'])); zind3[:] = np.nan

    #### ------------------------------------------------------------------------------
    #### fill model arrays with height index of main inversion base / sfml height
    #### ------------------------------------------------------------------------------
    for i in range(0, np.size(data1['scaledTime'])):        ### all can go in this loop, data1['scaledTime'] == 6-hourly data

        ### main inversion base assignments
        ###         (1) find where UM height array matches the invbase index for index i
        ###         (2) find where UM height array matches the invbase index for index i
        ###         (3) find where IFS height array is less than or equal to the UM-gridded invbase index for index i
        if np.size(np.where(um_data['height'][i,:].data == data1['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            zind1[i] = np.where(um_data['height'][i,:].data == data1['scaled' + var]['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(misc_data['height'][i,:].data == data2['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            zind2[i] = np.where(misc_data['height'][i,:].data == data2['scaled' + var]['inversionForCloudnet'][i])[0][0]
        if np.size(np.where(ifs_data['height'][i,:].data <= data3['scaled' + var]['inversionForCloudnet'][i])) > 0.0:
            temp = ifs_data['height'][i,:].data <= data3['scaled' + var]['inversionForCloudnet'][i]
            zind3[i] = np.where(temp == True)[0][-1]

    #### assign height indices to dictionary for later use
    data1['inversions']['invbase_kIndex'] = zind1
    data2['inversions']['invbase_kIndex'] = zind2
    data3['inversions']['invbase_kIndex'] = zind3

    fig = plt.figure(figsize=(9,10))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(311)
    plt.title(label1)
    for i in range(0, np.size(zind1)):
        if ~np.isnan(zind1[i]): plt.plot(data1['scaledTime'][i],um_data['height'][i,int(zind1[i])],'o')
    plt.plot(tim1,inv1)
    plt.ylabel('Z [m]')
    plt.ylim([0,3e3])
    plt.subplot(312)
    plt.title(label2)
    for i in range(0, np.size(zind2)):
        if ~np.isnan(zind2[i]): plt.plot(data2['scaledTime'][i],misc_data['height'][i,int(zind2[i])],'o')
    plt.plot(tim2, inv2)
    plt.ylim([0,3e3])
    plt.ylabel('Z [m]')
    plt.subplot(313)
    plt.title(label3)
    for i in range(0, np.size(zind3)):
        if ~np.isnan(zind3[i]): plt.plot(data3['scaledTime'][i],ifs_data['height'][i,int(zind3[i])],'o')
    plt.plot(tim3, inv3)
    plt.ylim([0,3e3])
    plt.ylabel('Z [m]')
    plt.xlabel('DOY')
    # plt.savefig('FIGS/' + var + '_model_inversionDetection_timeseries.png')
    plt.show()

    ### set 6 hourly cloudnet Cv arrays as tempvars
    if var == 'Cv':
        ra2m_var = um_data['model_Cv_filtered'][::6,:]
        casim_var = misc_data['model_Cv_filtered'][::6,:]
        ifs_var = ifs_data['model_snow_Cv_filtered'][::6,:]
    elif var == 'lwc':
        ra2m_var = um_data['model_lwc'][::6,:]
        casim_var = misc_data['model_lwc'][::6,:]
        ifs_var = ifs_data['model_lwc'][::6,:]
    elif var == 'iwc':
        ra2m_var = um_data['model_iwc_filtered'][::6,:]
        casim_var = misc_data['model_iwc_filtered'][::6,:]
        ifs_var = ifs_data['model_snow_iwc_filtered'][::6,:]

    ### find all Cv data below identified inversion
    for i in range(0,np.size(data1['scaledTime'])):     ## loop over time
        print ()
        print(str(i) + 'th timestep (model data):')

        ### create new dictionary entry for i-th timestep
        data1['scaled' + var]['binned']['t' + str(i)] = {}
        data2['scaled' + var]['binned']['t' + str(i)] = {}
        data3['scaled' + var]['binned']['t' + str(i)] = {}

        ###-----------------------------------------------------------------------------------------
        ### for main inversion
        ###-----------------------------------------------------------------------------------------
        ### create array of height points under the identified inversion
        if data1['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts1 = um_data['height'][i,:int(data1['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data2['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts2 = misc_data['height'][i,:int(data2['inversions']['invbase_kIndex'][i])]
        else:
            continue
        if data3['inversions']['invbase_kIndex'][i] >= 0.0:
            hgts3 = ifs_data['height'][i,:int(data3['inversions']['invbase_kIndex'][i])]
        else:
            continue

        ### scale BL height array by the inversion depth to give range Z 0 to 1 (1 = inversion height) (temporary variable)
        scaled_hgts1 = hgts1 / um_data['height'][i,int(data1['inversions']['invbase_kIndex'][i])]
        scaled_hgts2 = hgts2 / misc_data['height'][i,int(data2['inversions']['invbase_kIndex'][i])]
        scaled_hgts3 = hgts3 / ifs_data['height'][i,int(data3['inversions']['invbase_kIndex'][i])]

        # find Cv values below the BL inversion
        data1['bl' + var][i,:int(data1['inversions']['invbase_kIndex'][i]+1)] = ra2m_var[i,:int(data1['inversions']['invbase_kIndex'][i]+1)]
        data2['bl' + var][i,:int(data2['inversions']['invbase_kIndex'][i]+1)] = casim_var[i,:int(data2['inversions']['invbase_kIndex'][i]+1)]
        data3['bl' + var][i,:int(data3['inversions']['invbase_kIndex'][i]+1)] = ifs_var[i,:int(data3['inversions']['invbase_kIndex'][i]+1)]

        ## bin scaled BL heights into pre-set Zpts array so every timestep can be compared
        for k in range(len(Zpts)):
            tempvar1 = np.where(np.logical_and(scaled_hgts1 >= Zpts[k] - binres/2.0, scaled_hgts1 < Zpts[k] + binres/2.0))
            tempvar2 = np.where(np.logical_and(scaled_hgts2 >= Zpts[k] - binres/2.0, scaled_hgts2 < Zpts[k] + binres/2.0))
            tempvar3 = np.where(np.logical_and(scaled_hgts3 >= Zpts[k] - binres/2.0, scaled_hgts3 < Zpts[k] + binres/2.0))

            ### find mean and stdev of points within given Zpts range
            data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data1['bl' + var][i,tempvar1]
            if np.size(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data1['scaled' + var]['mean'][i,k] = np.nanmean(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data1['scaled' + var]['stdev'][i,k] = np.nanstd(data1['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

            data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data2['bl' + var][i,tempvar2]
            if np.size(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data2['scaled' + var]['mean'][i,k] = np.nanmean(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data2['scaled' + var]['stdev'][i,k] = np.nanstd(data2['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

            data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]] = data3['bl' + var][i,tempvar3]
            if np.size(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]]) > 0:
                data3['scaled' + var]['mean'][i,k] = np.nanmean(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]])
            data3['scaled' + var]['stdev'][i,k] = np.nanstd(data3['scaled' + var]['binned']['t' + str(i)][Zpts[k]])

    ### save working data for debug
    # np.save('working_data1', data1)

    ##################################################
    ##################################################
    #### figures
    ##################################################
    ##################################################

    ### timeseries
    fig = plt.figure(figsize=(8,12))
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.98, left = 0.1,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(411)
    plt.title('Obs')
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(412)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(413)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.subplot(414)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaled' + var]['mean']), vmin = 0, vmax = 1)
    plt.ylabel('Z [scaled]')
    plt.xlabel('DOY')
    # plt.savefig('FIGS/' + var + '_ALL_scaledZ_timeseries.png')
    plt.show()

    ### obs
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title('Obs - 6hourly because inversions from radiosondes')
    plt.pcolor(obs_data['time_6hrly'].data,obs_data['height_6hrly'][0,:].data,np.transpose(obs['inversions']['bl' + var])); plt.ylim([0,3e3])
    plt.plot(np.squeeze(obs['inversions']['thetaE']['time']),np.squeeze(obs['inversions']['thetaE']['invbase']),'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(obs['inversions']['scaledTime'],obs['inversions']['scaledZ'],np.transpose(obs['inversions']['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_obs_scaledZ_timeseries.png')
    plt.show()

    ### um_ra2m
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label1)
    plt.pcolor(data1['scaledTime'],um_data['height'][0,:].data,np.transpose(data1['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data1['inversions']['time'],data1['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data1['scaledTime'],data1['scaledZ'],np.transpose(data1['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    # plt.savefig('FIGS/' + var + '_RA2M_scaledZ_timeseries.png')
    plt.show()

    ### um_casim-100
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label2)
    plt.pcolor(data2['scaledTime'],misc_data['height'][0,:].data,np.transpose(data2['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data2['inversions']['time'],data2['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data2['scaledTime'],data2['scaledZ'],np.transpose(data2['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    plt.savefig('FIGS/' + var + '_CASIM-100_scaledZ_timeseries.png')
    plt.show()

    ### ecmwf_ifs
    fig = plt.figure(figsize=(8,6))
    plt.subplot(211)
    plt.title(label3)
    plt.pcolor(data3['scaledTime'],ifs_data['height'][0,:].data,np.transpose(data3['bl' + var])); plt.ylim([0,3e3])
    plt.plot(data3['inversions']['time'],data3['inversions']['invbase'],'r')
    plt.xlim([226,258])
    plt.subplot(212)
    plt.pcolor(data3['scaledTime'],data3['scaledZ'],np.transpose(data3['scaled' + var]['mean'])); plt.ylim([0,1])
    plt.xlim([226,258])
    # plt.savefig('FIGS/' + var + '_IFS_scaledZ_timeseries.png')
    plt.show()

    ####================================================================
    ### profiles
    ###         loop through and set all zeros to nans
    obsmean = obs['inversions']['scaled' + var]['mean']
    ra2mmean = data1['scaled' + var]['mean']
    casimmean = data2['scaled' + var]['mean']
    ifsmean = data3['scaled' + var]['mean']
    # for i in range(0, len(data1['scaledTime'])):
    #     obsmean[i,obs['inversions']['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     ra2mmean[i,data1['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     casimmean[i,data2['scaled' + var]['mean'][i,:] == 0.0] = np.nan
    #     ifsmean[i,data3['scaled' + var]['mean'][i,:] == 0.0] = np.nan

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
    ### Build figure (profiles)
    ### -------------------------------
    fig = plt.figure(figsize=(5,6))
    ax1  = fig.add_axes([0.2,0.1,0.7,0.8])   # left, bottom, width, height
    plt.plot(np.nanmean(obsmean,0),obs['inversions']['scaledZ'], '--', color = 'k', linewidth = 2, label = 'Obs')
    ax1.fill_betweenx(data1['scaledZ'],np.nanmean(obsmean,0) - np.nanstd(obsmean,0),
        np.nanmean(obsmean,0) + np.nanstd(obsmean,0), color = 'lightgrey', alpha = 0.5)
    # ax1.fill_betweenx(data1['scaledZ'],np.nanmean(obsmean,0) - np.nanstd(obs['inversions']['scaled' + var]['stdev'],0),
    #     np.nanmean(obsmean,0) + np.nanstd(obs['inversions']['scaled' + var]['stdev'],0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(ra2mmean,0),data1['scaledZ'], '^-',
        color = 'darkblue', markeredgecolor = 'midnightblue', linewidth = 2, label = label1)
    ax1.fill_betweenx(data1['scaledZ'],np.nanmean(ra2mmean,0) - np.nanstd(ra2mmean,0),
        np.nanmean(ra2mmean,0) + np.nanstd(ra2mmean,0), color = 'lightblue', alpha = 0.4)
    # ax1.fill_betweenx(data1['scaledZ'],np.nanmean(ra2mmean,0) - np.nanstd(data1['scaled' + var]['stdev'],0),
    #     np.nanmean(ra2mmean,0) + np.nanstd(data1['scaled' + var]['stdev'],0), color = 'lightblue', alpha = 0.4)
    plt.plot(np.nanmean(ifsmean,0),data3['scaledZ'], 'd-',
        color = 'gold', markeredgecolor = 'saddlebrown', linewidth = 2, label = label3)
    ax1.fill_betweenx(data3['scaledZ'],np.nanmean(ifsmean,0) - np.nanstd(ifsmean,0),
        np.nanmean(ifsmean,0) + np.nanstd(ifsmean,0), color = 'navajowhite', alpha = 0.35)
    # ax1.fill_betweenx(data3['scaledZ'],np.nanmean(ifsmean,0) - np.nanstd(data3['scaled' + var]['stdev'],0),
    #     np.nanmean(ifsmean,0) + np.nanstd(data3['scaled' + var]['stdev'],0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(casimmean,0),data2['scaledZ'], 'v-',
        color = 'mediumseagreen', markeredgecolor = 'darkslategrey', linewidth = 2, label = label2)
    ax1.fill_betweenx(data2['scaledZ'],np.nanmean(casimmean,0) - np.nanstd(casimmean,0),
        np.nanmean(casimmean,0) + np.nanstd(casimmean,0), color = 'mediumaquamarine', alpha = 0.15)
    # ax1.fill_betweenx(data2['scaledZ'],np.nanmean(casimmean,0) - np.nanstd(data2['scaled' + var]['stdev'],0),
    #     np.nanmean(casimmean,0) + np.nanstd(data2['scaled' + var]['stdev'],0), color = 'mediumaquamarine', alpha = 0.15)
    if var == 'Cv': plt.xlim([0,1])
    if var == 'lwc': plt.xlim([0,0.2])
    if var == 'iwc': plt.xlim([0,0.02])
    plt.ylim([0,1])
    if var == 'Cv': plt.xlabel('Cv')
    plt.ylabel('scaled Z \n (0 = lowest level; 1 = inversion base height)')
    plt.legend()
    plt.savefig('FIGS/' + var + '_Obs-IFSgrid-QF10_scaledZ.svg')
    plt.show()

def period_Selection(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, um_out_dir, doy, obs_switch, obs, data1, data2, data3, data4, nanind, wcind):

    # #### ---------------------------------------------------------------
    # #### Define meteorological periods from Jutta's paper
    # #### ---------------------------------------------------------------
    #
    # ## Meteorological period definitions from Jutta's paper:
    # ##     "Period 1 covers the time in the MIZ until 4 August 06:00 UTC. Period 2 encompasses
    # ##     the journey into the ice towards the North Pole until 12 August 00:00 UTC. Since cloud
    # ##     radar measurements were not possible during heavy ice breaking because of excessive
    # ##     vibration, cloud characteristics and fog heights are not available during period 2.
    # ##     Period 3 (12 to 17 August) includes the 'North Pole' station and the beginning of
    # ##     the ice drift. Period 4 (18 to 27 August) covers the end of the melt and the transition
    # ##     period into the freeze up. The freeze up is covered by period 5 (28 August to 3 September),
    # ##     6 (4 to 7 September) and 7 (8 to 12 September 12:00 UTC). Finally, period 8 (12 September
    # ##     12:00 UTC to 21 September 06:00 UTC) covers the end of the ice drift period and the transit
    # ##     out to the ice edge. "
    # #######         from Jutta: so if there is no time stamp mentioned it is eg. P4 18.08 0000UTC - 27.08 23:59 UTC , P5 then is 28.08 00UTC until 03.09 23:59 UTC...
    #
    # ### define met periods wrt cloudnet timestamps for ease (all runs should be able to use same indexing)
    # p3 = np.where(um_data['time'] < 230.0)
    # p4 = np.where(np.logical_and(um_data['time'] >= 230.0, um_data['time'] < 240.0))
    # p5 = np.where(np.logical_and(um_data['time'] >= 240.0, um_data['time'] < 247.0))
    # p6 = np.where(np.logical_and(um_data['time'] >= 247.0, um_data['time'] < 251.0))
    # p7 = np.where(np.logical_and(um_data['time'] >= 251.0, um_data['time'] < 255.5))
    # p8 = np.where(um_data['time'] >= 255.5)

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
        # from matplotlib.patches import Polygon

    ###################################
    ## PLOT MAP
    ###################################

    print ('******')
    print ('')
    print ('Plotting statistics based on periods 3,4 and 5,6:')
    print ('')

    # print um_data.keys()

    ###----------------------------------------------------------------
    ###         Set flagged Cvs to nans
    ###----------------------------------------------------------------
    # obs_data['Cv'][obs_data['Cv'] == -999] = np.nan
    # um_data['model_Cv_filtered'][um_data['model_Cv_filtered'] < 0.0] = np.nan
    # ifs_data['model_snow_Cv_filtered'][ifs_data['model_snow_Cv_filtered'] < 0.0] = np.nan
    # misc_data['model_Cv_filtered'][misc_data['model_Cv_filtered'] < 0.0] = np.nan
    # ra2t_data['model_Cv_filtered'][ra2t_data['model_Cv_filtered'] < 0.0] = np.nan

    ###----------------------------------------------------------------
    ###         1) LWC/IWC Method
    ###----------------------------------------------------------------
    # #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] < 0] = 0.0
    # um_data['model_lwc'][um_data['model_lwc'] < 0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 0] = 0.0
    # misc_data['model_lwc'][misc_data['model_lwc'] < 0] = 0.0
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0] = 0.0
    #
    # obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0] = 0.0
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0] = 0.0
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0] = 0.0
    #
    # # ###----------------------------------------------------------------
    # # ###         Calculate total water content
    # # ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    # um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    # misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    # ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    # ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']
    #
    # #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] < 1e-6] = np.nan
    # um_data['model_lwc'][um_data['model_lwc'] < 1e-6] = np.nan
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 1e-6] = np.nan
    # ifs_data['model_lwc'][ifs_data['model_lwc'] >= 20.0] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 1e-6] = np.nan
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 1e-6] = np.nan
    #
    # #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] < 1e-6] = np.nan
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 1e-6] = np.nan
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 1e-6] = np.nan
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 20.0] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 1e-6] = np.nan
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 1e-6] = np.nan
    #
    # #### set flagged um_data to nans
    # obs_data['twc'][obs_data['twc'] < 1e-6] = np.nan
    # um_data['model_twc'][um_data['model_twc'] < 1e-6] = np.nan
    # ifs_data['model_twc'][ifs_data['model_twc'] < 1e-6] = np.nan
    # ifs_data['model_twc'][ifs_data['model_twc'] >= 20.0] = np.nan
    # misc_data['model_twc'][misc_data['model_twc'] < 1e-6] = np.nan
    # ra2t_data['model_twc'][ra2t_data['model_twc'] < 1e-6] = np.nan

    ###----------------------------------------------------------------
    ###         2) TWC Method
    ###----------------------------------------------------------------
    #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 5.0e-3] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] == -999] = 0.0
    # # obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] == -999] = 0.0
    # um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### initialise TWC masks
    mask0 = np.zeros([np.size(obs_data['twc'],0), np.size(obs_data['twc'],1)])
    mask1 = np.zeros([np.size(um_data['model_twc'],0), np.size(um_data['model_twc'],1)])
    mask2 = np.zeros([np.size(misc_data['model_twc'],0), np.size(misc_data['model_twc'],1)])
    mask3 = np.zeros([np.size(ifs_data['model_twc'],0), np.size(ifs_data['model_twc'],1)])
    mask4 = np.zeros([np.size(ra2t_data['model_twc'],0), np.size(ra2t_data['model_twc'],1)])

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                    obs_data['iwc'][t,k] = np.nan
                    obs_data['lwc_adiabatic'][t,k] = np.nan
                    obs_data['lwc'][t,k] = np.nan
                elif obs_data['twc'][t,k] >= twc_thresh_um[k]:
                    mask0[t,k] = 1.0
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
                um_data['model_iwc_filtered'][t,k] = np.nan
                um_data['model_lwc'][t,k] = np.nan
            elif um_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask1[t,k] = 1.0
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
                misc_data['model_iwc_filtered'][t,k] = np.nan
                misc_data['model_lwc'][t,k] = np.nan
            elif misc_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask2[t,k] = 1.0
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
                ra2t_data['model_iwc_filtered'][t,k] = np.nan
                ra2t_data['model_lwc'][t,k] = np.nan
            elif ra2t_data['model_twc'][t,k] >= twc_thresh_um[k]:
                mask4[t,k] = 1.0
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
                ifs_data['model_snow_iwc_filtered'][t,k] = np.nan
                ifs_data['model_lwc'][t,k] = np.nan
            elif ifs_data['model_twc'][t,k] >= twc_thresh_ifs[k]:
                mask3[t,k] = 1.0

    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    ####    DEFINE PERIODS
    #### ---------------------------------------------------------------------------------------------------
    #### ---------------------------------------------------------------------------------------------------
    p3 = np.where(np.logical_and(obs_data['time'] >= doy[0], obs_data['time'] < 230.0))
    p4 = np.where(np.logical_and(obs_data['time'] >= 230.0, obs_data['time'] < 240.0))
    p5 = np.where(np.logical_and(obs_data['time'] >= 240.0, obs_data['time'] < 247.0))
    p6 = np.where(np.logical_and(obs_data['time'] >= 247.0, obs_data['time'] < 251.0))

    mask0[nanind] = np.nan
    mask1[nanind] = np.nan
    mask2[nanind] = np.nan
    mask3[nanind] = np.nan
    mask4[nanind] = np.nan

    mask0[wcind] = np.nan
    mask1[wcind] = np.nan
    mask2[wcind] = np.nan
    mask3[wcind] = np.nan
    mask4[wcind] = np.nan

    ##################################################
    ##################################################
    #### 	P3 AND P4
    ##################################################
    ##################################################

    SMALL_SIZE = 12
    MED_SIZE = 15
    LARGE_SIZE = 16

    plt.rc('font',size=LARGE_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)
    plt.figure(figsize=(10,10))
    plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
            hspace = 0.22, wspace = 0.19)

    fraction0p3 = np.squeeze(mask0[p3,:]) # np.squeeze(obs_data['Cv'][p3,:]) #
    fraction1p3 = np.squeeze(mask1[p3,:]) # np.squeeze(um_data['model_Cv_filtered'][p3,:]) #
    fraction2p3 = np.squeeze(mask2[p3,:]) # np.squeeze(misc_data['model_Cv_filtered'][p3,:]) #
    fraction3p3 = np.squeeze(mask3[p3,:]) # np.squeeze(ifs_data['model_snow_Cv_filtered'][p3,:]) #
    fraction4p3 = np.squeeze(mask4[p3,:]) # np.squeeze(ra2t_data['model_Cv_filtered'][p3,:]) #

    fraction0p4 = np.squeeze(mask0[p4,:]) # np.squeeze(obs_data['Cv'][p4,:]) #
    fraction1p4 = np.squeeze(mask1[p4,:]) # np.squeeze(um_data['model_Cv_filtered'][p4,:]) #
    fraction2p4 = np.squeeze(mask2[p4,:]) # np.squeeze(misc_data['model_Cv_filtered'][p4,:]) #
    fraction3p4 = np.squeeze(mask3[p4,:]) # np.squeeze(ifs_data['model_snow_Cv_filtered'][p4,:]) #
    fraction4p4 = np.squeeze(mask4[p4,:]) # np.squeeze(ra2t_data['model_Cv_filtered'][p4,:]) #

    plt.subplot(231)
    ax2 = plt.gca()
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),np.nanmean(fraction0p3,0) - np.nanstd(fraction0p3,0),
        np.nanmean(fraction0p3,0) + np.nanstd(fraction0p3,0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(fraction0p3,0) - np.nanstd(fraction0p3,0), np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(fraction0p3,0) + np.nanstd(fraction0p3,0), np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p3,:]),0),np.nanmean(fraction1p3,0) - np.nanstd(fraction1p3,0),
        np.nanmean(fraction1p3,0) + np.nanstd(fraction1p3,0), color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(fraction1p3,0) - np.nanstd(fraction1p3,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(fraction1p3,0) + np.nanstd(fraction1p3,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0),np.nanmean(fraction4p3,0) - np.nanstd(fraction4p3,0),
        np.nanmean(fraction4p3,0) + np.nanstd(fraction4p3,0), color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(fraction4p3,0) - np.nanstd(fraction4p3,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(fraction4p3,0) + np.nanstd(fraction4p3,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p3,:]),0),np.nanmean(fraction2p3,0) - np.nanstd(fraction2p3,0),
        np.nanmean(fraction2p3,0) + np.nanstd(fraction2p3,0), color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(fraction2p3,0) - np.nanstd(fraction2p3,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(fraction2p3,0) + np.nanstd(fraction2p3,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),np.nanmean(fraction3p3,0) - np.nanstd(fraction3p3,0),
        np.nanmean(fraction3p3,0) + np.nanstd(fraction3p3,0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(fraction3p3,0) - np.nanstd(fraction3p3,0), np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(fraction3p3,0) + np.nanstd(fraction3p3,0), np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmean(fraction0p3,0),np.nanmean(np.squeeze(obs_data['height'][p3,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 5)
    plt.plot(np.nanmean(fraction3p3,0), np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 4)
    plt.plot(np.nanmean(fraction2p3,0),np.nanmean(np.squeeze(misc_data['height'][p3,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    plt.plot(np.nanmean(fraction4p3,0),np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 2)
    plt.plot(np.nanmean(fraction1p3,0),np.nanmean(np.squeeze(um_data['height'][p3,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 3)

    # plt.xlabel('C$_{V}$')
    plt.xlabel('TWC Cloud mask')
    plt.ylabel('Z [km]')
    # plt.title('Melt')
    plt.ylim([0,9000])
    axmajor = np.arange(0,9.01e3,1.0e3)
    axminor = np.arange(0,9.01e3,0.5e3)
    plt.yticks(axmajor)
    ax2.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,1.0])

    #----------------------
    plt.subplot(232)
    ax2 = plt.gca()
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p3,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p3,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0),np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p3,:])*1e3,0),
        np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p3,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p3,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p3,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adb', zorder = 3)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p3,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-AeroProf', zorder = 1)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p3,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)

    plt.xlabel('LWC [g m$^{-3}$]')
    # plt.title('Melt')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax2.set_yticklabels([])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,0.3])
    plt.legend()

    #-------------------------
    plt.subplot(233)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p3,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p3,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p3,:]),0),
        '--', color = 'k', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p3,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p3,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p3,:]),0),np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:])*1e3,0),
        np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p3,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p3,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p3,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p3,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    plt.xlabel('IWC [g m$^{-3}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax2.set_yticklabels([])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,0.06])

    # #-------------------------
    # plt.subplot(234)
    # ax2 = plt.gca()
    # plt.plot(np.nanmean(fraction0p4,0),np.nanmean(np.squeeze(obs_data['height'][p4,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),np.nanmean(fraction0p4,0) - np.nanstd(fraction0p4,0),
    #     np.nanmean(fraction0p4,0) + np.nanstd(fraction0p4,0), color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(fraction0p4,0) - np.nanstd(fraction0p4,0), np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction0p4,0) + np.nanstd(fraction0p4,0), np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(fraction1p4,0),np.nanmean(np.squeeze(um_data['height'][p4,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p4,:]),0),np.nanmean(fraction1p4,0) - np.nanstd(fraction1p4,0),
    #     np.nanmean(fraction1p4,0) + np.nanstd(fraction1p4,0), color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(fraction1p4,0) - np.nanstd(fraction1p4,0), np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction1p4,0) + np.nanstd(fraction1p4,0), np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(fraction4p4,0),np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0),np.nanmean(fraction4p4,0) - np.nanstd(fraction4p4,0),
    #     np.nanmean(fraction4p4,0) + np.nanstd(fraction4p4,0), color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(fraction4p4,0) - np.nanstd(fraction4p4,0), np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction4p4,0) + np.nanstd(fraction4p4,0), np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(fraction2p4,0),np.nanmean(np.squeeze(misc_data['height'][p4,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p4,:]),0),np.nanmean(fraction2p4,0) - np.nanstd(fraction2p4,0),
    #     np.nanmean(fraction2p4,0) + np.nanstd(fraction2p4,0), color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(fraction2p4,0) - np.nanstd(fraction2p4,0), np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction2p4,0) + np.nanstd(fraction2p4,0), np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(fraction3p4,0), np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0),np.nanmean(fraction3p4,0) - np.nanstd(fraction3p4,0),
    #     np.nanmean(fraction3p4,0) + np.nanstd(fraction3p4,0), color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(fraction3p4,0) - np.nanstd(fraction3p4,0), np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction3p4,0) + np.nanstd(fraction3p4,0), np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # # plt.xlabel('C$_{V}$')
    # plt.xlabel('TWC Cloud mask')
    # # plt.title('Freeze up')
    # plt.ylabel('Z [km]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,1.0])
    #
    # #--------------
    # plt.subplot(235)
    # ax2 = plt.gca()
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p4,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p4,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p4,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p4,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0),np.nanmean(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p4,:])*1e3,0),
    #     np.nanmean(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p4,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p4,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p4,:])*1e3,0),
    #     np.nanmean(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p4,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.xlabel('LWC [g m$^{-3}$]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,0.3])
    #
    # #---------------------------
    # plt.subplot(236)
    # ax2 = plt.gca()
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p4,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p4,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p4,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p4,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p4,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p4,:]),0),np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:])*1e3,0),
    #     np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p4,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p4,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p4,:])*1e3,0),
    #     np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p4,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3,
    #     np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p4,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p3,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.xlabel('IWC [g m$^{-3}$]')
    # plt.ylim([0,9e3])
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,0.06])
    #
    # print ('******')
    # print ('')
    # print ('Finished plotting! :)')
    # print ('')
    #
    # fileout = 'FIGS/Obs-' + obs_switch + 'grid_IFS_RA2M_CASIM-100_RA2T_TWCMask-LWC-IWC_p3-p4_MTThresholding-wLWCadiabatic-noOffsetLWP_226-257DOY_fixedRA2T_wSTDEV_newColours_wSetFlags.svg'
    # # plt.savefig(fileout)
    # plt.show()


    ##################################################
    ##################################################
    #### 	P5 AND P6
    ##################################################
    ##################################################
    #
    # SMALL_SIZE = 12
    # MED_SIZE = 15
    # LARGE_SIZE = 16
    #
    # plt.rc('font',size=LARGE_SIZE)
    # plt.rc('axes',titlesize=LARGE_SIZE)
    # plt.rc('axes',labelsize=LARGE_SIZE)
    # plt.rc('xtick',labelsize=LARGE_SIZE)
    # plt.rc('ytick',labelsize=LARGE_SIZE)
    # plt.rc('legend',fontsize=MED_SIZE)
    # plt.figure(figsize=(10,10))
    # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.97, left = 0.08,
    #         hspace = 0.22, wspace = 0.19)

    fraction0p5 = np.squeeze(mask0[p5,:]) # np.squeeze(obs_data['Cv'][p5,:]) #
    fraction1p5 = np.squeeze(mask1[p5,:]) # np.squeeze(um_data['model_Cv_filtered'][p5,:]) #
    fraction2p5 = np.squeeze(mask2[p5,:]) # np.squeeze(misc_data['model_Cv_filtered'][p5,:]) #
    fraction3p5 = np.squeeze(mask3[p5,:]) # np.squeeze(ifs_data['model_snow_Cv_filtered'][p5,:]) #
    fraction4p5 = np.squeeze(mask4[p5,:]) # np.squeeze(ra2t_data['model_Cv_filtered'][p5,:]) #

    fraction0p6 = np.squeeze(mask0[p6,:]) # np.squeeze(obs_data['Cv'][p6,:]) #
    fraction1p6 = np.squeeze(mask1[p6,:]) # np.squeeze(um_data['model_Cv_filtered'][p6,:]) #
    fraction2p6 = np.squeeze(mask2[p6,:]) # np.squeeze(misc_data['model_Cv_filtered'][p6,:]) #
    fraction3p6 = np.squeeze(mask3[p6,:]) # np.squeeze(ifs_data['model_snow_Cv_filtered'][p6,:]) #
    fraction4p6 = np.squeeze(mask4[p6,:]) # np.squeeze(ra2t_data['model_Cv_filtered'][p6,:]) #

    # plt.subplot(231)
    # ax2 = plt.gca()
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),np.nanmean(fraction0p5,0) - np.nanstd(fraction0p5,0),
    #     np.nanmean(fraction0p5,0) + np.nanstd(fraction0p5,0), color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(fraction0p5,0) - np.nanstd(fraction0p5,0), np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction0p5,0) + np.nanstd(fraction0p5,0), np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p5,:]),0),np.nanmean(fraction1p5,0) - np.nanstd(fraction1p5,0),
    #     np.nanmean(fraction1p5,0) + np.nanstd(fraction1p5,0), color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(fraction1p5,0) - np.nanstd(fraction1p5,0), np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction1p5,0) + np.nanstd(fraction1p5,0), np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0),np.nanmean(fraction4p5,0) - np.nanstd(fraction4p5,0),
    #     np.nanmean(fraction4p5,0) + np.nanstd(fraction4p5,0), color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(fraction4p5,0) - np.nanstd(fraction4p5,0), np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction4p5,0) + np.nanstd(fraction4p5,0), np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p5,:]),0),np.nanmean(fraction2p5,0) - np.nanstd(fraction2p5,0),
    #     np.nanmean(fraction2p5,0) + np.nanstd(fraction2p5,0), color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(fraction2p5,0) - np.nanstd(fraction2p5,0), np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction2p5,0) + np.nanstd(fraction2p5,0), np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),np.nanmean(fraction3p5,0) - np.nanstd(fraction3p5,0),
    #     np.nanmean(fraction3p5,0) + np.nanstd(fraction3p5,0), color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(fraction3p5,0) - np.nanstd(fraction3p5,0), np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(fraction3p5,0) + np.nanstd(fraction3p5,0), np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(fraction0p5,0),np.nanmean(np.squeeze(obs_data['height'][p5,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 5)
    # plt.plot(np.nanmean(fraction3p5,0), np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 4)
    # plt.plot(np.nanmean(fraction2p5,0),np.nanmean(np.squeeze(misc_data['height'][p5,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # plt.plot(np.nanmean(fraction4p5,0),np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 2)
    # plt.plot(np.nanmean(fraction1p5,0),np.nanmean(np.squeeze(um_data['height'][p5,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 3)
    #
    # # plt.xlabel('C$_{V}$')
    # plt.xlabel('TWC Cloud mask')
    # plt.ylabel('Z [km]')
    # # plt.title('Melt')
    # plt.ylim([0,9000])
    # axmajor = np.arange(0,9.01e3,1.0e3)
    # axminor = np.arange(0,9.01e3,0.5e3)
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,1.0])
    #
    # #----------------------
    # plt.subplot(232)
    # ax2 = plt.gca()
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p5,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p5,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0),np.nanmean(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p5,:])*1e3,0),
    #     np.nanmean(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p5,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p5,:])*1e3,0),
    #     np.nanmean(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(obs_data['lwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p5,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid-adb', zorder = 3)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p5,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p5,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    #
    # plt.xlabel('LWC [g m$^{-3}$]')
    # # plt.title('Melt')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,0.3])
    # plt.legend()
    #
    # #-------------------------
    # plt.subplot(233)
    # ax2 = plt.gca()
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p5,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p5,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p5,:]),0),
    #     '--', color = 'k', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p5,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p5,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3, color = 'blue', alpha = 0.05)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(um_data['height'],0),
    #     '--', color = 'darkblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p5,:]),0),np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:])*1e3,0),
    #     np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
    #     '--', color = 'steelblue', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p5,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p5,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p5,:])*1e3,0),
    #     np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(misc_data['height'],0),
    #     '--', color = 'mediumseagreen', linewidth = 0.5)
    #
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    # ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3,
    #     np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    # plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p5,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p5,:]),0),
    #     '--', color = 'gold', linewidth = 0.5)
    #
    # plt.xlabel('IWC [g m$^{-3}$]')
    # plt.ylim([0,9000])
    # plt.yticks(axmajor)
    # ax2.set_yticklabels([])
    # ax2.set_yticks(axminor, minor = True)
    # ax2.grid(which = 'major', alpha = 0.5)
    # plt.xlim([0,0.06])

    #-------------------------
    plt.subplot(234)
    ax2 = plt.gca()
    plt.plot(np.nanmean(fraction0p6,0),np.nanmean(np.squeeze(obs_data['height'][p6,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),np.nanmean(fraction0p6,0) - np.nanstd(fraction0p6,0),
        np.nanmean(fraction0p6,0) + np.nanstd(fraction0p6,0), color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(fraction0p6,0) - np.nanstd(fraction0p6,0), np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(fraction0p6,0) + np.nanstd(fraction0p6,0), np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)

    plt.plot(np.nanmean(fraction1p6,0),np.nanmean(np.squeeze(um_data['height'][p6,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p6,:]),0),np.nanmean(fraction1p6,0) - np.nanstd(fraction1p6,0),
        np.nanmean(fraction1p6,0) + np.nanstd(fraction1p6,0), color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(fraction1p6,0) - np.nanstd(fraction1p6,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(fraction1p6,0) + np.nanstd(fraction1p6,0), np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(fraction4p6,0),np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0),np.nanmean(fraction4p6,0) - np.nanstd(fraction4p6,0),
        np.nanmean(fraction4p6,0) + np.nanstd(fraction4p6,0), color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(fraction4p6,0) - np.nanstd(fraction4p6,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(fraction4p6,0) + np.nanstd(fraction4p6,0), np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(fraction2p6,0),np.nanmean(np.squeeze(misc_data['height'][p6,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p6,:]),0),np.nanmean(fraction2p6,0) - np.nanstd(fraction2p6,0),
        np.nanmean(fraction2p6,0) + np.nanstd(fraction2p6,0), color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(fraction2p6,0) - np.nanstd(fraction2p6,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(fraction2p6,0) + np.nanstd(fraction2p6,0), np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmean(fraction3p6,0), np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),np.nanmean(fraction3p6,0) - np.nanstd(fraction3p6,0),
        np.nanmean(fraction3p6,0) + np.nanstd(fraction3p6,0), color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(fraction3p6,0) - np.nanstd(fraction3p6,0), np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(fraction3p6,0) + np.nanstd(fraction3p6,0), np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    # plt.xlabel('C$_{V}$')
    plt.xlabel('TWC Cloud mask')
    # plt.title('Freeze up')
    plt.ylabel('Z [km]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax2.set_yticklabels([0,1,2,3,4,5,6,7,8,9])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,1.0])

    #--------------
    plt.subplot(235)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p6,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['lwc_adiabatic'][p6,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p6,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p6,:]),0),np.nanmean(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0),np.nanmean(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p6,:])*1e3,0),
        np.nanmean(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p6,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p6,:]),0),np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p6,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_lwc'][p6,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    plt.xlabel('LWC [g m$^{-3}$]')
    plt.ylim([0,9000])
    plt.yticks(axmajor)
    ax2.set_yticklabels([])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,0.3])

    #---------------------------
    plt.subplot(236)
    ax2 = plt.gca()
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3,np.nanmean(np.squeeze(obs_data['height'][p6,:]),0), 'k', linewidth = 3, label = 'Obs_' + obs_switch + 'grid', zorder = 3)
    ax2.fill_betweenx(np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),np.nanmean(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3, color = 'lightgrey', alpha = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(obs_data['iwc'][p6,:]),0)*1e3, np.nanmean(np.squeeze(obs_data['height'][p6,:]),0),
        '--', color = 'k', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3,np.nanmean(np.squeeze(um_data['height'][p6,:]),0), color = 'darkblue', linewidth = 2, label = 'UM_RA2M', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(um_data['height'][p6,:]),0),np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3, color = 'blue', alpha = 0.05)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(um_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(um_data['height'],0),
        '--', color = 'darkblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3,np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0), color = 'steelblue', linewidth = 2, label = 'UM_RA2T', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ra2t_data['height'][p6,:]),0),np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:])*1e3,0),
        np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3, color = 'lightblue', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ra2t_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(ra2t_data['height'],0),
        '--', color = 'steelblue', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3,np.nanmean(np.squeeze(misc_data['height'][p6,:]),0), color = 'mediumseagreen', linewidth = 2, label = 'UM_CASIM-100', zorder = 1)
    ax2.fill_betweenx(np.nanmean(np.squeeze(misc_data['height'][p6,:]),0),np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p6,:])*1e3,0),
        np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3, color = 'mediumaquamarine', alpha = 0.15)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(misc_data['model_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(misc_data['height'],0),
        '--', color = 'mediumseagreen', linewidth = 0.5)

    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3,np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0), color = 'gold', linewidth = 2, label = 'ECMWF_IFS', zorder = 2)
    ax2.fill_betweenx(np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3,
        np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3, color = 'navajowhite', alpha = 0.35)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3 - np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)
    plt.plot(np.nanmean(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3 + np.nanstd(np.squeeze(ifs_data['model_snow_iwc_filtered'][p6,:]),0)*1e3, np.nanmean(np.squeeze(ifs_data['height'][p6,:]),0),
        '--', color = 'gold', linewidth = 0.5)

    plt.xlabel('IWC [g m$^{-3}$]')
    plt.ylim([0,9e3])
    plt.yticks(axmajor)
    ax2.set_yticklabels([])
    ax2.set_yticks(axminor, minor = True)
    ax2.grid(which = 'major', alpha = 0.5)
    plt.xlim([0,0.06])

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    fileout = 'FIGS/Obs-' + obs_switch + 'grid_IFS_RA2M_CASIM-AeroProf_RA2T_TWCMask-LWC-IWC_p3-p6_MTThresholding-wLWCadiabatic-noOffsetLWP_226-257DOY_fixedRA2T_wSTDEV_newColours_wSetFlags.svg'
    # plt.savefig(fileout)
    plt.show()

    Zindex1 = np.where(np.round(np.nanmean(um_data['height'],0),-2) <= 4000)#<= 2e3) # == 5.e+02)#
    Zindex3 = np.where(np.round(np.nanmean(ifs_data['height'],0),-2) <= 4000)# <= 2e3) #== 5.e+02)#
    print ('Zindex1 = ')
    print (np.nanmean(um_data['height'][:,Zindex1[0][13]],0))
    print ('Zindex3 = ')
    print (np.nanmean(ifs_data['height'][:,Zindex3[0][17]],0))
    # print ('UM_RA2M = ')
    # print (np.nanmean(np.squeeze(um_data['model_lwc'][p3,:]),0)*1e3)
    # print ('UM_RA2T = ')
    # print (np.nanmean(np.squeeze(ra2t_data['model_lwc'][p3,:]),0)*1e3)
    # print ('UM_CASIM-100 = ')
    # print (np.nanmean(np.squeeze(misc_data['model_lwc'][p6,:]),0)*1e3)
    print ('ECMWF_IFS = ')
    print (np.round(np.nanmean(np.squeeze(ifs_data['model_lwc'][p6,13]),0)*1e3,4))
    print (np.nanmean(np.squeeze(fraction3p6[:,17]),0))
    print ('Obs = ')
    print (np.round(np.nanmean(np.squeeze(obs_data['lwc_adiabatic'][p6,13]),0)*1e3,4))
    print (np.nanmean(np.squeeze(fraction0p6[:,17]),0))

def plot_BiasCorrelation(obs_data, um_data, misc_data, ifs_data, ra2t_data, doy, obs_switch):

    ###################################
    ## Directly compare LWC and T biases
    ###################################

    print ('******')
    print ('')
    print ('Comparing radiosonde (T) and Cloudnet (LWC) biases:')
    print ('')


    ######################################################################################################
    ######################################################################################################
    #####           load sonde biases first
    ######################################################################################################
    ######################################################################################################

    data1 = np.load('UM_RA2M_SondeBiases.npy', allow_pickle=True, encoding='latin1').item()
    data2 = np.load('UM_CASIM-100_SondeBiases.npy', allow_pickle=True, encoding='latin1').item()
    data3 = np.load('ECMWF_IFS_SondeBiases.npy', allow_pickle=True, encoding='latin1').item()
    data4 = np.load('UM_RA2T_SondeBiases.npy', allow_pickle=True, encoding='latin1').item()
    data5 = np.load('UM_GLM_SondeBiases.npy', allow_pickle=True, encoding='latin1').item()
    obs = np.load('SondeData_reGridded.npy', allow_pickle=True, encoding='latin1').item()

    plt.figure()
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'])
    # plt.plot(np.nanmedian(obs['temp_driftSondes_UM'] + 273.15,0),data1['universal_height'])
    plt.show()

    ######################################################################################################
    ######################################################################################################
    #####           next, calculate LWC biases from Cloudnet
    ######################################################################################################
    ######################################################################################################

    # print (data1['temp_anomalies'].shape)
    # print (um_data['model_lwc'][::6,data1['universal_height_UMindex']].shape)


    ###----------------------------------------------------------------
    ###         2) TWC Method
    ###----------------------------------------------------------------
    #### set flagged um_data to nans
    # obs_data['iwc'][obs_data['iwc'] < 0] = 0.0
    # um_data['model_iwc_filtered'][um_data['model_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] < 0.0] = 0.0
    # ifs_data['model_snow_iwc_filtered'][ifs_data['model_snow_iwc_filtered'] >= 5.0e-3] = np.nan
    # misc_data['model_iwc_filtered'][misc_data['model_iwc_filtered'] < 0.0] = 0.0
    # ra2t_data['model_iwc_filtered'][ra2t_data['model_iwc_filtered'] < 0.0] = 0.0

    #### set flagged um_data to nans
    # obs_data['lwc'][obs_data['lwc'] == -999] = 0.0
    # obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] == -999] = 0.0
    # obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] == -999] = 0.0
    # um_data['model_lwc'][um_data['model_lwc'] < 0.0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] < 0.0] = 0.0
    # ifs_data['model_lwc'][ifs_data['model_lwc'] >= 0.4] = np.nan
    # misc_data['model_lwc'][misc_data['model_lwc'] < 0.0] = 0.0
    # ra2t_data['model_lwc'][ra2t_data['model_lwc'] < 0.0] = 0.0

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    um_data['model_twc'] = um_data['model_lwc'] + um_data['model_iwc_filtered']
    misc_data['model_twc'] = misc_data['model_lwc'] + misc_data['model_iwc_filtered']
    ifs_data['model_twc'] = ifs_data['model_lwc'] + ifs_data['model_snow_iwc_filtered']
    ra2t_data['model_twc'] = ra2t_data['model_lwc'] + ra2t_data['model_iwc_filtered']

    ####-------------------------------------------------------------------------
    ### from Michael's paper:
    ####    "we consider a grid point cloudy when cloud water exceeded
    ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    ####    interpolation in between."
    ####-------------------------------------------------------------------------
    twc_thresh_um = np.zeros([np.size(um_data['model_twc'],1)])
    twc_thresh_ifs = np.zeros([np.size(ifs_data['model_twc'],1)])

    ####-------------------------------------------------------------------------
    ### first look at values below 1 km
    ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    um_lt1km = np.where(um_data['height'][0,:]<=1e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]<=1e3)
    twc_thresh_um[um_lt1km] = 1e-6
    twc_thresh_ifs[ifs_lt1km] = 1e-6

    ####-------------------------------------------------------------------------
    ### next look at values above 4 km
    ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    um_lt1km = np.where(um_data['height'][0,:]>=4e3)
    ifs_lt1km = np.where(ifs_data['height'][0,:]>=4e3)
    twc_thresh_um[um_lt1km] = 1e-7
    twc_thresh_ifs[ifs_lt1km] = 1e-7

    ### find heights not yet assigned
    um_intZs = np.where(twc_thresh_um == 0.0)
    ifs_intZs = np.where(twc_thresh_ifs == 0.0)

    ### interpolate for twc_thresh_um
    x = [1e-6, 1e-7]
    y = [1e3, 4e3]
    f = interp1d(y, x)
    twc_thresh_um[um_intZs] = f(um_data['height'][0,um_intZs].data)

    ### interpolate for twc_thresh_ifs
    twc_thresh_ifs[ifs_intZs] = f(ifs_data['height'][0,ifs_intZs].data)

    ### plot profile of threshold as sanity check
    # plt.plot(twc_thresh_um, um_data['height'][0,:])
    # plt.plot(twc_thresh_ifs, ifs_data['height'][0,:]); plt.show()
    # print (twc_thresh_um)

    for t in range(0,np.size(um_data['model_twc'],0)):
        for k in range(0,np.size(um_data['model_twc'],1)):
            if obs_switch == 'UM':
                if obs_data['twc'][t,k] < twc_thresh_um[k]:
                    obs_data['twc'][t,k] = np.nan
                    obs_data['lwc'][t,k] = np.nan
                if obs_data['twc_ad'][t,k] < twc_thresh_um[k]:
                    obs_data['twc_ad'][t,k] = np.nan
                    obs_data['lwc_adiabatic'][t,k] = np.nan
                if obs_data['twc_ad_nolwp'][t,k] < twc_thresh_um[k]:
                    obs_data['twc_ad_nolwp'][t,k] = np.nan
                    obs_data['lwc_adiabatic_inc_nolwp'][t,k] = np.nan
            if um_data['model_twc'][t,k] < twc_thresh_um[k]:
                um_data['model_twc'][t,k] = np.nan
                um_data['model_lwc'][t,k] = np.nan
            if misc_data['model_twc'][t,k] < twc_thresh_um[k]:
                misc_data['model_twc'][t,k] = np.nan
                misc_data['model_lwc'][t,k] = np.nan
            if ra2t_data['model_twc'][t,k] < twc_thresh_um[k]:
                ra2t_data['model_twc'][t,k] = np.nan
                ra2t_data['model_lwc'][t,k] = np.nan
        for k in range(0,np.size(ifs_data['model_twc'],1)):
            if ifs_data['model_twc'][t,k] < twc_thresh_ifs[k]:
                ifs_data['model_twc'][t,k] = np.nan
                ifs_data['model_lwc'][t,k] = np.nan

    um_data['lwc_bias'] = (um_data['model_lwc'][::6,data1['universal_height_UMindex']] - obs_data['lwc_adiabatic'][::6,data1['universal_height_UMindex']])*1e3
    misc_data['lwc_bias'] = (misc_data['model_lwc'][::6,data1['universal_height_UMindex']] - obs_data['lwc_adiabatic'][::6,data1['universal_height_UMindex']])*1e3
    ifs_data['lwc_bias'] = (ifs_data['model_lwc'][::6,data1['universal_height_UMindex']] - obs_data['lwc_adiabatic'][::6,data1['universal_height_UMindex']])*1e3

    um_data['lwc_bias'][um_data['lwc_bias'] == 0] = np.nan
    misc_data['lwc_bias'][misc_data['lwc_bias'] == 0] = np.nan
    ifs_data['lwc_bias'][ifs_data['lwc_bias'] == 0] = np.nan

    ###############################
    ### TEST FIGS
    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16
    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=LARGE_SIZE)
    plt.rc('axes',labelsize=LARGE_SIZE)
    plt.rc('xtick',labelsize=LARGE_SIZE)
    plt.rc('ytick',labelsize=LARGE_SIZE)
    plt.rc('legend',fontsize=LARGE_SIZE)

    plt.figure(figsize=(6,6))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.9, left = 0.1,
            hspace = 0.4, wspace = 0.2)
    plt.subplot(121)
    plt.plot(np.nanmedian(data1['temp_anomalies'],1),data1['universal_height'],'.-' ,color = 'darkblue', zorder = 3)
    plt.ylim([0, 1e4])
    plt.subplot(122)
    plt.plot(np.nanmean(um_data['lwc_bias'],0),um_data['height'][0,data1['universal_height_UMindex']])
    plt.ylim([0, 1e4])
    plt.xlim([0, 0.1])
    plt.show()

    Z1index = np.where(data1['universal_height'] <= 3e3)
    Z2index = np.where(um_data['height'][0,:] <= 3.2e3)

    print (Z1index[0].shape)
    print (Z2index[0].shape)

    print (data1['universal_height'][Z1index])
    print (um_data['height'][0,2:Z2index[0][-1]])

    plt.figure()
    plt.subplots_adjust(top = 0.9, bottom = 0.15, right = 0.9, left = 0.2,
            hspace = 0.4, wspace = 0.2)
    plt.scatter(np.ndarray.flatten(data1['temp_anomalies'][Z1index,:]),np.ndarray.flatten(um_data['lwc_bias'][:,2:Z2index[0][-1]]), s = 3, color = 'darkblue')
    plt.scatter(np.ndarray.flatten(data2['temp_anomalies'][Z1index,:]),np.ndarray.flatten(misc_data['lwc_bias'][:,2:Z2index[0][-1]]), s = 3, color = 'mediumseagreen')
    plt.scatter(np.ndarray.flatten(data3['temp_anomalies'][Z1index,:]),np.ndarray.flatten(ifs_data['lwc_bias'][:,2:Z2index[0][-1]]), s = 3, color = 'gold')
    plt.ylabel('LWC bias [g m$^{-3}$]')
    plt.xlabel('T bias [K]')
    plt.show()

def interpCloudnet(obs_data):
    #interpolates missing times up to 1 hour
    from scipy.interpolate import interp1d
    from manipFuncts import nanhelper
    print ('*******')
    print ('Interpolate obs cloudnet field for continuous array:')

    varlist = ['Cv', 'lwc', 'iwc']
    for var in varlist:
        ### remove bad and flagged data
        obs_data[var][obs_data[var] < 0.0] = np.nan

        ### save relevant fields as tempvars for ease
        cv = np.copy(obs_data[var].data)
        times = np.copy(obs_data['time'].data)
        height = np.copy(obs_data['height'][0,:])        ### height array constant in time, so just take first column
        nans,id=nanhelper(cv)
    #    from IPython import embed; embed()
        #for i in range(0,len(height)):
        for i in range(0,len(height)):
            tmp=id(~nans[:,i])
            idtmp=np.squeeze(np.nonzero(np.diff(np.append([0],tmp))>4)) #only interp ,gaps of max 3 timesteps
            nanint=(nans[:,i])
            if idtmp.tolist():
                ide=int2list(tmp[idtmp])
                ids=int2list(tmp[idtmp-1]+1)
                for m in range(0,len(ids)):
                    nanint[ids[m]:ide[m]]=np.False_
            if any(nanint==True) and any(nanint==False):
                cv[nanint,i]=np.interp(id(nanint),id(~nanint),cv[~nanint,i])
        ### check if interpolation worked
        # fig = plt.figure(figsize=(9,6))
        # ax=plt.gca()
        # plt.subplot(211)
        # plt.pcolor(obs_data['time'],np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['Cv']))
        # plt.clim(vmin=0,vmax=1)
        # plt.subplot(212)
        # plt.pcolor(obs_data['time'],np.squeeze(obs_data['height'][0,:]), np.transpose(cv))
        # plt.clim(vmin=0,vmax=1)
        #
        #
        ### save back to dictionary after completion of updates
        obs_data[var] = cv
    print('done')
    print ('*******')
    print ('')

    return obs_data

def buildNaNMask(obs_data):

    print ('*******')
    print ('Only include model data where we have observations:')
    print ('*******')
    print ('')

    ### build nanmask
    nanmask = np.zeros([np.size(obs_data['time']), np.size(obs_data['Cv'],1)])
    nanindex = np.zeros([np.size(obs_data['time'])])
    wcindex = np.zeros([np.size(obs_data['time'])])
    wc0index = np.zeros([np.size(obs_data['time'])])
    lwpindex = np.zeros([np.size(obs_data['time'])])
    # print(nanindex.shape)
    for i in range(len(obs_data['time'])):
        if np.isnan(np.nanmean(obs_data['Cv'][i,:], 0)):       ## if there are only nans in the time profile
            nanmask[i,:] = 1.0
            #nanmask[i-1,:] = 1.0   ##why setting -1 and +1 to nan???
            #nanmask[i+1,:] = 1.0
            nanindex[i] = 1
        if np.logical_and(np.isnan(np.nanmean(obs_data['lwc_adiabatic'][i,:], 0)), np.isnan(np.nanmean(obs_data['iwc'][i,:], 0))):       ## if both wc profiles contain only nans
            wcindex[i] = 1
        elif np.logical_and(np.nanmean(obs_data['lwc'][i,:], 0) == 0, np.nanmean(obs_data['iwc'][i,:], 0) == 0):       ## if both wc profiles contain only zeros
            wc0index[i] = 1
        elif np.isnan(obs_data['lwp'][i,0]):       ## if there are only nans in the lwp timeseries
            lwpindex[i] = 1
    nanmask[nanmask == 0.0] = np.nan
    nanind = np.where(nanindex == 1)
    wcind = np.where(wcindex == 1)
    wc0ind = np.where(wc0index == 1)
    lwpind = np.where(lwpindex == 1)

    return nanind, nanmask, wcind, wc0ind, lwpind

def buildNaNMaskadv(obs_data):

    print ('*******')
    print ('Only include model data where we have observations:')
    print ('*******')
    print ('')

    ### build nanmask
    nanmask = np.zeros([np.size(obs_data['time']), np.size(obs_data['Cv_adv'],1)])
    nanindex = np.zeros([np.size(obs_data['time'])])
    wcindex = np.zeros([np.size(obs_data['time'])])
    wc0index = np.zeros([np.size(obs_data['time'])])
    lwpindex = np.zeros([np.size(obs_data['time'])])
    # print(nanindex.shape)
    for i in range(len(obs_data['time'])):
        if np.isnan(np.nanmean(obs_data['Cv_adv'][i,:], 0)):       ## if there are only nans in the time profile
            nanmask[i,:] = 1.0
            #nanmask[i-1,:] = 1.0   ##why setting -1 and +1 to nan???
            #nanmask[i+1,:] = 1.0
            nanindex[i] = 1
        if np.logical_and(np.isnan(np.nanmean(obs_data['lwc_adv_adiabatic'][i,:], 0)), np.isnan(np.nanmean(obs_data['iwc_adv'][i,:], 0))):       ## if both wc profiles contain only nans
            wcindex[i] = 1
        elif np.logical_and(np.nanmean(obs_data['lwc_adv'][i,:], 0) == 0, np.nanmean(obs_data['iwc_adv'][i,:], 0) == 0):       ## if both wc profiles contain only zeros
            wc0index[i] = 1
        elif np.isnan(obs_data['lwp_adv'][i,0]):       ## if there are only nans in the lwp timeseries
            lwpindex[i] = 1
    nanmask[nanmask == 0.0] = np.nan
    nanind = np.where(nanindex == 1)
    wcind = np.where(wcindex == 1)
    wc0ind = np.where(wc0index == 1)
    lwpind = np.where(lwpindex == 1)

    return nanind, nanmask, wcind, wc0ind, lwpind

def setFlags(obs_data, um_data, misc_data, ra2t_data, obs_var_list, um_var_list, misc_var_list,  ra2t_var_list):

    print ('*******')
    print ('Set flagged data to nan:')
    print ('*******')
    print ('')

    ###----------------------------------------------------------------
    ###         Set flagged data to nans
    ###----------------------------------------------------------------
    for c in range(0,3):
        for j in range(0,len(obs_var_list[c])):
            obs_data[obs_var_list[c][j]][obs_data[obs_var_list[c][j]] == -999] = np.nan
            obs_data[obs_var_list[c][j]][obs_data[obs_var_list[c][j]] < 0] = 0.0
        for j in range(0,len(um_var_list[c])):
            um_data[um_var_list[c][j]][um_data[um_var_list[c][j]] == -999] = np.nan
            um_data[um_var_list[c][j]][um_data[um_var_list[c][j]] < 0] = 0.0
        for j in range(0,len(misc_var_list[c])):
            misc_data[misc_var_list[c][j]][misc_data[misc_var_list[c][j]] == -999] = np.nan
            misc_data[misc_var_list[c][j]][misc_data[misc_var_list[c][j]] < 0] = 0.0
        for j in range(0,len(ra2t_var_list[c])):
            ra2t_data[ra2t_var_list[c][j]][ra2t_data[ra2t_var_list[c][j]] == -999] = np.nan
            ra2t_data[ra2t_var_list[c][j]][ra2t_data[ra2t_var_list[c][j]] < 0] = 0.0

    # for c in range(0,2):
    #     for j in range(0,len(monc_var_list[c])):
    #         monc_data[monc_var_list[c][j]][monc_data[monc_var_list[c][j]]== -999] = np.nan

    return obs_data, um_data, misc_data,  ra2t_data

################################################################################
################################################################################
def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    names = ['20180912_oden_','20180913_oden_']
    # names = ['20180814_oden_','20180815_oden_','20180816_oden_',
    #         '20180817_oden_','20180818_oden_','20180819_oden_','20180820_oden_',
    #         '20180821_oden_','20180822_oden_','20180823_oden_','20180824_oden_',
    #         '20180825_oden_','20180826_oden_','20180827_oden_','20180828_oden_',
    #         '20180829_oden_','20180830_oden_','20180831_oden_','20180901_oden_',
    #         '20180902_oden_','20180903_oden_','20180904_oden_','20180905_oden_',
    #         '20180906_oden_','20180907_oden_','20180908_oden_','20180909_oden_',
    #         '20180911_oden_','20180912_oden_','20180913_oden_',
    #         '20180914_oden_']

    sdate = dtime.datetime.strptime('2018091301','%Y%m%d%H')
    edate = dtime.datetime.strptime('2018091314','%Y%m%d%H')
    # sdate = dtime.datetime.strptime('2018081400','%Y%m%d%H')
    # edate = dtime.datetime.strptime('2018091500','%Y%m%d%H')
    dates = [date2datenum(sdate),date2datenum(edate)]
    moccha_missing_files = ['20180813_oden_','20180910_oden_']   ### cloud radar not working    #,'20180914_oden_'

    #---- MONC SPIN UP TIME
    monc_spin = 6 *60 *60

    ### Set output directory for plots
    plots_out_dir='/nfs/a96/MOCCHA/working/jutta/plots/CaseStudies/ModelComparison/nointerp/'
    if not os.path.exists(plots_out_dir):
        os.makedirs(plots_out_dir)

    ### Choose observations vertical gridding used in Cloudnet processing (UM/IFS/RADAR)
    obs_switch = 'UM'

    um_root_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/data/calibrated/metum/'
    cn_obs_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/Halo/measurements_V7/'
    cn_um_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/metum/V7/'
    cn_monc_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/1_control/'
    monc_filename = cn_monc_dir + 'moccha_casim_dg_50400.nc'
    ### -----------------------------------------------------------------
    ### CHOSEN RUN - MODEL DATA
    out_dir1 = '25_u-cc568_RA2M_CON/'
    out_dir3 = '24_u-cc324_RA2T_CON/'
    out_dir2 = '23_u-cc278_RA1M_CASIM/'

    ### -----------------------------------------------------------------
    ### CHOSEN RUN - CLOUDNET DATA
    cn_um_out_dir = ['cloud-fraction-metum-grid/2018/',
                     'lwc-scaled-metum-grid/2018/',
                     'iwc-Z-T-metum-grid/2018/']
    cn_obs_out_dir = ['cloud-fraction-metum-grid/2018/qual70/',
                        'lwc-adiabatic-metum-grid/2018/',
                        'iwc-Z-T-metum-grid/2018/']


    ######## lwc-adiabatic-metum-grid/2018/
    ########-> liquid water content derived using measurements averaged on to model grid
    ### cloud-fraction-metum-grid/2018/ + cloud-fraction-ecmwf-grid/2018/
    ###             -> cloud fraction both from a forecast model and derived from the high-resolution observations on the grid of that model.
    ### lwc-scaled-metum-grid/2018/ + lwc-scaled-ecmwf-grid/2018/ + lwc-scaled-adiabatic/2018/
    ###             -> dataset contains liquid water content derived using radar/lidar cloud boundaries and liquid water path from dual-wavelength
    ###                 microwave radiometers, averaged on to the grid of a forecast model.
    ###                 It also contains the liquid water content and liquid water path from that model, so may be used to calculate statistics
    ###                 quantifying the model performance.
    ### iwc-Z-T-metum-grid/2018/ + iwc-Z-T-ecmwf-grid/2018/ + iwc-Z-T-method/2018/
    ###             -> dataset contains ice water content derived using radar reflectivity and temperature, averaged on to the grid of a forecast
    ###                 model. It also contains the ice water content from that model, so may be used to calculate statistics quantifying the
    ###                 model performance.

    # # -------------------------------------------------------------
    # # Load cloudnet data
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin nc read in at ' + time.strftime("%c"))
    print (' ')

    for i in range(0,len(names)):
        dstr=names[i][0:8]
        datenum = date2datenum(dtime.datetime.strptime(dstr,'%Y%m%d'))

        print ('load cloudnet data:')
        cn_filename_um = [cn_um_dir + out_dir1 + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + out_dir1 + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + out_dir1 + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc']
        cn_filename_misc = [cn_um_dir + out_dir2 + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + out_dir2 + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + out_dir2 + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc']
        cn_filename_ra2t = [cn_um_dir + out_dir3 + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + out_dir3 + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + out_dir3 + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc']
        cn_filename_obs = [cn_obs_dir + cn_obs_out_dir[0] + names[i] + cn_obs_out_dir[0][:-13] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[1] + names[i] + cn_obs_out_dir[1][:-6] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[2] + names[i] + cn_obs_out_dir[2][:-6] + '.nc']
        #print (cn_filename_um)
        #print (cn_filename_misc)
        #print (cn_filename_ra2t)
        #print (cn_filename_obs)
        print ('')

    ### --------------------------------------------------------------------
    ###     READ IN ALL CLOUDNET FILES: reinitialise diagnostic dictionaries
    ### --------------------------------------------------------------------
        print ('Loading multiple diagnostics from cloudnet results:')
        cn_nc1 = {}
        cn_nc2 = {}
        cn_nc3 = {}
        cn_nc0 = {}
        for c in range(0,3):
            cn_nc1[c] = Dataset(cn_filename_um[c],'r')
            cn_nc2[c] = Dataset(cn_filename_misc[c],'r')
            cn_nc3[c] = Dataset(cn_filename_ra2t[c],'r')
            cn_nc0[c] = Dataset(cn_filename_obs[c],'r')

        # -------------------------------------------------------------
        print ('')

        ### --------------------------------------------------------------------
        ###     LOAD CLOUDNET DIAGS INTO DICTIONARY
        ### --------------------------------------------------------------------
        #### LOAD IN SPECIFIC DIAGNOSTICS
        obs_var_list = [['Cv', 'Cv_adv'],
                    ['lwc','lwc_adv','lwp','lwp_adv','lwc_adiabatic','lwc_adv_adiabatic','lwc_adiabatic_inc_nolwp','lwc_adv_adiabatic_inc_nolwp'],
                    ['height','iwc','iwc_adv']]

        um_var_list = [['Cv','model_Cv_filtered','model_temperature'],
                ['lwc','lwp','model_lwc','model_lwp'],
                ['height','iwc','model_iwc','model_iwc_filtered']]   ### time always read in separately
        misc_var_list=um_var_list
        ra2t_var_list=um_var_list
        ### --------------------------------------------------------------------
        ### create arrays for all cloudnet data
        ### --------------------------------------------------------------------
        if i == 0:
            ### --------------------------------------------------------------------
            ### initialise cloudnet data dictionaries
            ### --------------------------------------------------------------------
            um_data = {}
            misc_data = {}
            ra2t_data = {}
            obs_data = {}
            ### --------------------------------------------------------------------
            ### create time arrays for all cloudnet data
            ### --------------------------------------------------------------------
            if obs_switch == 'RADAR':
                time_obs = datenum + np.float64((cn_nc0[1].variables['time'][:])/24.0)
            else:
                time_obs = datenum + np.float64((cn_nc0[0].variables['time'][:])/24.0)
            time_um = datenum + np.float64((cn_nc1[0].variables['time'][:])/24.0)
            time_misc = datenum + np.float64((cn_nc2[0].variables['time'][:])/24.0)
            time_ra2t = datenum + np.float64((cn_nc3[0].variables['time'][:])/24.0)
            ### --------------------------------------------------------------------
            ### loop over each Cloudnet class
            ### --------------------------------------------------------------------
            for c in range(0,3):
                ### --------------------------------------------------------------------
                ### load in initial obs data
                ### --------------------------------------------------------------------
                if obs_switch == 'RADAR':
                    for j in range(0,len(obs_var_list[c])):
                        if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:  # 1d timeseries only
                            obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                else:
                    for j in range(0,len(obs_var_list[c])):
                        if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:  # 1d timeseries only
                            obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            obs_data[obs_var_list[c][j]] = cn_nc0[c].variables[obs_var_list[c][j]][:]
                    # lwp = []
                    # lwp = cn_nc0[1].variables['lwp'][:]
                ### --------------------------------------------------------------------
                ### load in initial UM_RA2M data
                ### --------------------------------------------------------------------
                for j in range(0,len(um_var_list[c])):
                    if np.ndim(cn_nc1[c].variables[um_var_list[c][j]]) == 1:  # 1d timeseries only
                        um_data[um_var_list[c][j]] = cn_nc1[c].variables[um_var_list[c][j]][:]
                    else:                                   # 2d column um_data
                        um_data[um_var_list[c][j]] = cn_nc1[c].variables[um_var_list[c][j]][:]
                ### --------------------------------------------------------------------
                ### load in initial UM_CASIM-100 data
                ### --------------------------------------------------------------------
                for j in range(0,len(misc_var_list[c])):
                    if np.ndim(cn_nc2[c].variables[misc_var_list[c][j]]) == 1:  # 1d timeseries only
                        misc_data[misc_var_list[c][j]] = cn_nc2[c].variables[misc_var_list[c][j]][:]
                    else:                                   # 2d column um_data
                        misc_data[misc_var_list[c][j]] = cn_nc2[c].variables[misc_var_list[c][j]][:]
                ### --------------------------------------------------------------------
                ### load in initial UM_RA2T data
                ### --------------------------------------------------------------------
                for j in range(0,len(ra2t_var_list[c])):
                    if np.ndim(cn_nc3[c].variables[ra2t_var_list[c][j]]) == 1:  # 1d timeseries only
                        ra2t_data[ra2t_var_list[c][j]] = cn_nc3[c].variables[ra2t_var_list[c][j]][:]
                    else:                                   # 2d column um_data
                        ra2t_data[ra2t_var_list[c][j]] = cn_nc3[c].variables[ra2t_var_list[c][j]][:]
        ### --------------------------------------------------------------------
        ### fill arrays with remaining data
        ### --------------------------------------------------------------------
        else:
            if obs_switch == 'RADAR':
                time_obs = np.append(time_obs,  datenum +  np.float64((cn_nc0[1].variables['time'][:])/24.0))
            else:
                time_obs = np.append(time_obs,  datenum + np.float64((cn_nc0[0].variables['time'][:])/24.0))
            time_um = np.append(time_um, datenum + np.float64((cn_nc1[0].variables['time'][:])/24.0))
            time_misc = np.append(time_misc,  datenum +  np.float64((cn_nc2[0].variables['time'][:])/24.0))
            time_ra2t = np.append(time_ra2t,  datenum +  np.float64((cn_nc3[0].variables['time'][:])/24.0))

            ### --------------------------------------------------------------------
            ### loop over all Cloudnet classes
            ### --------------------------------------------------------------------
            for c in range(0,3):
                ### --------------------------------------------------------------------
                ### append rest of obs data
                ### --------------------------------------------------------------------
                if obs_switch == 'RADAR':
                    for j in range(0,len(obs_var_list[c])):
                        if obs_var_list[c][j] == 'height':
                            continue
                        elif np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:
                            obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:])
                        elif np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 2:
                            obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:],0)
                else:
                    for j in range(0,len(obs_var_list[c])):
                        if np.ndim(cn_nc0[c].variables[obs_var_list[c][j]]) == 1:
                            obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:])
                        elif np.sum(cn_nc0[c].variables[obs_var_list[c][j]].shape) == 71:
                            continue
                        else:
                            obs_data[obs_var_list[c][j]] = np.append(obs_data[obs_var_list[c][j]],cn_nc0[c].variables[obs_var_list[c][j]][:],0)
                ### --------------------------------------------------------------------
                ### append rest of UM_RA2M data
                ### --------------------------------------------------------------------
                for j in range(0,len(um_var_list[c])):
                    # print (cn_nc1[c])
                    # print (um_var_list[c][j])
                    if np.ndim(cn_nc1[c].variables[um_var_list[c][j]]) == 1:
                        um_data[um_var_list[c][j]] = np.append(um_data[um_var_list[c][j]],cn_nc1[c].variables[um_var_list[c][j]][:])
                    else:
                        um_data[um_var_list[c][j]] = np.append(um_data[um_var_list[c][j]],cn_nc1[c].variables[um_var_list[c][j]][:],0)
                ### --------------------------------------------------------------------
                ### append rest of UM_CASIM-100 data
                ### --------------------------------------------------------------------
                for j in range(0,len(misc_var_list[c])):
                    if np.ndim(cn_nc2[c].variables[misc_var_list[c][j]]) == 1:
                        misc_data[misc_var_list[c][j]] = np.append(misc_data[misc_var_list[c][j]],cn_nc2[c].variables[misc_var_list[c][j]][:])
                    else:
                        misc_data[misc_var_list[c][j]] = np.append(misc_data[misc_var_list[c][j]],cn_nc2[c].variables[misc_var_list[c][j]][:],0)
                ### --------------------------------------------------------------------
                ### append rest of UM_RA2T data
                ### --------------------------------------------------------------------
                for j in range(0,len(ra2t_var_list[c])):
                    if np.ndim(cn_nc3[c].variables[ra2t_var_list[c][j]]) == 1:
                        ra2t_data[ra2t_var_list[c][j]] = np.append(ra2t_data[ra2t_var_list[c][j]],cn_nc3[c].variables[ra2t_var_list[c][j]][:])
                    else:
                        ra2t_data[ra2t_var_list[c][j]] = np.append(ra2t_data[ra2t_var_list[c][j]],cn_nc3[c].variables[ra2t_var_list[c][j]][:],0)
            # lwp = np.append(lwp, cn_nc0[1].variables['lwp'][:],0)

        ### --------------------------------------------------------------------
        ### Close cloudnet netCDF files
        ### --------------------------------------------------------------------
        for c in range(0,3): cn_nc0[c].close()
        for c in range(0,3): cn_nc1[c].close()
        for c in range(0,3): cn_nc2[c].close()
        for c in range(0,3): cn_nc3[c].close()

        # print (time_obs.shape)
        # print (time_ifs.shape)
        # print (obs_data['lwp'].shape)
        # print (time_obs)
        # print ((len(doy)*24.0))
        # print (len(moccha_names))
        # obs_data['lwp'] = lwp

        print ('Loaded!')
        print ('')
        print ('*************** NEXT:')
        print ('')
    #################################################################
    ## save time to dictionaries now we're not looping over all diags anymore
    #################################################################
    um_data['time'] = time_um
    misc_data['time'] = time_misc
    ra2t_data['time'] = time_ra2t
    obs_data['time'] = time_obs

    #################################################################
    ### if RADAR flag used, force 2D height array for function compatibility
    #################################################################
    # if obs_switch == 'RADAR':
    #     tmp_height = np.zeros([np.size(obs_data['time']), np.size(obs_data['height'])])
    #     for t in range(0,len(obs_data['time'])): tmp_height[t,:] = obs_data['height'][:]
    #     obs_data['height'] = tmp_height

    #################################################################
    ###     READ IN MONC DATA
    #################################################################
    print ('Loading monc results:')
    print ('')
    ###1d variables, 2d variables ( 200x84)
    monc_var_list =[['time_series_2_60','time_series_20_600' ,'z', 'LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
                    ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
                    'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
                    'graupel_mmr_mean']]

    print(monc_var_list)


    ncm = {}
    monc_data = {}
    ncm = Dataset(monc_filename,'r')
    for var in ncm.variables: print (var)
    for c in range(0,2):
        for j in range(0,len(monc_var_list[c])):
            monc_data[monc_var_list[c][j]] = ncm.variables[monc_var_list[c][j]][:]

    monc_data['time1']=monc_data['time_series_2_60'] #1d data
    monc_data['time2']=monc_data['time_series_20_600'] #2d data
    monc_data.pop('time_series_2_60')
    monc_data.pop('time_series_20_600')

    print ('Loaded!')

    monc_var_list = list(monc_data.keys())
    ## remove spin up time from monc data1
    id =np.argwhere(monc_data['time1']<=monc_spin)
    for j in range(0,len(monc_var_list)):
        if any(monc_data[monc_var_list[j].shape == len(monc_data['time1']))
            monc_data[monc_var_list[c][j]] = monc_data[monc_var_list[c][j]]



    ##################################################################################################################################
    ## -------------------------------------------------------------
    ## maximise obs data available and build mask for available data
    ## -------------------------------------------------------------
    print('setting missing data to nan and interpolate missing obs')
    obs_data, um_data, misc_data, ra2t_data = setFlags(obs_data, um_data, misc_data, ra2t_data, obs_var_list, um_var_list, misc_var_list, ra2t_var_list)
    #obs_data = interpCloudnet(obs_data)
    nanind, nanmask, wcind, wc0ind, lwpind = buildNaNMask(obs_data)
    nanindadv, nanmaskadv, wcindadv, wc0indadv, lwpindadv = buildNaNMaskadv(obs_data)

    varlist_obs = ['Cv', 'lwc_adiabatic', 'iwc', 'lwp']
    varlist_um = ['model_Cv_filtered', 'model_lwc', 'model_iwc_filtered', 'model_lwp']
    # print(um_data.keys())

    ## remove missing Cv obs timesteps (remove from all)
    print(' remove missing Cv obs timesteps (remove from all)')
    for c in range(0, 1):
        um_data[varlist_um[c]][nanind, :] = np.nan
        misc_data[varlist_um[c]][nanind, :] = np.nan
        ra2t_data[varlist_um[c]][nanind, :] = np.nan
    ### remove missing water content obs timestep (only remove from water contents)
    for c in range(1, 3):
        um_data[varlist_um[c]][wcind, :] = np.nan
        misc_data[varlist_um[c]][wcind, :] = np.nan
        ra2t_data[varlist_um[c]][wcind, :] = np.nan


    # ### remove zeroed water content on obs timestep (only remove from water contents)
    # for c in range(1, 3):
    #     obs_data[varlist_obs[c]][wc0ind, :] = np.nan
    #     um_data[varlist_um[c]][wc0ind, :] = np.nan
    #     ifs_data[varlist_ifs[c]][wc0ind, :] = np.nan
    #     misc_data[varlist_um[c]][wc0ind, :] = np.nan
    #     ra2t_data[varlist_um[c]][wc0ind, :] = np.nan
    # ## remove missing lwpo obs timestep (only remove from water contents)
    # for c in range(1, 3):
    #     um_data[varlist_um[c]][lwpind, :] = np.nan
    #     ifs_data[varlist_ifs[c]][lwpind, :] = np.nan
    #     misc_data[varlist_um[c]][lwpind, :] = np.nan
    #     ra2t_data[varlist_um[c]][lwpind, :] = np.nan
    #
    ## lwp only 1d
    um_data['model_lwp'][lwpind] = np.nan
    misc_data['model_lwp'][lwpind] = np.nan
    ra2t_data['model_lwp'][lwpind] = np.nan


    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '25_u-cc568': label1 = 'UM_RA2M'
    if out_dir1[:10] == '24_u-cc324': label1 = 'UM_RA2T_' + out_dir1[-4:-1]
    if out_dir1[:10] == '23_u-cc278': label1 = 'UM_CASIM-100'

    label2 = 'undefined_label'
    if out_dir2[:10] == '25_u-cc568': label2 = 'UM_RA2M'
    if out_dir2[:10] == '24_u-cc324': label2 = 'UM_RA2T_' + out_dir2[-4:-1]
    if out_dir2[:10] == '23_u-cc278': label2 = 'UM_CASIM-100'

    label3 = 'undefined_label'
    if out_dir3 == 'OUT_25H/': label3 = 'ECMWF_IFS'
    if out_dir3[:10] == '25_u-cc568': label3 = 'UM_RA2M'
    if out_dir3[:10] == '24_u-cc324': label3 = 'UM_RA2T_' + out_dir3[-4:-1]
    if out_dir3[:10] == '23_u-cc278': label3 = 'UM_CASIM-100'

    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    ### model/measurement data
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    # np.save('working_dataObs', obs['sondes'])
    #
    # ### cloudnet
    # np.save('working_um_data', um_data)
    # np.save('working_ifs_data', ifs_data)
    # if cn_misc_flag != -1: np.save('working_misc_data', misc_data)

###################################################################################################################
###################################################################################################################
################################################ FIGURES ##########################################################
###################################################################################################################
###################################################################################################################

    # -------------------------------------------------------------
    # Cloudnet plot: Plot Cv statistics from drift period
    # -------------------------------------------------------------
    # figure = plot_CvProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs, obs_switch)
    # figure = plot_lwcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)
    # figure = plot_iwcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)
    # figure = plot_twcProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)

    # -------------------------------------------------------------
    # Cloudnet plot: Plot contour timeseries
    # -------------------------------------------------------------
    figure = plot_CvTimeseries(um_data, misc_data, ra2t_data, obs_data, dates,plots_out_dir)# ,monc_data=monc_data)
    #figure = plot_LWCTimeseries(um_data, misc_data, ra2t_data,obs_data, plots_out_dir, dates, obs_switch)
    #figure = plot_IWCTimeseries(um_data, misc_data, ra2t_data,obs_data, plots_out_dir, dates, obs_switch)
    #figure = plot_TWCTimeseries(um_data, misc_data, ra2t_data, obs_data, plots_out_dir, dates, obs_switch,nanind,wcind)
    # figure = plot_TWCTesting(um_data, ifs_data, misc_data, obs_data, data1, data2, data3, obs, month_flag, missing_files, doy)

    # -------------------------------------------------------------
    # Model plots
    # -------------------------------------------------------------
    # figure = plot_line_TSa(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_RadiosondesThetaE(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)

    # -------------------------------------------------------------
    # plot LWP timeseries with missing files accounted for
    # -------------------------------------------------------------
    # if obs_switch == 'RADAR': lwp = []
    # figure = plot_LWP(um_data, ifs_data, misc_data, ra2t_data, obs_data, obs, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)#, lwp) #, lon, lat):


    # -------------------------------------------------------------
    # plot bivariate distributions
    # -------------------------------------------------------------
    #figure = plot_BiVAR(um_data, ifs_data, misc_data, ra2t_data, obs_data, obs, month_flag, missing_files, cn_um_out_dir, doy, obs_switch, data1, data2, data3, data4, nanind, wcind)

    # -------------------------------------------------------------
    # make obs comparison fig between um and ifs grids
    # -------------------------------------------------------------
    # figure = plot_ObsGridComparison(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)

    # -------------------------------------------------------------
    # plot cloudnet split season figures with missing files accounted for
    # -------------------------------------------------------------
    # figure = plot_CvProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy)
    # figure = plot_lwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)
    # figure = plot_iwcProfiles_SplitSeason(um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)
    # figure = plot_profiles_SplitSeason(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch)

    # -------------------------------------------------------------
    # look closer at specific periods
    # -------------------------------------------------------------
    # figure = period_Selection(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs_switch, obs, data1, data2, data3, data4, nanind, wcind)

    # -------------------------------------------------------------
    # look closer at biases
    # -------------------------------------------------------------
    # figure = plot_BiasCorrelation(obs_data, um_data, misc_data, ifs_data, ra2t_data, doy, obs_switch)

    # -------------------------------------------------------------
    # cloud properties scaled by BL depth
    # -------------------------------------------------------------
    # -------------------------------------------------------------
    # Identify inversions in model / radiosonde data
    # -------------------------------------------------------------
    #################################################################
    ## load calculated model inversion heights
    #################################################################
    # print ('Load calculated model inversion heights (JV algorithm)...')
    # data1['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_RA2M_inversion_results.mat')
    # data2['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/UM_CASIM-100_inversion_results.mat')
    # data3['inversions'] = readMatlabStruct(obs_root_dir + 'radiosondes/ECMWF_IFS_inversion_results.mat')

    # print ('Load calculated model inversion heights (GY algorithm)...')
    # obs['inversions']['thetaE'] = np.load(um_root_dir[:-5] + 'obs_inversions_v2.npy').item()
    # data1['inversions'] = np.load(um_root_dir[:-5] + 'um_ra2m_inversions_v2.npy').item()
    # data2['inversions'] = np.load(um_root_dir[:-5] + 'um_casim-100_inversions_v2.npy').item()
    # data3['inversions'] = np.load(um_root_dir[:-5] + 'ecmwf_ifs_inversions_v2.npy').item()

    # -------------------------------------------------------------
    ### use IFS named directory to allocate variable to plot
    # -------------------------------------------------------------
    # if cn_ifs_out_dir[0] == 'cloud-fraction-ecmwf-grid/2018/': var = 'Cv'
    # if cn_ifs_out_dir[0] == 'lwc-scaled-ecmwf-grid/2018/': var = 'lwc'
    # if cn_ifs_out_dir[0] == 'iwc-Z-T-ecmwf-grid/2018/': var = 'iwc'
    # var = 'lwc'

    # obs_data = interpCloudnet(obs_data, month_flag, missing_files, doy)
    # figure = plot_scaledBL_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3, var)
    # figure = plot_scaledBLCv_thetaE(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_scaledBLCv_JVInv(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_scaledBLlwc(data1, data2, data3, um_data, ifs_data, misc_data, obs_data, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)


    # -------------------------------------------------------------
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    ### model/measurement data
    # np.save('working_data1', data1)
    # np.save('working_data2', data2)
    # np.save('working_data3', data3)
    # np.save('working_data4', data4)
    # np.save('working_dataObs', obs['inversions'])
    #
    # ### cloudnet
    # np.save('working_um_data', um_data)
    # np.save('working_ifs_data', ifs_data)
    # if cn_misc_flag != -1: np.save('working_misc_data', misc_data)
    # np.save('working_ra2t_data', ra2t_data)
    # np.save('working_obsV6_data', obs_data)

    # print (data1['height'][:])
    # print (data1['height'][data1['height'] <= 1e3])

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
