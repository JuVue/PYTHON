### SCRIPT TO READ IN UM, IFS, and UM-CASIM model data
### loads in the Cloudnet data and the raw model data for the cloudnet figures

from __future__ import print_function
import time
import datetime as dtime

import pandas as pd
from netCDF4 import Dataset
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import matplotlib.ticker as ticker
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,LogNorm
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
from physFuncts import calcThetaE, calcThetaVL,calcAirDensity, calcT
from manipFuncts import int2list
# from conversionFuncts import reGrid_Sondes


def plot_CvTimeseries(um_data,  obs_data, label,outstr, plots_out_dir,dates,  **args):
    numsp=len(um_data)+1
    if bool(args):
        monc_data=args[list(args.keys())[0]]
        numsp += len(monc_data)
        if len(args)>1:
            for n in range(1,len(args)):
                if list(args.keys())[n] == 'mlabel':
                    mlabel = args[list(args.keys())[n]]
                elif list(args.keys())[n] == 'moutdir':
                    moutdir= args[list(args.keys())[n]]

    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting Cv timeseries:')
    print ('')
    print('plotting' , numsp, 'subplots')

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:20, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)

    #####Plot Cv###############################################
    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(numsp,1,1)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['Cv']),
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
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 50)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('C$_{V}$')
    # plt.colorbar()
    for m in range(0,len(um_data)):
        plt.subplot(numsp,1,m+2)
        ax = plt.gca()
        plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_Cv_filtered']),
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
        if m == numsp:
            plt.xlabel('Date')
        nans = ax.get_ylim()
        ax2 = ax.twinx()
        ax2.set_ylabel(label[m], rotation = 270, labelpad = 27)
        ax2.set_yticks([])

    if bool(args):
        for m in range(0,len(monc_data)):
            plt.subplot(numsp,1,numsp-len(monc_data)+1+m)
            ax = plt.gca()
            plt.contourf(monc_data[m]['time2']/60/60, np.squeeze(monc_data[m]['z'][:]), np.transpose(monc_data[m]['total_cloud_fraction']),
                np.arange(0,1.1,0.1),
                cmap = newcmp,
                zorder = 1)
            # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
            # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
            plt.ylabel('Z [km]')
            plt.ylim(ylims)
            plt.yticks(yticks)
            ax.set_yticklabels(ytlabels)
            #plt.xlim(monc_data['time2']/60/60])
            xticks=np.arange(np.floor(monc_data[m]['time2'][0]/60/60),np.ceil(monc_data[m]['time2'][-1]/60/60)+1,2,dtype=int)
            xlabs=[x*100 for x in xticks]
            xlabs=["{0:-04d}".format(t) for t in xlabs]
            plt.xticks(xticks)
            ax.set_xticklabels(xlabs)
            plt.xlabel('Time (UTC)')
            nans = ax.get_ylim()
            ax2 = ax.twinx()
            embed()
            ax2.set_ylabel(mlabel[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if bool(args):
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_'.join(moutstr) + '_CvTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_CvTimeseries.png'
    print(fileout)
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_LWCTimeseries(um_data, obs_data, label,outstr,plots_out_dir, dates,obs_switch, **args): #, lon, lat):

    numsp=len(um_data)+1
    if bool(args):
        monc_data=args[list(args.keys())[0]]
        numsp += 1

    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ytlabels=yticks/1e3

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:1, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)

    obs_data['lwc'][obs_data['lwc'] <= 0] = np.nan
    obs_data['lwc_adiabatic'][obs_data['lwc_adiabatic'] <= 0] = np.nan
    for m in range(0,len(um_data)):
            um_data[m]['model_lwc'][um_data[m]['model_lwc'] <= 0.0] = np.nan

    if bool(args):
        monc_data['model_lwc']= monc_data['liquid_mmr_mean']*monc_data['rho']
        monc_data['model_lwc'][monc_data['model_lwc'] <= 0.0] = np.nan

    print ('******')
    print ('')
    print ('Plotting LWC timeseries for whole drift period:')
    print ('')

    cmax=0.3
    clev=np.arange(0.0,0.45,0.05)


    #####PlotLwc###############################################
    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)


    plt.subplot(numsp,1,1)
    ax = plt.gca()
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['lwc_adiabatic'])*1e3,
        levels=clev,cmap = newcmp)
    nans = ax.get_ylim()
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    plt.xlim([dates[0], dates[1]])
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    #nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 50)
    ax2.set_yticks([])
    #plt.title('Obs-' + obs_switch + 'grid')
    cbaxes = fig.add_axes([0.225, 0.95, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('LWC [g m$^{-3}$]')

    for m in range(0,len(um_data)):
        plt.subplot(numsp,1,m+2)
        ax = plt.gca()
        plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_lwc'])*1e3,
            levels=clev,cmap = newcmp)
        nans = ax.get_ylim()
        plt.ylabel('Z [km]')
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax.set_yticklabels(ytlabels)
        plt.xlim([dates[0], dates[1]])
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
        ax2 = ax.twinx()
        ax2.set_ylabel('UM_CASIM', rotation = 270, labelpad = 17)
        ax2.set_yticks([])
        # plt.colorbar()
        if m== numsp:
            plt.xlabel('Date')


    if bool(args):
        plt.subplot(numsp,1,numsp)
        ax = plt.gca()
        # ax.set_facecolor('aliceblue')
        plt.contourf(monc_data['time2']/60/60, np.squeeze(monc_data['z'][:]), np.transpose(monc_data['model_lwc'])*1e3,
        levels=clev,cmap = newcmp)
        plt.ylabel('Z [km]')
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax.set_yticklabels(ytlabels)
        xticks=np.arange(np.floor(monc_data['time2'][0]/60/60),np.ceil(monc_data['time2'][-1]/60/60)+1,2,dtype=int)
        xlabs=[x*100 for x in xticks]
        xlabs=["{0:-04d}".format(t) for t in xlabs]
        plt.xticks(xticks)
        ax.set_xticklabels(xlabs)
        plt.xlabel('Time (UTC)')
        ax2 = ax.twinx()
        ax2.set_ylabel('MONC', rotation = 270, labelpad = 17)
        ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if bool(args):
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_'  + '_'.join(outstr) + '_MONC_LWCTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_'  + '_'.join(outstr) + '_LWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_IWCTimeseries(um_data, obs_data,label,outstr, plots_out_dir, dates, obs_switch,**args): #, lon, lat):

    numsp=len(um_data)+1
    if bool(args):
        monc_data=args[list(args.keys())[0]]
        numsp += 1


    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ytlabels=yticks/1e3

    obs_data['iwc'][obs_data['iwc'] <= 0.0] = np.nan
    for m in range(0,len(um_data)):
        um_data[m]['model_iwc_filtered'][um_data[m]['model_iwc_filtered'] <= 0.0] = np.nan

    if bool(args):
        monc_data['model_iwc']= (monc_data['ice_mmr_mean']+monc_data['graupel_mmr_mean']+monc_data['snow_mmr_mean'])*monc_data['rho']
        monc_data['model_iwc'][monc_data['model_iwc'] <= 0.0] = np.nan

    embed()
    print ('******')
    print ('')
    print ('Plotting IWC timeseries for whole drift period:')
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

    viridis = mpl_cm.get_cmap('viridis', 256) # nice colormap purple to yellow
    newcolors = viridis(np.linspace(0, 1, 256)) #assgin new colormap with viridis colors
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:1, :] = greyclr   # make first 20 colors greyclr
    newcmp = ListedColormap(newcolors)

    #####PlotIwc###############################################
    cmax = 0.05
    clev=np.arange(0.0,0.05,0.001)
    clev=[1e-5,0.5e-4,1e-4,0.5e-3, 1e-3, 0.5e-2,1e-2,0.5e-1,1e-1]

    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(numsp,1,1)
    ax = plt.gca()
    if obs_switch == 'RADAR':
        img = plt.pcolor(obs_data['time'], obs_data['height'], np.transpose(obs_data['iwc'])*1e3,
            cmap=newcmp,vmin = 0.0, vmax = cmax)
            #cmap = newcmp)
    else:
        img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data['iwc'])*1e3,
            levels=clev, norm = LogNorm(),cmap = newcmp)
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
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 50)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.96, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('IWC')

    for m in range(0,len(um_data)):
        plt.subplot(numsp,1,m+2)
        ax = plt.gca()
        plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_iwc_filtered'])*1e3,
            levels=clev,norm = LogNorm(),cmap = newcmp)
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
        ax2.set_ylabel(label[m], rotation = 270, labelpad = 27)
        ax2.set_yticks([])

        if m == numsp:
            plt.xlabel('Date')


    if bool(args):
        plt.subplot(numsp,1,numsp)
        ax = plt.gca()
        # ax.set_facecolor('aliceblue')
        plt.contourf(monc_data['time2']/60/60, np.squeeze(monc_data['z'][:]), np.transpose(monc_data['model_iwc'])*1e3,
        levels=clev,norm = LogNorm(),cmap = newcmp)
#        cmap=newcmp,vmin = 0.0, vmax = cmax)
        plt.ylabel('Z [km]')
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax.set_yticklabels(ytlabels)
        xticks=np.arange(np.floor(monc_data['time2'][0]/60/60),np.ceil(monc_data['time2'][-1]/60/60)+1,2,dtype=int)
        xlabs=[x*100 for x in xticks]
        xlabs=["{0:-04d}".format(t) for t in xlabs]
        plt.xticks(xticks)
        ax.set_xticklabels(xlabs)
        plt.xlabel('Time (UTC)')
        ax2 = ax.twinx()
        ax2.set_ylabel('MONC', rotation = 270, labelpad = 27)
        ax2.set_yticks([])


    dstr=datenum2date(dates[1])
    if bool(args):
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_MONC_IWCTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_TWCTimeseries(um_data,  obs_data,label,outstr, plots_out_dir, dates, obs_switch, **args):

    numsp=len(um_data)+1
    if bool(args):
        monc_data=args[list(args.keys())[0]]
        numsp += 1

    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting TWC timeseries for whole drift period:')
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    for m in range(0,len(um_data)):
        um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc_filtered']

    if bool(args):
        monc_data['model_iwc']= (monc_data['ice_mmr_mean']+monc_data['graupel_mmr_mean']+monc_data['snow_mmr_mean'])*monc_data['rho']
        monc_data['model_lwc']= monc_data['liquid_mmr_mean']*monc_data['rho']
        monc_data['model_twc'] = monc_data['model_lwc'] +monc_data['model_iwc']

    twc0 = np.transpose(obs_data['twc'])*1e3

    cmax = 0.3

    viridis = mpl_cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    greyclr = np.array([0.1, 0.1, 0.1, 0.1])
    newcolors[:10, :] = greyclr
    newcmp = ListedColormap(newcolors)

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    #####Plot Twc###############################################
    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(numsp,1,1)
    ax = plt.gca()
    # ax.set_facecolor('aliceblue')
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), twc0,
        #levels=clev,cmap = newcmp)
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
    #nans = ax.get_ylim()
    ax2 = ax.twinx()
    ax2.set_ylabel('Measurements \n (1 hour sampling)', rotation = 270, labelpad = 50)
    ax2.set_yticks([])
    #plt.title('Obs-' + obs_switch + 'grid')
    cbaxes = fig.add_axes([0.225, 0.95, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')
    plt.title('TWC [g m$^{-3}$]')

    for m in range(0,len(um_data)):
        plt.subplot(numsp,1,m+2)
        ax = plt.gca()
        plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_twc'])*1e3,
            levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
            cmap = newcmp)
        nans = ax.get_ylim()
        plt.ylabel('Z [km]')
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax.set_yticklabels(ytlabels)
        plt.xlim([dates[0], dates[1]])
        ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
        ax2 = ax.twinx()
        ax2.set_ylabel(label[m], rotation = 270, labelpad = 27)
        ax2.set_yticks([])

        if m== numsp:
            plt.xlabel('Time UTC')

    if bool(args):
        plt.subplot(numsp,1,numsp)
        ax = plt.gca()
        # ax.set_facecolor('aliceblue')
        plt.contourf(monc_data['time2']/60/60, np.squeeze(monc_data['z'][:]), np.transpose(monc_data['model_twc'])*1e3,
        levels=[1e-4, 1e-3, 1e-2, 1e-1, 1e0], norm = LogNorm(),
        cmap = newcmp)
        # )
        plt.ylabel('Z [km]')
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax.set_yticklabels(ytlabels)
        xticks=np.arange(np.floor(monc_data['time2'][0]/60/60),np.ceil(monc_data['time2'][-1]/60/60)+1,2,dtype=int)
        xlabs=[x*100 for x in xticks]
        xlabs=["{0:-04d}".format(t) for t in xlabs]
        plt.xticks(xticks)
        ax.set_xticklabels(xlabs)
        plt.xlabel('Time (UTC)')
        ax2 = ax.twinx()
        ax2.set_ylabel('MONC', rotation = 270, labelpad = 27)
        ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if bool(args):
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_MONC_TWCTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid' + '_'.join(outstr) + '_TWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()
    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')


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
    print('nanmask')
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

def setFlags(obs_data, um_data,obs_var_list, um_var_list):

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
    for m in range(0,len(um_data)):
        for c in range(0,3):
            for j in range(0,len(um_var_list[c])):
                um_data[m][um_var_list[c][j]][um_data[m][um_var_list[c][j]] == -999] = np.nan
                um_data[m][um_var_list[c][j]][um_data[m][um_var_list[c][j]] < 0] = 0.0

    return obs_data, um_data

def removeSpinUp(monc_data,monc_spin):
    for m in range(0,len(monc_data)):
        monc_var_list = list(monc_data[m].keys())
        monc_var_list.remove('time1')
        monc_var_list.remove('time2')

        id1 = np.squeeze(np.argwhere(monc_data[m]['time1']<=monc_spin)) #1D data
        id2 = np.squeeze(np.argwhere(monc_data[m]['time2']<=monc_spin))
        for j in range(0,len(monc_var_list)):
            if any(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time1'])):
                monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id1,0)
            elif any(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time2'])):
                tmp2=np.argwhere(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time2']))
                if tmp2 == 0:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,0)
                elif tmp2 == 1:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,1)
                elif tmp2 == 2:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,2)
                elif tmp2 == 3:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,3)
            else:
                print(monc_var_list[j], ': not time depend')
        monc_data[m]['time1']=np.delete(monc_data[m]['time1'],id1,0)
        monc_data[m]['time2']=np.delete(monc_data[m]['time2'],id2,0)
        monc_data[m]['time1']=monc_data[m]['time1']-monc_data[m]['time1'][0]
        monc_data[m]['time2']=monc_data[m]['time2']-monc_data[m]['time2'][0]

    return monc_data
################################################################################
################################################################################
def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    names = ['20180912_oden_','20180913_oden_']

    sdate = dtime.datetime.strptime('2018091300','%Y%m%d%H')
    edate = dtime.datetime.strptime('2018091314','%Y%m%d%H')

    dates = [date2datenum(sdate),date2datenum(edate)]
    moccha_missing_files = ['20180813_oden_','20180910_oden_']   ### cloud radar not working    #,'20180914_oden_'

    #---- MONC SPIN UP TIME
    monc_spin = 6 *60 *60

    ### Set output directory for plots
    plots_out_dir='/nfs/a96/MOCCHA/working/jutta/plots/CaseStudies/ModelComparison/'
    if not os.path.exists(plots_out_dir):
        os.makedirs(plots_out_dir)

    ### Choose observations vertical gridding used in Cloudnet processing (UM/IFS/RADAR)
    obs_switch = 'UM' #RADAR
    #/nfs/a96/MOCCHA/working/gillian/UM/DATA
    um_root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
    cn_obs_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/Halo/measurements_V7/'
    cn_um_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/metum/V7/'
    monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'

    ### -----------------------------------------------------------------
    ### CHOOSE UM RUNS - MODEL DATA
    out_dir = ['23_u-cc278_RA1M_CASIM/',
               '24_u-cc324_RA2T_CON/',
               '25_u-cc568_RA2M_CON/']

    ### CHOOSE MONC RUNS
    m_out_dir = ['4_control_20180913T0000Z_Wsub-1.5/',
                '5_control_20180913T0000Z_Wsub-1.5_Fletcher/']
    monc_filename=[]
    for m in range(0, len(m_out_dir)):
        monc_filename.append(monc_root_dir + m_out_dir[m] + 'moccha_casim_dg_72000.nc')

    ### -----------------------------------------------------------------
    ### CHOSEN RUN - CLOUDNET DATA
    cn_um_out_dir = ['cloud-fraction-metum-grid/2018/',
                     'lwc-scaled-metum-grid/2018/',
                     'iwc-Z-T-metum-grid/2018/']
    cn_obs_out_dir = ['cloud-fraction-metum-grid/2018/qual70/',
                        'lwc-adiabatic-metum-grid/2018/',
                        'iwc-Z-T-metum-grid/2018/']


    # # -------------------------------------------------------------
    # # Load cloudnet data
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin nc read in')
    print (' ')

    for i in range(0,len(names)):
        dstr=names[i][0:8]
        datenum = date2datenum(dtime.datetime.strptime(dstr,'%Y%m%d'))
        cn_filename_um=[]
        for m in range(0,len(out_dir)):
            cn_filename_um.append([cn_um_dir + out_dir[m] + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + out_dir[m] + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + out_dir[m] + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc'])
        cn_filename_obs = [cn_obs_dir + cn_obs_out_dir[0] + names[i] + cn_obs_out_dir[0][:-13] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[1] + names[i] + cn_obs_out_dir[1][:-6] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[2] + names[i] + cn_obs_out_dir[2][:-6] + '.nc']
        print ('')

    ### --------------------------------------------------------------------
    ###     READ IN ALL CLOUDNET FILES: reinitialise diagnostic dictionaries
    ### --------------------------------------------------------------------
        print ('')
        print ('*************** LOAD DAY:', i)
        print ('Loading multiple diagnostics from cloudnet results:')
        cn_nc0 = {} # observations
        for c in range(0,3):
            cn_nc0[c] = Dataset(cn_filename_obs[c],'r')
        cn_nc1={}
        for m in range(0,len(out_dir)): #UM model data
            cn_nc1[m]={}
            for c in range(0,3):
                cn_nc1[m][c] = Dataset(cn_filename_um[m][c],'r')

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
        ### --------------------------------------------------------------------
        ### create arrays for all cloudnet data
        ### --------------------------------------------------------------------
        if i == 0:
            obs_data = {}  ### OBS data
            if obs_switch == 'RADAR':
                time_obs = datenum + np.float64((cn_nc0[1].variables['time'][:])/24.0)
            else:
                time_obs = datenum + np.float64((cn_nc0[0].variables['time'][:])/24.0)

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
        else:
            if obs_switch == 'RADAR':
                time_obs = np.append(time_obs,  datenum +  np.float64((cn_nc0[1].variables['time'][:])/24.0))
            else:
                time_obs = np.append(time_obs,  datenum + np.float64((cn_nc0[0].variables['time'][:])/24.0))
            for c in range(0,3):
                ### append rest of obs data
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

        print ('   observation data loaded!')

        if i == 0:
            um_data={}
            time_um={}

        for m in range(0,len(out_dir)): #UM model data
            if i == 0:
                um_data[m] = {}             ### initialise cloudnet data dictionaries
                time_um[m] = {}
                ### create time arrays for all cloudnet data
                time_um[m] = datenum + np.float64((cn_nc1[m][0].variables['time'][:])/24.0)
                ### loop over each Cloudnet class
                for c in range(0,3):
                    ### load in initial UM data
                    for j in range(0,len(um_var_list[c])):
                        if np.ndim(cn_nc1[m][c].variables[um_var_list[c][j]]) == 1:  # 1d timeseries only
                            um_data[m][um_var_list[c][j]] = cn_nc1[m][c].variables[um_var_list[c][j]][:]
                        else:                                   # 2d column um_data
                            um_data[m][um_var_list[c][j]] = cn_nc1[m][c].variables[um_var_list[c][j]][:]
            ### fill arrays with remaining data
            else:
                time_um[m] = np.append(time_um[m], datenum + np.float64((cn_nc1[m][0].variables['time'][:])/24.0))
                ### loop over all Cloudnet classes
                for c in range(0,3):
                    ### append rest of UM data
                    for j in range(0,len(um_var_list[c])):
                        if np.ndim(cn_nc1[m][c].variables[um_var_list[c][j]]) == 1:
                            um_data[m][um_var_list[c][j]] = np.append(um_data[m][um_var_list[c][j]],cn_nc1[m][c].variables[um_var_list[c][j]][:])
                        else:
                            um_data[m][um_var_list[c][j]] = np.append(um_data[m][um_var_list[c][j]],cn_nc1[m][c].variables[um_var_list[c][j]][:],0)

        ### --------------------------------------------------------------------
        ### Close cloudnet netCDF files
        ### --------------------------------------------------------------------
        for c in range(0,3): cn_nc0[c].close()
        for m in range(0,len(out_dir)):
             for c in range(0,3): cn_nc1[m][c].close()

        print ('    UM data loaded!')
        print ('')
    #################################################################
    ## save time to dictionaries now we're not looping over all diags anymore
    #################################################################
    for m in range(0,len(out_dir)):
        um_data[m]['time'] = time_um[m]
    obs_data['time'] = time_obs

    #################################################################
    ###     READ IN MONC DATA
    #################################################################
    print ('Loading monc results:')
    print ('')
    ###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)
    monc_var_list =[['time_series_2_60','time_series_20_600' ,'z','rho', 'LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
                    ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
                    'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
                    'graupel_mmr_mean']]
                #    ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp'],
                #    ['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]

    ncm = {}
    monc_data = {}
    for m in range(0, len(m_out_dir)):
        ncm = Dataset(monc_filename[m],'r')
        monc_data[m]={}
        for c in range(0,len(monc_var_list)):
            for j in range(0,len(monc_var_list[c])):
                monc_data[m][monc_var_list[c][j]] = ncm.variables[monc_var_list[c][j]][:]

        monc_data[m]['time2']=monc_data[m]['time_series_20_600'] #2d data
        monc_data[m]['time1']=monc_data[m]['time_series_2_60'] #1d data
        monc_data[m].pop('time_series_2_60')
        monc_data[m].pop('time_series_20_600')

    print (' Monc data Loaded!')

    ## remove spin up time from monc data
    monc_data=removeSpinUp(monc_data,monc_spin)
    ##################################################################################################################################
    ## -------------------------------------------------------------
    ## maximise obs data available and build mask for available data
    ## -------------------------------------------------------------
    print ('*******************')
    print('setting missing data to nan and interpolate missing obs')
    obs_data, um_data = setFlags(obs_data, um_data, obs_var_list, um_var_list)
    obs_data = interpCloudnet(obs_data)
    nanind, nanmask, wcind, wc0ind, lwpind = buildNaNMask(obs_data)
    #nanindadv, nanmaskadv, wcindadv, wc0indadv, lwpindadv = buildNaNMaskadv(obs_data)

    varlist_obs = ['Cv', 'lwc_adiabatic', 'iwc', 'lwp']
    varlist_um = ['model_Cv_filtered', 'model_lwc', 'model_iwc_filtered', 'model_lwp']
    # print(um_data.keys())

    ## remove missing Cv obs timesteps (remove from all)
    print(' remove missing Cv obs timesteps')
    for c in range(0, 1):
        for m in range(0, len(out_dir)):
            um_data[m][varlist_um[c]][nanind, :] = np.nan

    ### remove missing water content obs timestep (only remove from water contents)
    for c in range(1, 3):
        for m in range(0, len(out_dir)):
            um_data[m][varlist_um[c]][wcind, :] = np.nan

    ## lwp only 1d
    for m in range(0, len(out_dir)):
        um_data[m]['model_lwp'][lwpind] = np.nan

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label=[]
    outstr=[]
    for m in range(0, len(out_dir)):
        if out_dir[m][:10] == '24_u-cc324':
            label.append('UM_RA2T_' + out_dir[m][-4:-1])
            outstr.append('RA2T')
        elif out_dir[m][:10] == '25_u-cc568':
            label.append('UM_RA2M')
            outstr.append('RA2M')
        elif out_dir[m][:10] == '23_u-cc278':
            label.append('UM_CASIM-100')
            outstr.append('CASIM100')
        elif out_dir[m][:10] == '26_u-cd847':
            label.append('UM_CASIM-AP')
            outstr.append('CASIM-AP')
        elif out_dir[m][:10] == '27_u-ce112':
            label.append('UM_CASIM-AP-PasProc')
            outstr.append('CASIM-AP-PasProc')
        else:
            label.append('undefined_label')
            outstr.append('')

    mlabel=[]
    moutstr=[]
    for m in range(0, len(m_out_dir)):
        if m_out_dir[m][:1] == '4':
            mlabel.append('MONC Wsub1.5' + out_dir[m][-4:-1])
            moutstr.append('MONC-wsub')
        elif out_dir[m][:1] == '5':
            mlabel.append('MONC Wsub1.5 \n Fletcher')
            moutstr.append('MONC-wsub-fle')
        else:
            label.append('undefined_label')
            outstr.append('')

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
    figure = plot_CvTimeseries(um_data, obs_data, label, outstr, plots_out_dir, dates, x=monc_data,mlabel=mlabel,moutsr=moutstr)
    #figure = plot_LWCTimeseries(um_data,obs_data, label, outstr, plots_out_dir, dates, obs_switch,x=monc_data)
    #figure = plot_TWCTimeseries(um_data, obs_data, label,outstr, plots_out_dir, dates, obs_switch,x=monc_data)
    #figure = plot_IWCTimeseries(um_data, obs_data, label, outstr,plots_out_dir, dates, obs_switch,x=monc_data)
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
