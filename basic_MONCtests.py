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
import glob

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


def plot_LWCTimeseries(obs_data,lwcvar,lwcstr, plots_out_dir, dates, **args): #, lon, lat):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
                pmonc=True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
                pum=True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    # ylims=[0,2.5]
    # yticks=np.arange(0,2.5e3,0.5e3)
    # ylims=[0,6]
    # yticks=np.arange(0,6e3,0.5e3)
    # ytlabels=yticks/1e3

    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
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
    obs_data['lwc_adiabatic_inc_nolwp'][obs_data['lwc_adiabatic_inc_nolwp'] <= 0] = np.nan
    if pum==True:
        for m in range(0,len(um_data)):
                um_data[m]['model_lwc'][um_data[m]['model_lwc'] <= 0.0] = np.nan

    if pmonc==True:
        lwc_tvar=[]
        lwc_zvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            monc_data[m]['model_lwc']=monc_data[m]['lwc_tot_mean'].copy()
            monc_data[m]['model_lwc'][monc_data[m]['model_lwc'] <= 0.0] = np.nan
            #lwc_tvar=monc_data[m]['tvar']['liquid_mmr_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_tvar+=[monc_data[m]['tvar']['lwc_tot_mean']]
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]
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
    img = plt.contourf(obs_data['time'], np.squeeze(obs_data['height'][0,:]), np.transpose(obs_data[lwcvar])*1e3,
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
    if pum==True:
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
            ax2.set_ylabel(label[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])
            # plt.colorbar()
            if m== numsp:
                plt.xlabel('Date')

    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.subplot(numsp,1,numsp-len(monc_data)+1+m)
            ax = plt.gca()
            # ax.set_facecolor('aliceblue')
            plt.contourf(monc_data[m][lwc_tvar[m]]/60/60, np.squeeze(monc_data[m][lwc_zvar[m]][:]), np.transpose(monc_data[m]['model_lwc'])*1e3,
            levels=clev,cmap = newcmp)
            plt.ylabel('Z [km]')
            plt.ylim(ylims)
            plt.yticks(yticks)
            ax.set_yticklabels(ytlabels)
            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
            # xticks=np.arange(np.floor(monc_data[m][lwc_tvar][0]/60/60),np.ceil(monc_data[m][lwc_tvar][-1]/60/60)+1,2,dtype=int)
            # xlabs=[x*100 for x in xticks]
            # xlabs=["{0:-04d}".format(t) for t in xlabs]
            # plt.xticks(xticks)
            # ax.set_xticklabels(xlabs)
            plt.xlabel('Time (UTC)')
            ax2 = ax.twinx()
            ax2.set_ylabel(mlabel[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if pmonc == True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) +'_' + '_'.join(moutstr) + '_LWCTimeseries'+ lwcstr + '.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWCTimeseries'+ lwcstr + '.png'
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_IWCTimeseries( obs_data, plots_out_dir, dates,**args): #, lon, lat):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
                pmonc=True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
                pum=True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    # ylims=[0,2.5]
    # yticks=np.arange(0,2.5e3,0.5e3)
    # ylims=[0,6]
    # yticks=np.arange(0,6e3,0.5e3)
    # ytlabels=yticks/1e3
    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    ytlabels=yticks/1e3

    obs_data['iwc'][obs_data['iwc'] <= 0.0] = np.nan
    if pum ==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_iwc_filtered'][um_data[m]['model_iwc_filtered'] <= 0.0] = np.nan

    if pmonc ==True :
        iwc_tvar=[]
        iwc_zvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            monc_data[m]['model_iwc']= (monc_data[m]['iwc_tot_mean'])
            monc_data[m]['model_iwc'][monc_data[m]['model_iwc'] <= 0.0] = np.nan
            #iwc_tvar=monc_data[m]['tvar']['ice_mmr_mean']
            #iwc_zvar=monc_data[m]['zvar']['ice_mmr_mean']
            iwc_tvar+=[monc_data[m]['tvar']['iwc_tot_mean']]
            iwc_zvar+=[monc_data[m]['zvar']['iwc_tot_mean']]

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
    clev=[1e-5,0.25e-4,0.5e-4,0.75e-4,
          1e-4,0.25e-3,0.5e-3,0.75e-3,
          1e-3,0.25e-2,0.5e-2,0.75e-2,
          1e-2,0.25e-1,0.5e-1,0.75e-1,1e-1]

    #clevs=1e-5,0.5e-3, 1e-3, 0.5e-2,1e-2,0.5e-1,1e-1]

    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)

    plt.subplot(numsp,1,1)
    ax = plt.gca()
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
    if pum ==True:
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


    if pmonc == True:
        for m in range(0,len(monc_data)):
            plt.subplot(numsp,1,numsp-len(monc_data)+1+m)
            ax = plt.gca()
            # ax.set_facecolor('aliceblue')
            plt.contourf(monc_data[m][iwc_tvar[m]]/60/60, np.squeeze(monc_data[m][iwc_zvar[m]][:]), np.transpose(monc_data[m]['model_iwc'])*1e3,
            levels=clev,norm = LogNorm(),cmap = newcmp)
    #        cmap=newcmp,vmin = 0.0, vmax = cmax)
            plt.ylabel('Z [km]')
            plt.ylim(ylims)
            plt.yticks(yticks)
            ax.set_yticklabels(ytlabels)
            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
            # xticks=np.arange(np.floor(monc_data[m][iwc_tvar][0]/60/60),np.ceil(monc_data[m][iwc_tvar][-1]/60/60)+1,2,dtype=int)
            # xlabs=[x*100 for x in xticks]
            # xlabs=["{0:-04d}".format(t) for t in xlabs]
            # plt.xticks(xticks)
            # ax.set_xticklabels(xlabs)
            plt.xlabel('Time (UTC)')
            ax2 = ax.twinx()
            ax2.set_ylabel(mlabel[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if pmonc ==True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) +'_' + '_'.join(moutstr) + '_IWCTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWCTimeseries.png'
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_TWCTimeseries(obs_data,twcvar,twcstr,plots_out_dir, dates,  **args):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
                pmonc=True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
                pum=True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    # ylims=[0,6]
    # yticks=np.arange(0,6e3,1e3)
    # ytlabels=yticks/1e3
    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting TWC timeseries for whole drift period:')
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
#    obs_data['twc'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc_filtered']

    if pmonc==True:
        twc_tvar=[]
        twc_zvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            twc_tvar+=[monc_data[m]['tvar']['twc_tot_mean']]
            twc_zvar+=[monc_data[m]['zvar']['twc_tot_mean']]

    twc0 = np.transpose(obs_data[twcvar])*1e3

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

    #clevs=[1e-5, 0.5e-4, 1e-4,0.5e-3, 1e-3,0.5e-2, 1e-2,0.5e-1, 1e-1,0.5e0, 1e0]
    clevs=[1e-4,0.25e-3,0.5e-3,0.75e-3,
            1e-3,0.25e-2,0.5e-2,0.75e-2,
            1e-2,0.25e-1,0.5e-1,0.75e-1,
            1e-1,0.25e-0,0.5e-0,0.75e-0,1e-0]
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
        levels=clevs, norm = LogNorm(),
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
    if pum==True:
        for m in range(0,len(um_data)):
            plt.subplot(numsp,1,m+2)
            ax = plt.gca()
            plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_twc'])*1e3,
                levels=clevs, norm = LogNorm(),
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

    if pmonc == True:
        for m in range(0,len(monc_data)):
            plt.subplot(numsp,1,numsp-len(monc_data)+1+m)
            ax = plt.gca()
            # ax.set_facecolor('aliceblue')
            plt.contourf(monc_data[m][twc_tvar[m]][:]/60/60, np.squeeze(monc_data[m][twc_zvar[m]][:]), np.transpose(monc_data[m]['model_twc'])*1e3,
            levels=clevs, norm = LogNorm(),
            cmap = newcmp)
            # )
            plt.ylabel('Z [km]')
            plt.ylim(ylims)
            plt.yticks(yticks)
            ax.set_yticklabels(ytlabels)
            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
            # xticks=np.arange(np.floor(monc_data[m][twc_tvar][0]/60/60),np.ceil(monc_data[m][twc_tvar][-1]/60/60)+1,2,dtype=int)
            # xlabs=[x*100 for x in xticks]
            # xlabs=["{0:-04d}".format(t) for t in xlabs]
            # plt.xticks(xticks)
            # ax.set_xticklabels(xlabs)
            plt.xlabel('Time (UTC)')
            ax2 = ax.twinx()
            ax2.set_ylabel(mlabel[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if pmonc==True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_TWCTimeseries' + twcstr + '.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_TWCTimeseries' + twcstr + '.png'
    plt.savefig(fileout)
    plt.close()
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')


def plot_lwcProfiles(obs_data,lwcvar,lwcstr, thresholding, plots_out_dir,dates, **args): #, lon, lat):

    obs_zorder = 1

    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                obs_zorder += len(monc_data)
                pmonc =True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                obs_zorder += len(um_data)
                pum =True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting LWC mean profiles:' + lwcvar)
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc_filtered']
    if pmonc==True:
        lwc_zvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]

    if thresholding == True:
        if pum==True:

            ####-------------------------------------------------------------------------
            ### from Michael's paper:
            ####    "we consider a grid point cloudy when cloud water exceeded
            ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
            ####    interpolation in between."
            ####-------------------------------------------------------------------------
            m=0 # use first um model run for height grid definition
            twc_thresh_um = np.zeros([np.size(um_data[m]['model_twc'],1)])
            ####-------------------------------------------------------------------------
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            um_lt1km = np.where(um_data[m]['height'][0,:]<=1e3)
            twc_thresh_um[um_lt1km] = 1e-6
            ####-------------------------------------------------------------------------
            ### next look at values above 4 km
            ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
            um_lt1km = np.where(um_data[m]['height'][0,:]>=4e3)
            twc_thresh_um[um_lt1km] = 1e-7

            ### find heights not yet assigned
            um_intZs = np.where(twc_thresh_um == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(um_data)):
                twc_thresh_um[um_intZs] = f(um_data[m]['height'][0,um_intZs].data)

            for t in range(0,np.size(um_data[m]['model_twc'],0)):
                for k in range(0,np.size(um_data[m]['model_twc'],1)):
                    if obs_data['twc'][t,k] < twc_thresh_um[k]:
                        obs_data['twc'][t,k] = np.nan
                        obs_data['lwc'][t,k] = np.nan
                    if obs_data['twc_ad'][t,k] < twc_thresh_um[k]:
                        obs_data['twc_ad'][t,k] = np.nan
                        obs_data['lwc_adiabatic'][t,k] = np.nan
                    if obs_data['twc_ad_nolwp'][t,k] < twc_thresh_um[k]:
                        obs_data['twc_ad_nolwp'][t,k] = np.nan
                        obs_data['lwc_adiabatic_inc_nolwp'][t,k] = np.nan
                    for m in range(0,len(um_data)):
                        if  um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                            um_data[m]['model_twc'][t,k] = np.nan
                            um_data[m]['model_lwc'][t,k] = np.nan
        if pmonc==True:
            m=0 # use first um model run for height grid definition
            twc_thresh_monc = np.zeros([np.size(monc_data[m]['model_twc'],1)])
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]<=1e3)
            twc_thresh_monc[monc_lt1km] = 1e-6
            monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]>=4e3)
            twc_thresh_monc[monc_lt1km] = 1e-7
            monc_intZs = np.where(twc_thresh_monc == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(monc_data)):
                twc_thresh_monc[monc_intZs] = f(monc_data[m][lwc_zvar[m]][monc_intZs].data)

            for t in range(0,np.size(monc_data[m]['model_twc'],0)):
                for k in range(0,np.size(monc_data[m]['model_twc'],1)):
                    for m in range(0,len(monc_data)):
                        if  monc_data[m]['model_twc'][t,k] < twc_thresh_monc[k]:
                            monc_data[m]['model_twc'][t,k] = np.nan
                            monc_data[m]['model_lwc'][t,k] = np.nan

            # ### plot profile of threshold as sanity check
            # plt.plot(twc_thresh_um, um_data[0]['height'][0,:])
            # plt.plot(twc_thresh_monc, monc_data[0]['z'][:]); plt.show()
            # print (twc_thresh_monc)

    ###----------------------------------------------------------------
    ###         Plot figure - Mean profiles
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
    plt.figure(figsize=(6,8.5))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)
    ### define axis instance
    ax1 = plt.gca()
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    #### ADIABATIC LWC
    plt.plot(np.nanmean(obs_data[lwcvar],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid'  + lwcstr, zorder = obs_zorder)
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[lwcvar],0)*1e3 - np.nanstd(obs_data[lwcvar],0)*1e3,
        np.nanmean(obs_data[lwcvar],0)*1e3 + np.nanstd(obs_data[lwcvar],0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.xlim([0,0.2])
    plt.plot(np.nanmean(obs_data[lwcvar],0)*1e3 - np.nanstd(obs_data[lwcvar],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(obs_data[lwcvar],0)*1e3 + np.nanstd(obs_data[lwcvar],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    if pum==True:
        for m in range(0,len(um_data)):
            ax1.fill_betweenx(np.nanmean(um_data[m]['height'],0),np.nanmean(um_data[m]['model_lwc'],0)*1e3 - np.nanstd(um_data[m]['model_lwc']*1e3,0),
                np.nanmean(um_data[m]['model_lwc'],0)*1e3 + np.nanstd(um_data[m]['model_lwc'],0)*1e3, color = fcols[m], alpha = 0.05)
            plt.plot(np.nanmean(um_data[m]['model_lwc'],0)*1e3 - np.nanstd(um_data[m]['model_lwc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color =lcols[m], linewidth = 0.5)
            plt.plot(np.nanmean(um_data[m]['model_lwc'],0)*1e3 + np.nanstd(um_data[m]['model_lwc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color = lcols[m], linewidth = 0.5)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            ax1.fill_betweenx(monc_data[m][lwc_zvar[m]],np.nanmean(monc_data[m]['model_lwc'],0)*1e3 - np.nanstd(monc_data[m]['model_lwc']*1e3,0),
                np.nanmean(monc_data[m]['model_lwc'],0)*1e3 + np.nanstd(monc_data[m]['model_lwc'],0)*1e3, color = fcolsmonc[m], alpha = 0.05)
            plt.plot(np.nanmean(monc_data[m]['model_lwc'],0)*1e3 - np.nanstd(monc_data[m]['model_lwc'],0)*1e3, monc_data[m][lwc_zvar[m]],
                '--', color =lcolsmonc[m], linewidth = 0.5)
            plt.plot(np.nanmean(monc_data[m]['model_lwc'],0)*1e3 + np.nanstd(monc_data[m]['model_lwc'],0)*1e3, monc_data[m][lwc_zvar[m]],
                '--', color = lcolsmonc[m], linewidth = 0.5)
    if pum==True:
        for m in range(0,len(um_data)):
            plt.plot(np.nanmean(um_data[m]['model_lwc'],0)*1e3,np.nanmean(um_data[m]['height'],0), color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.plot(np.nanmean(monc_data[m]['model_lwc'],0)*1e3,monc_data[m][lwc_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)

    plt.xlabel('Liquid water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    # plt.ylim([0,5e3])
    # plt.yticks(np.arange(0,5.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax1.set_yticklabels(ytlabels)
    plt.xlim([0,0.4])
    plt.xticks(np.arange(0,0.45,0.05))
    #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
    plt.legend()
    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_LWC-MTThresh' + lwcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWC-MTThresh'+ lwcstr + '.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_LWC' + lwcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWC'+ lwcstr + '.png'

    plt.savefig(fileout)

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')


def plot_iwcProfiles(obs_data, twcvar,twcstr, thresholding,plots_out_dir,dates, **args): #, lon, lat):

    obs_zorder = 1

    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                obs_zorder += len(monc_data)
                pmonc=True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                obs_zorder += len(um_data)
                pum=True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    ytlabels=yticks/1e3

    print ('******')
    print ('')
    print ('Plotting IWC mean profiles:')
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc_filtered']
    if pmonc==True:
        twc_zvar=[]
        for m in range(0,len(monc_data)):
            # monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            # monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            # monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            # twc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_iwc'] = monc_data[m]['iwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            twc_zvar+=[monc_data[m]['zvar']['twc_tot_mean']]

    if thresholding ==True:
        if pum==True:
            ####-------------------------------------------------------------------------
            ### from Michael's paper:
            ####    "we consider a grid point cloudy when cloud water exceeded
            ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
            ####    interpolation in between."
            ####-------------------------------------------------------------------------
            m=0 # use first um model run for height grid definition
            twc_thresh_um = np.zeros([np.size(um_data[m]['model_twc'],1)])
            ####-------------------------------------------------------------------------
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            um_lt1km = np.where(um_data[m]['height'][0,:]<=1e3)
            twc_thresh_um[um_lt1km] = 1e-6
            ####-------------------------------------------------------------------------
            ### next look at values above 4 km
            ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
            um_lt1km = np.where(um_data[m]['height'][0,:]>=4e3)
            twc_thresh_um[um_lt1km] = 1e-7

            ### find heights not yet assigned
            um_intZs = np.where(twc_thresh_um == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(um_data)):
                twc_thresh_um[um_intZs] = f(um_data[m]['height'][0,um_intZs].data)

            for t in range(0,np.size(um_data[m]['model_twc'],0)):
                for k in range(0,np.size(um_data[m]['model_twc'],1)):
                    if obs_data[twcvar][t,k] < twc_thresh_um[k]:
                        obs_data[twcvar][t,k] = np.nan
                        obs_data['iwc'][t,k] = np.nan
                    for m in range(0,len(um_data)):
                        if  um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                            um_data[m]['model_twc'][t,k] = np.nan
                            um_data[m]['model_iwc_filtered'][t,k] = np.nan
        if pmonc==True:
            m=0 # use first um model run for height grid definition
            twc_thresh_monc = np.zeros([np.size(monc_data[m]['model_twc'],1)])
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            monc_lt1km = np.where(monc_data[m][twc_zvar[m]][:]<=1e3)
            twc_thresh_monc[monc_lt1km] = 1e-6
            monc_lt1km = np.where(monc_data[m][twc_zvar[m]][:]>=4e3)
            twc_thresh_monc[monc_lt1km] = 1e-7
            monc_intZs = np.where(twc_thresh_monc == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(monc_data)):
                twc_thresh_monc[monc_intZs] = f(monc_data[m][twc_zvar[m]][monc_intZs].data)

            for t in range(0,np.size(monc_data[m]['model_twc'],0)):
                for k in range(0,np.size(monc_data[m]['model_twc'],1)):
                    for m in range(0,len(monc_data)):
                        if  monc_data[m]['model_twc'][t,k] < twc_thresh_monc[k]:
                            monc_data[m]['model_twc'][t,k] = np.nan
                            monc_data[m]['model_iwc'][t,k] = np.nan

            # ### plot profile of threshold as sanity check
            # plt.plot(twc_thresh_um, um_data[0]['height'][0,:])
            # plt.plot(twc_thresh_monc, monc_data[0]['z'][:]); plt.show()
            # print (twc_thresh_monc)

    ###----------------------------------------------------------------
    ###         Plot figure - Mean profiles
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
    plt.figure(figsize=(6,8.5))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)
    ### define axis instance
    ax1 = plt.gca()
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']

    #### ADIABATIC IWC (all times)
    plt.plot(np.nanmean(obs_data['iwc'],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid', zorder = obs_zorder)
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3,
        np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.xlim([0,0.2])
    plt.plot(np.nanmean(obs_data['iwc'],0)*1e3 - np.nanstd(obs_data['iwc'],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(obs_data['iwc'],0)*1e3 + np.nanstd(obs_data['iwc'],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    if pum==True:
        for m in range(0,len(um_data)):
            ax1.fill_betweenx(np.nanmean(um_data[m]['height'],0),np.nanmean(um_data[m]['model_iwc'],0)*1e3 - np.nanstd(um_data[m]['model_iwc']*1e3,0),
                np.nanmean(um_data[m]['model_iwc'],0)*1e3 + np.nanstd(um_data[m]['model_iwc'],0)*1e3, color = fcols[m], alpha = 0.05)
            plt.plot(np.nanmean(um_data[m]['model_iwc'],0)*1e3 - np.nanstd(um_data[m]['model_iwc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color =lcols[m], linewidth = 0.5)
            plt.plot(np.nanmean(um_data[m]['model_iwc'],0)*1e3 + np.nanstd(um_data[m]['model_iwc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color = lcols[m], linewidth = 0.5)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            ax1.fill_betweenx(monc_data[m][twc_zvar[m]],np.nanmean(monc_data[m]['model_iwc'],0)*1e3 - np.nanstd(monc_data[m]['model_iwc']*1e3,0),
                np.nanmean(monc_data[m]['model_iwc'],0)*1e3 + np.nanstd(monc_data[m]['model_iwc'],0)*1e3, color = fcolsmonc[m], alpha = 0.05)
            plt.plot(np.nanmean(monc_data[m]['model_iwc'],0)*1e3 - np.nanstd(monc_data[m]['model_iwc'],0)*1e3, monc_data[m][twc_zvar[m]],
                '--', color =lcolsmonc[m], linewidth = 0.5)
            plt.plot(np.nanmean(monc_data[m]['model_iwc'],0)*1e3 + np.nanstd(monc_data[m]['model_iwc'],0)*1e3, monc_data[m][twc_zvar[m]],
                '--', color = lcolsmonc[m], linewidth = 0.5)
    for m in range(0,len(um_data)):
        plt.plot(np.nanmean(um_data[m]['model_iwc'],0)*1e3,np.nanmean(um_data[m]['height'],0), color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.plot(np.nanmean(monc_data[m]['model_iwc'],0)*1e3,monc_data[m][twc_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)

    plt.xlabel('Ice water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    # plt.ylim([0,5e3])
    # plt.yticks(np.arange(0,5.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax1.set_yticklabels(ytlabels)
    plt.xlim([0,0.05])
    plt.xticks(np.arange(0,0.051,0.015))
    #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.0075))
    plt.legend()
    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if  pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_IWC-MTThresh' + twcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWC-MTThresh' + twcstr + '.png'
    else:
        if  pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_IWC' + twcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWC' + twcstr + '.png'


    plt.savefig(fileout)

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')


def plot_twcProfiles( obs_data,twcvar,twcstr, thresholding, plots_out_dir,dates, **args): #, lon, lat):

    obs_zorder = 1

    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                obs_zorder += len(monc_data)
                pmonc=True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                obs_zorder += len(um_data)
                pum=True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting TWC mean profiles:' + twcstr)
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc_filtered']
    if  pmonc==True:
        for m in range(0,len(monc_data)):
            # monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            # monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            # monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            # monc_data[m]['twc_height'] = monc_data[m][monc_data[m]['zvar']['ice_mmr_mean']]
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['twc_height'] = monc_data[m][monc_data[m]['zvar']['twc_tot_mean']]

    if thresholding == True:
        if pum==True:
            ####-------------------------------------------------------------------------
            ### from Michael's paper:
            ####    "we consider a grid point cloudy when cloud water exceeded
            ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
            ####    interpolation in between."
            ####-------------------------------------------------------------------------
            m=0 # use first um model run for height grid definition
            twc_thresh_um = np.zeros([np.size(um_data[m]['model_twc'],1)])
            ####-------------------------------------------------------------------------
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            um_lt1km = np.where(um_data[m]['height'][0,:]<=1e3)
            twc_thresh_um[um_lt1km] = 1e-6
            ####-------------------------------------------------------------------------
            ### next look at values above 4 km
            ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
            um_lt1km = np.where(um_data[m]['height'][0,:]>=4e3)
            twc_thresh_um[um_lt1km] = 1e-7

            ### find heights not yet assigned
            um_intZs = np.where(twc_thresh_um == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(um_data)):
                twc_thresh_um[um_intZs] = f(um_data[m]['height'][0,um_intZs].data)

            for t in range(0,np.size(um_data[m]['model_twc'],0)):
                for k in range(0,np.size(um_data[m]['model_twc'],1)):
                    if obs_data['twc'][t,k] < twc_thresh_um[k]:
                        obs_data['twc'][t,k] = np.nan
                        obs_data['lwc'][t,k] = np.nan
                    if obs_data['twc_ad'][t,k] < twc_thresh_um[k]:
                        obs_data['twc_ad'][t,k] = np.nan
                        obs_data['lwc_adiabatic'][t,k] = np.nan
                    if obs_data['twc_ad_nolwp'][t,k] < twc_thresh_um[k]:
                        obs_data['twc_ad_nolwp'][t,k] = np.nan
                        obs_data['lwc_adiabatic_inc_nolwp'][t,k] = np.nan
                    for m in range(0,len(um_data)):
                        if  um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                            um_data[m]['model_twc'][t,k] = np.nan
                            um_data[m]['model_iwc_filtered'][t,k] = np.nan
        if pmonc==True:
            m=0 # use first um model run for height grid definition
            twc_thresh_monc = np.zeros([np.size(monc_data[m]['model_twc'],1)])
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            monc_lt1km = np.where(monc_data[m]['twc_height'][:]<=1e3)
            twc_thresh_monc[monc_lt1km] = 1e-6
            monc_lt1km = np.where(monc_data[m]['twc_height'][:]>=4e3)
            twc_thresh_monc[monc_lt1km] = 1e-7
            monc_intZs = np.where(twc_thresh_monc == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            for m in range(0,len(monc_data)):
                twc_thresh_monc[monc_intZs] = f(monc_data[m]['twc_height'][monc_intZs].data)

            for t in range(0,np.size(monc_data[m]['model_twc'],0)):
                for k in range(0,np.size(monc_data[m]['model_twc'],1)):
                    for m in range(0,len(monc_data)):
                        if  monc_data[m]['model_twc'][t,k] < twc_thresh_monc[k]:
                            monc_data[m]['model_twc'][t,k] = np.nan
                            monc_data[m]['model_iwc'][t,k] = np.nan

            # ### plot profile of threshold as sanity check
            # plt.plot(twc_thresh_um, um_data[0]['height'][0,:])
            # plt.plot(twc_thresh_monc, monc_data[0][zvar][:]); plt.show()
            # print (twc_thresh_monc)

    ###----------------------------------------------------------------
    ###         Plot figure - Mean profiles
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
    plt.figure(figsize=(6,8.5))
    plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
            hspace = 0.4, wspace = 0.1)
    ### define axis instance
    ax1 = plt.gca()
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    #### ADIABATIC LWC (where there are HATPRO LWP data available)
    plt.plot(np.nanmean(obs_data[twcvar],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid' + twcstr, zorder = obs_zorder)
    ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[twcvar],0)*1e3 - np.nanstd(obs_data[twcvar],0)*1e3,
        np.nanmean(obs_data[twcvar],0)*1e3 + np.nanstd(obs_data[twcvar],0)*1e3, color = 'lightgrey', alpha = 0.5)
    # plt.xlim([0,0.2])
    plt.plot(np.nanmean(obs_data[twcvar],0)*1e3 - np.nanstd(obs_data[twcvar],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    plt.plot(np.nanmean(obs_data[twcvar],0)*1e3 + np.nanstd(obs_data[twcvar],0)*1e3, np.nanmean(obs_data['height'],0),
        '--', color = 'k', linewidth = 0.5)
    if pum==True:
        for m in range(0,len(um_data)):
            ax1.fill_betweenx(np.nanmean(um_data[m]['height'],0),np.nanmean(um_data[m]['model_twc'],0)*1e3 - np.nanstd(um_data[m]['model_twc']*1e3,0),
                np.nanmean(um_data[m]['model_twc'],0)*1e3 + np.nanstd(um_data[m]['model_twc'],0)*1e3, color = fcols[m], alpha = 0.05)
            plt.plot(np.nanmean(um_data[m]['model_twc'],0)*1e3 - np.nanstd(um_data[m]['model_twc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color =lcols[m], linewidth = 0.5)
            plt.plot(np.nanmean(um_data[m]['model_twc'],0)*1e3 + np.nanstd(um_data[m]['model_twc'],0)*1e3, np.nanmean(um_data[m]['height'],0),
                '--', color = lcols[m], linewidth = 0.5)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            ax1.fill_betweenx(monc_data[m]['twc_height'],np.nanmean(monc_data[m]['model_twc'],0)*1e3 - np.nanstd(monc_data[m]['model_twc']*1e3,0),
                np.nanmean(monc_data[m]['model_twc'],0)*1e3 + np.nanstd(monc_data[m]['model_twc'],0)*1e3, color = fcolsmonc[m], alpha = 0.05)
            plt.plot(np.nanmean(monc_data[m]['model_twc'],0)*1e3 - np.nanstd(monc_data[m]['model_twc'],0)*1e3, monc_data[m]['twc_height'],
                '--', color =lcolsmonc[m], linewidth = 0.5)
            plt.plot(np.nanmean(monc_data[m]['model_twc'],0)*1e3 + np.nanstd(monc_data[m]['model_twc'],0)*1e3, monc_data[m]['twc_height'],
                '--', color = lcolsmonc[m], linewidth = 0.5)
    if pum==True:
        for m in range(0,len(um_data)):
            plt.plot(np.nanmean(um_data[m]['model_twc'],0)*1e3,np.nanmean(um_data[m]['height'],0), color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.plot(np.nanmean(monc_data[m]['model_twc'],0)*1e3,monc_data[m]['twc_height'], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)

    plt.xlabel('Total water content [g m$^{-3}$]')
    plt.ylabel('Z [km]')
    # plt.ylim([0,5e3])
    # plt.yticks(np.arange(0,5.01e3,0.5e3))
    # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
    plt.ylim(ylims)
    plt.yticks(yticks)
#    plt.yticks(np.arange(0,6.01e3,0.5e3))
    ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax1.set_yticklabels(ytlabels)
    plt.xlim([0,0.45])
    plt.xticks(np.arange(0,0.41,0.05))
    #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))
    plt.legend()

    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_TWC-MTThresh' + twcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_TWC-MTThresh' +twcstr +'.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_TWC' + twcstr + '.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_TWC' +twcstr +'.png'

    plt.savefig(fileout)
    print(twcstr)
    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')


def removeSpinUp(monc_data,monc_spin):
    print('')
    print('remove MONC spinup time')
    for m in range(0,len(monc_data)):
        monc_var_list = list(monc_data[m].keys())
        monc_var_list.remove('time1')
        monc_var_list.remove('time2')
        if 'time3' in monc_data:
            monc_var_list.remove('time3')
        monc_var_list.remove('tvar')
        monc_var_list.remove('zvar')

        id1 = np.squeeze(np.argwhere(monc_data[m]['time1']<=monc_spin)) #1D data
        id2 = np.squeeze(np.argwhere(monc_data[m]['time2']<=monc_spin))
        if 'time3' in monc_data:
            id3 = np.squeeze(np.argwhere(monc_data[m]['time3']<=monc_spin))

        for j in range(0,len(monc_var_list)):
            if monc_data[m]['tvar'][monc_var_list[j]]=='time1':
            #if any(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time1'])):
                monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id1,0)
            elif monc_data[m]['tvar'][monc_var_list[j]]=='time2':
            #elif any(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time2'])):
                tmp2=np.argwhere(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time2']))
                if tmp2 == 0:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,0)
                elif tmp2 == 1:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,1)
                elif tmp2 == 2:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,2)
                elif tmp2 == 3:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id2,3)
            elif monc_data[m]['tvar'][monc_var_list[j]]=='time3':
            #elif any(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time2'])):
                tmp2=np.argwhere(np.array(monc_data[m][monc_var_list[j]].shape) == len(monc_data[m]['time3']))
                if tmp2 == 0:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id3,0)
                elif tmp2 == 1:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id3,1)
                elif tmp2 == 2:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id3,2)
                elif tmp2 == 3:
                    monc_data[m][monc_var_list[j]]=np.delete(monc_data[m][monc_var_list[j]],id3,3)
        monc_data[m]['time1']=np.delete(monc_data[m]['time1'],id1,0)
        monc_data[m]['time2']=np.delete(monc_data[m]['time2'],id2,0)
        if 'time3' in monc_data:
            monc_data[m]['time3']=np.delete(monc_data[m]['time3'],id3,0)
        monc_data[m]['time1']=monc_data[m]['time1']-monc_data[m]['time1'][0]
        monc_data[m]['time2']=monc_data[m]['time2']-monc_data[m]['time2'][0]
        if 'time3' in monc_data:
            monc_data[m]['time3']=monc_data[m]['time3']-monc_data[m]['time3'][0]

    return monc_data
    print('done')
    print('*************')
    print ('')



################################################################################
################################################################################
def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    names = ['20180913_oden_']

    sdate = dtime.datetime.strptime('2018091300','%Y%m%d%H')
    edate = dtime.datetime.strptime('2018091314','%Y%m%d%H')

    dates = [date2datenum(sdate),date2datenum(edate)]
    moccha_missing_files = ['20180813_oden_','20180910_oden_']   ### cloud radar not working    #,'20180914_oden_'

    #---- MONC SPIN UP TIME
    monc_spin = 6 *60 *60

    #---- SPLIT PROFILES IN TIME JUNKS
    prof_times=[[dates[0], dates[0]+4/24],
                [dates[0]+4/24, dates[0]+8/24],
                [dates[0]+8/24, dates[0]+14/24]]

    #which server
    machine = 'JASMIN'
    ### Choose observations vertical gridding used in Cloudnet processing (UM/IFS/RADAR)
    obs_switch = 'UM' #RADAR

    if machine =='JASMIN':
        plots_out_dir='/gws/nopw/j04/ncas_radar_vol1/jutta/PLOTS/CaseStudy/'
        if not os.path.exists(plots_out_dir):
            os.makedirs(plots_out_dir)

        um_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/UM/DATA/'
        cn_obs_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/Cloudnet/Observations/measurements_V7/'
        cn_um_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/Cloudnet/UM/'
        monc_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
        #monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/MONC/output/'
        monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'

    elif machine =='LEEDS':
        ### Set output directory for plots
        plots_out_dir='/nfs/a96/MOCCHA/working/jutta/plots/CaseStudies/ModelComparison/'
        if not os.path.exists(plots_out_dir):
            os.makedirs(plots_out_dir)

        um_root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        cn_obs_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/Halo/measurements_V7/'
        cn_um_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/output/metum/V7/'
        monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'

    ### CHOOSE MONC RUNS
    # m_out_dir = ['5_control_20180913T0000Z_Wsub-1.5_Fletcher/',
    #              '6_control_20180913T0000Z_Wsub-1.5-1km/',
    #              '8_control_20180913T0000Z_Wsub-1.0-1km/',
    #              '9_control_20180913T0000Z_Wsub-0.5-1km/']
    m_out_dir = [#'22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/']
               '23_20180913T0000Z_6hSpin-up_12h0600-0000thTend_20h1200-0600thTend/']
            #'4_control_20180913T0000Z_Wsub-1.5/',
    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    mlabel=[]
    moutstr=[]
    for m in range(0, len(m_out_dir)):
        if m_out_dir[m][:1] == '3':
            mlabel.append('MONC nosub')
            moutstr.append('Mnowsub')
        elif m_out_dir[m][:1] == '4':
            mlabel.append('MONC Wsub1.5')
            moutstr.append('Mwsub')
        elif m_out_dir[m][:1] == '5':
            mlabel.append('MONC Wsub1.5 \n Fletcher')
            moutstr.append('Mwsubfle')
        elif m_out_dir[m][:1] == '6':
            mlabel.append('MONC Wsub1.5-1km')
            moutstr.append('Mwsub1.5-1km')
        elif m_out_dir[m][:1] == '7':
            mlabel.append('MONC Wsub1.5-1km \n solACC-100')
            moutstr.append('Mwsub1kmsolACC100')
        elif m_out_dir[m][:1] == '8':
            mlabel.append('MONC Wsub1.0-1km')
            moutstr.append('Mwsub1.0-1km')
        elif m_out_dir[m][:1] == '9':
            mlabel.append('MONC Wsub0.5-1km')
            moutstr.append('Mwsub0.5-1km')
        elif m_out_dir[m][:2] == '20':
            mlabel.append('MONC qinit2 800m \n thqvTend noice')
            moutstr.append('qin2-thqvTend-noice')
        elif m_out_dir[m][:2] == '22':
            mlabel.append('MONC qinit2 800m \n thForcing-0000-0600')
            moutstr.append('MONC-22')
        elif m_out_dir[m][:2] == '23':
            mlabel.append('MONC thForcing-0600-0000')
            moutstr.append('MONC-23')
        else:
            label.append('undefined_label')
            moutstr.append('')

    ### -----------------------------------------------------------------
    ### create monc filenames
    monc_filename=[]
    for m in range(0, len(m_out_dir)):
        fname=glob.glob(monc_root_dir + m_out_dir[m] +'*.nc')
        monc_filename.append(fname)
    monc_3d_filename=[]
    for m in range(0, len(m_out_dir)):
        fname=glob.glob(monc_avg_dir + m_out_dir[m] +'3d*npy')
        monc_3d_filename.append(fname)


    #################################################################
    ###     READ IN MONC DATA
    #################################################################
    print ('Loading MONC data:')
    print ('')
    ###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)
    ### time1 = 'time_series_30_600',time2='time_series_30_60'
    monc_var_list =[['z', 'zn','LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
                    ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction'],
                    ['liquid_mmr_mean','ice_mmr_mean','graupel_mmr_mean','snow_mmr_mean']]
                #    ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp'],
                #    ['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]


    monc_var_3d_list =['T_mean','p_mean','th_mean','rho_mean','q_vapour_mean','q_cloud_liquid_mass_mean',
                        'q_ice_mass_mean','q_snow_mass_mean','q_graupel_mass_mean',
                        'twc_tot_mean','iwc_tot_mean','lwc_tot_mean']
                        #'z', 'zn', 'u_mean', 'v_mean', 'w_mean', 'q_vapour_mean',
                        # 'time1', 'time2', 'p_mean', 'T_mean', 'th_mean', 'rho_mean',
                        # 'q_cloud_liquid_mass_mean', 'q_ice_mass_mean', 'q_snow_mass_mean',
                        # 'q_graupel_mass_mean', 'iwc_tot_mean', 'lwc_tot_mean', 'twc_tot_mean', 'zvar', 'tvar']
    ncm = {}
    monc_data = {}
    for m in range(0, len(m_out_dir)):
        for n in range(0, len(monc_filename[m])):
            print(monc_filename[m][n])
            ncm = Dataset(monc_filename[m][n],'r')
            if n == 0:
                monc_data[m]={}
                zvar={}
                tvar={}
            full_var_list=[]
            time_var_list=[]
            for var in ncm.variables:
                if 'time' in str(var):
                    print(var)
                    time_var_list=time_var_list+[var]
            full_var_list=monc_var_list.copy()
            full_var_list[0]=time_var_list+monc_var_list[0]
            for c in range(0,len(full_var_list)):
                for j in range(0,len(full_var_list[c])):
                    var = full_var_list[c][j]
                    if n == 0:
                        monc_data[m][var] = ncm.variables[var][:]
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
                        if time_var_list[0] in str(tmp):
                            tvar[var]='time1'        #time_series_30_600
                        elif time_var_list[1] in str(tmp):
                            tvar[var]='time2'        #time_series_30_60
                        if len(time_var_list)>2:
                            if time_var_list[2] in str(tmp):
                                tvar[var]='time3'
                    else:
                        if var =='z': continue
                        elif var =='zn': continue
                        else:
                            print('appending ' + var)
                            monc_data[m][var]=np.append(monc_data[m][var],ncm.variables[var][:],axis=0)
    #loading 3d variables
        for n in range(0,len(monc_3d_filename[m])):
            print(monc_3d_filename[m][n])
            #afile = open(monc_3d_filename[m][0], "rb")
            pyd = np.load(monc_3d_filename[m][n],allow_pickle=True).item()   #pickle.load(afile)
            #pyd['zvar']['q_vapour_mean']=['zn']
            #pyd['tvar']['q_vapour_mean']=['time1']
            for c in range(0,len(monc_var_3d_list)):
                var = monc_var_3d_list[c]
                if n == 0:
                    zvar[var]=pyd['zvar'][var]
                    tvar[var]=pyd['tvar'][var]
                    monc_data[m][var] = pyd[var]
                else:
                    monc_data[m][var] =np.append(monc_data[m][var],pyd[var],axis=0)

        # for v1 in monc_data[m]:
        #     print(v1 + ':')
        #     print(monc_data[m][v1].shape)
        #

        monc_data[m]['zvar']=zvar
        monc_data[m]['tvar']=tvar
        monc_data[m]['time1']=monc_data[m][time_var_list[0]] #1d data
        monc_data[m]['time2']=monc_data[m][time_var_list[1]] #2d data
        monc_data[m].pop(time_var_list[0])
        monc_data[m].pop(time_var_list[1])
        if len(time_var_list)>2:
            monc_data[m]['time3']=monc_dat[m][time_var_list[2]] #2d data
            monc_data[m].pop(time_var_list[2])

    print (' Monc data Loaded!')
    ##################################################################################################################################
    ## -------------------------------------------------------------
    ## remove spin up time from monc data
    ## -------------------------------------------------------------
    # monc_data=removeSpinUp(monc_data,monc_spin)

    ## -------------------------------------------------------------
    ## convert monc time to datenum
    ## -------------------------------------------------------------
    for m in range(0,len(monc_data)):
        monc_data[m]['time2']=dates[0] + monc_data[m]['time2']/60/60/24
        monc_data[m]['time1']=dates[0] + monc_data[m]['time1']/60/60/24
        if 'time3' in monc_data:
            monc_data[m]['time3']=dates[0] + monc_data[m]['time3']/60/60/24

###################################################################################################################
###################################################################################################################
################################################ FIGURES ##########################################################
###################################################################################################################
###################################################################################################################
    ## DEFINE which Obs variable should be plotted
    ### options: lwc                     = measured lwc
    ###          lwc_adiabatic           = adiabatic lwc calculated
    ###          lwc_adiabatic_inc_nolwp = adiabatic lwc calculated incl. times when no HATPRO LWP available
    ### options: twc                     = using lwc
    ###          twc_ad                  = using lwc_adiabatic
    ###          lwc_ad_nolwp            = using lwc_adiabatic_inc_nolwp

    # lwcvar='lwc_adiabatic_inc_nolwp'
    # lwcstr='-adincNoLWP'   #''/'-adiabatic' / '-adincNoLWP'
    # twcvar='twc_ad_nolwp'
    # twcstr='-adincNoLWP'  #''/'-adiabatic' / '-adincNoLWP'
    # #
    lwcvar='lwc'
    lwcstr=''   #''/'-adiabatic' / '-adincNoLWP'
    twcvar='twc'
    twcstr=''  #''/'-adiabatic' / '-adincNoLWP'
    #
    ### select if Michaels cloud mask should be used for plotting profiles
    thresholding =True
        # -------------------------------------------------------------
    # Cloudnet plot: Plot Cv statistics from drift period
    # -------------------------------------------------------------
    #figure = plot_CvProfiles(um_data, ifs_data, misc_data, ra2t_data, obs_data, month_flag, missing_files, cn_um_out_dir, doy, obs, obs_switch)
    figure = plot_lwcProfiles(obs_data, lwcvar,lwcstr,thresholding, plots_out_dir,dates,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_iwcProfiles(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates, um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_twcProfiles(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_lwcProfiles_split(obs_data, lwcvar,lwcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_iwcProfiles_split(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)


    # -------------------------------------------------------------
    # Cloudnet plot: Plot contour timeseries
    # -------------------------------------------------------------
    # figure = plot_CvTimeseries(obs_data,plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_LWCTimeseries(obs_data,  lwcvar,lwcstr, plots_out_dir, dates, um_data=um_data,label=label, outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_TWCTimeseries( obs_data, twcvar,twcstr, plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr,monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_IWCTimeseries(obs_data, plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr,monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_monc_comparison(obs_data,  lwcvar,lwcstr, plots_out_dir, dates, um_data=um_data,label=label, outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)


    # -------------------------------------------------------------
    # plot LWP timeseries
    # -------------------------------------------------------------
    #figure = plot_lwp(obs_data, plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr,monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)

    # -------------------------------------------------------------
    # FIN.
    # -------------------------------------------------------------
    END_TIME = time.time()
    print ('******')
    print ('')
    print ('End: ' + time.strftime("%c"))
    print ('')

if __name__ == '__main__':

    main()
