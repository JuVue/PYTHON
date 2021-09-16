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


def calc_TWC(obs_data,thresholding, **args):

    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                pmonc =True
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                    pum =True

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    if pmonc==True:
        lwc_zvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            monc_data[m]['model_iwc'] = monc_data[m]['iwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]

    if thresholding == True:
            ####-------------------------------------------------------------------------
            ### from Michael's paper:
            ####    "we consider a grid point cloudy when cloud water exceeded
            ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
            ####    interpolation in between."
            ####-------------------------------------------------------------------------
            twc_thresh_um = np.zeros([np.size(obs_data['twc'],1)])
            ####-------------------------------------------------------------------------
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            um_lt1km = np.where(obs_data['height'][0,:]<=1e3)
            twc_thresh_um[um_lt1km] = 1e-6
            ####-------------------------------------------------------------------------
            ### next look at values above 4 km
            ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
            um_lt1km = np.where(obs_data[m]['height'][0,:]>=4e3)
            twc_thresh_um[um_lt1km] = 1e-7

            ### find heights not yet assigned
            um_intZs = np.where(twc_thresh_um == 0.0)

            ### interpolate for twc_thresh_um
            x = [1e-6, 1e-7]
            y = [1e3, 4e3]
            f = interp1d(y, x)
            twc_thresh_um[um_intZs] = f(obs_data['height'][0,um_intZs].data)

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
                    if pum==True:
                        for m in range(0,len(um_data)):
                            if  um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                                um_data[m]['model_twc'][t,k] = np.nan
                                um_data[m]['model_lwc'][t,k] = np.nan
        if pmonc==True:
            m=0 # use first monc model run for height grid definition
            twc_thresh_monc = np.zeros([np.size(monc_data[m]['model_twc'],1)])
            ### first look at values below 1 km
            ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
            monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]<=1e3)
            twc_thresh_monc[monc_lt1km] = 1e-6
            monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]>=4e3)
            twc_thresh_monc[monc_lt1km] = 1e-7
            monc_intZs = np.where(twc_thresh_monc == 0.0)

            ### interpolate for twc_thresh_monc
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

    if ((pum==True) and (pmonc==True)):
        return obs_data,um_data,monc_data
    elif (pum==True)):
        return obs_data,um_data
    elif (pmonc==True)):
        return obs_data,monc_data
    else:
        return obs_data



def plot_CvTimeseries(obs_data, plots_out_dir,dates,  **args):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
                pmonc = True
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
                pum = True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ylims=[0,3]
    yticks=np.arange(0,3e3,0.5e3)
    # ylims=[0,6]
    # yticks=np.arange(0,6e3,0.5e3)
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
    if pum ==True:
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

    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.subplot(numsp,1,numsp-len(monc_data)+1+m)
            ax = plt.gca()
            plt.contourf(monc_data[m][monc_data[m]['tvar']['total_cloud_fraction']],
                np.squeeze(monc_data[m][monc_data[m]['zvar']['total_cloud_fraction']][:]),
                np.transpose(monc_data[m]['total_cloud_fraction']),
                np.arange(0,1.1,0.1),
                cmap = newcmp,
                zorder = 1)
            # plt.plot(np.squeeze(obs['inversions']['doy']),np.squeeze(obs['inversions']['invbase']), 'k', linewidth = 1.0)
            # plt.plot(data1['time_hrly'][::6], bldepth1[::6], 'k', linewidth = 1.0)
            plt.ylabel('Z [km]')
            plt.ylim(ylims)
            plt.yticks(yticks)
            ax.set_yticklabels(ytlabels)
            ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
            #plt.xlim(monc_data['time2']/60/60])
            # xticks=np.arange(np.floor(monc_data[m]['time2'][0]/60/60),np.ceil(monc_data[m]['time2'][-1]/60/60)+1,2,dtype=int)
            # xlabs=[x*100 for x in xticks]
            # xlabs=["{0:-04d}".format(t) for t in xlabs]
            # plt.xticks(xticks)
            # ax.set_xticklabels(xlabs)
            plt.xlabel('Time (UTC)')
            nans = ax.get_ylim()
            ax2 = ax.twinx()
            ax2.set_ylabel(mlabel[m], rotation = 270, labelpad = 27)
            ax2.set_yticks([])

    dstr=datenum2date(dates[1])
    if pmonc==True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + 'CvTimeseries.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + 'CvTimeseries.png'
    print(fileout)
    plt.savefig(fileout)

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

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
            um_data[m]['model_iwc'][um_data[m]['model_iwc'] <= 0.0] = np.nan

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
            plt.contourf(um_data[m]['time'], np.squeeze(um_data[m]['height'][0,:]), np.transpose(um_data[m]['model_iwc'])*1e3,
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
    print ('Plotting TWC timeseries for case study:')
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
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']

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


    obs_data,um_data,monc_data=calc_TWC(obs_data, thresholding um_data=um_data,monc_data=monc_data)
    embed()
    # ###----------------------------------------------------------------
    # ###         Calculate total water content
    # ###----------------------------------------------------------------
    # obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    # obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    # obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    # if pum==True:
    #     for m in range(0,len(um_data)):
    #         um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    # if pmonc==True:
    #     lwc_zvar=[]
    #     for m in range(0,len(monc_data)):
    #         #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
    #         #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
    #         #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
    #         monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
    #         monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
    #         #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
    #         lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]
    #
    # if thresholding == True:
    #     if pum==True:
    #
    #         ####-------------------------------------------------------------------------
    #         ### from Michael's paper:
    #         ####    "we consider a grid point cloudy when cloud water exceeded
    #         ####    0.001 (0.0001) g kg-1 below 1 km (above 4 km), with linear
    #         ####    interpolation in between."
    #         ####-------------------------------------------------------------------------
    #         m=0 # use first um model run for height grid definition
    #         twc_thresh_um = np.zeros([np.size(um_data[m]['model_twc'],1)])
    #         ####-------------------------------------------------------------------------
    #         ### first look at values below 1 km
    #         ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    #         um_lt1km = np.where(um_data[m]['height'][0,:]<=1e3)
    #         twc_thresh_um[um_lt1km] = 1e-6
    #         ####-------------------------------------------------------------------------
    #         ### next look at values above 4 km
    #         ###     find Z indices >= 4km, then set twc_thresh values to 1e-7
    #         um_lt1km = np.where(um_data[m]['height'][0,:]>=4e3)
    #         twc_thresh_um[um_lt1km] = 1e-7
    #
    #         ### find heights not yet assigned
    #         um_intZs = np.where(twc_thresh_um == 0.0)
    #
    #         ### interpolate for twc_thresh_um
    #         x = [1e-6, 1e-7]
    #         y = [1e3, 4e3]
    #         f = interp1d(y, x)
    #         for m in range(0,len(um_data)):
    #             twc_thresh_um[um_intZs] = f(um_data[m]['height'][0,um_intZs].data)
    #
    #         for t in range(0,np.size(um_data[m]['model_twc'],0)):
    #             for k in range(0,np.size(um_data[m]['model_twc'],1)):
    #                 if obs_data['twc'][t,k] < twc_thresh_um[k]:
    #                     obs_data['twc'][t,k] = np.nan
    #                     obs_data['lwc'][t,k] = np.nan
    #                 if obs_data['twc_ad'][t,k] < twc_thresh_um[k]:
    #                     obs_data['twc_ad'][t,k] = np.nan
    #                     obs_data['lwc_adiabatic'][t,k] = np.nan
    #                 if obs_data['twc_ad_nolwp'][t,k] < twc_thresh_um[k]:
    #                     obs_data['twc_ad_nolwp'][t,k] = np.nan
    #                     obs_data['lwc_adiabatic_inc_nolwp'][t,k] = np.nan
    #                 for m in range(0,len(um_data)):
    #                     if  um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
    #                         um_data[m]['model_twc'][t,k] = np.nan
    #                         um_data[m]['model_lwc'][t,k] = np.nan
    #     if pmonc==True:
    #         m=0 # use first um model run for height grid definition
    #         twc_thresh_monc = np.zeros([np.size(monc_data[m]['model_twc'],1)])
    #         ### first look at values below 1 km
    #         ###     find Z indices <= 1km, then set twc_thresh values to 1e-6
    #         monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]<=1e3)
    #         twc_thresh_monc[monc_lt1km] = 1e-6
    #         monc_lt1km = np.where(monc_data[m][lwc_zvar[m]][:]>=4e3)
    #         twc_thresh_monc[monc_lt1km] = 1e-7
    #         monc_intZs = np.where(twc_thresh_monc == 0.0)
    #
    #         ### interpolate for twc_thresh_um
    #         x = [1e-6, 1e-7]
    #         y = [1e3, 4e3]
    #         f = interp1d(y, x)
    #         for m in range(0,len(monc_data)):
    #             twc_thresh_monc[monc_intZs] = f(monc_data[m][lwc_zvar[m]][monc_intZs].data)
    #
    #         for t in range(0,np.size(monc_data[m]['model_twc'],0)):
    #             for k in range(0,np.size(monc_data[m]['model_twc'],1)):
    #                 for m in range(0,len(monc_data)):
    #                     if  monc_data[m]['model_twc'][t,k] < twc_thresh_monc[k]:
    #                         monc_data[m]['model_twc'][t,k] = np.nan
    #                         monc_data[m]['model_lwc'][t,k] = np.nan
    #
    #         # ### plot profile of threshold as sanity check
    #         # plt.plot(twc_thresh_um, um_data[0]['height'][0,:])
    #         # plt.plot(twc_thresh_monc, monc_data[0]['z'][:]); plt.show()
    #         # print (twc_thresh_monc)


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
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
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
                            um_data[m]['model_iwc'][t,k] = np.nan
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
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
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
                            um_data[m]['model_iwc'][t,k] = np.nan
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

def plot_NdropProfiles_split(obs_data, lwcvar, lwcstr, thresholding, plots_out_dir, dates, prof_time, **args): #, lon, lat):

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
            elif  list(args.keys())[n] == 'raw_data':
                raw_data=args[list(args.keys())[n]]
                obs_zorder += len(raw_data)
                praw = True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,2]
    yticks=np.arange(0,2e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting Ndrop mean profiles split times:')
    print ('')
    # print (raw_data[0].keys())

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    if pmonc==True:
        lwc_zvar=[]
        lwc_tvar=[]
        ndrop_zvar=[]
        ndrop_tvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]
            lwc_tvar+=[monc_data[m]['tvar']['lwc_tot_mean']]
            ndrop_zvar+=[monc_data[m]['zvar']['ndrop_tot_mean']]
            ndrop_tvar+=[monc_data[m]['tvar']['ndrop_tot_mean']]

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
                        if um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                            um_data[m]['model_twc'][t,k] = np.nan
                            um_data[m]['model_lwc'][t,k] = np.nan
                            if praw == True:
                                raw_data[m]['qliq'][t,k] = np.nan       ## apply nan mask to qliq based on LWC
                                raw_data[m]['qnliq'][t,k] = np.nan       ## apply nan mask to ndrop based on LWC
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
                            monc_data[m]['model_lwc'][t,k] = np.nan
                            monc_data[m]['ndrop_tot_mean'][t,k] = np.nan
                            monc_data[m]['nisg_tot_mean'][t,k] = np.nan

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
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    ### define axis instance
    ####LWC
    plt.figure(figsize=(14,7))
    plt.subplots_adjust(top = 0.8, bottom = 0.1, right = 0.98, left = 0.08)

    for pt in range(0,len(prof_time)):
        plt.subplot(1,len(prof_time),pt+1)
        ax1 = plt.gca()
        sstr=datenum2date(prof_time[pt][0])
        estr=datenum2date(prof_time[pt][1])
        plt.title(sstr.strftime("%H") +'-' + estr.strftime("%H") + ' UTC')
        # obsid= np.squeeze(np.argwhere((raw_data['time']>=prof_time[pt][0]) & (raw_data['time']<prof_time[pt][1])))
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid'  + lwcstr, zorder = obs_zorder)
        # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3,
        #     np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)
        if praw==True:
            for m in range(0,len(raw_data)):
                id=  np.squeeze(np.argwhere((raw_data[m]['time']>=prof_time[pt][0]) & (raw_data[m]['time']<prof_time[pt][1])))
                ax1.fill_betweenx(raw_data[m]['height'],np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 - np.nanstd(raw_data[m]['qnliq'][id,:]/1e6,0),
                    np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 + np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, color = fcols[m], alpha = 0.05)
                plt.plot(np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 - np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, raw_data[m]['height'],
                    '--', color =lcols[m], linewidth = 0.5)
                plt.plot(np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 + np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, raw_data[m]['height'],
                    '--', color = lcols[m], linewidth = 0.5)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][ndrop_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][ndrop_tvar[m]]<prof_time[pt][1])))
                ax1.fill_betweenx(monc_data[m][ndrop_zvar[m]],np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 - np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:]/1e6,0),
                    np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 + np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, color = fcolsmonc[m], alpha = 0.05)
                plt.plot(np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 - np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, monc_data[m][ndrop_zvar[m]],
                    '--', color =lcolsmonc[m], linewidth = 0.5)
                plt.plot(np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 + np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, monc_data[m][ndrop_zvar[m]],
                    '--', color = lcolsmonc[m], linewidth = 0.5)
        if praw==True:
            for m in range(0,len(raw_data)):
                id= np.squeeze(np.argwhere((raw_data[m]['time']>=prof_time[pt][0]) & (raw_data[m]['time']<prof_time[pt][1])))
                plt.plot(np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6,raw_data[m]['height'], color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][ndrop_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][ndrop_tvar[m]]<prof_time[pt][1])))
                plt.plot(np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6,monc_data[m][ndrop_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)
        if pt == 1:
            plt.legend(bbox_to_anchor=(1.5, 1.05), loc=4, ncol=2)


        plt.xlabel('N$_{drop}$ [cm$^{-3}$]')
        plt.ylabel('Z [km]')
        # plt.ylim([0,5e3])
        # plt.yticks(np.arange(0,5.01e3,0.5e3))
        # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.set_yticklabels(ytlabels)
        # plt.xlim([0,0.4])
        # plt.xticks(np.arange(0,0.45,0.1))
        # #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
        # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))


    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_Ndrop-MTThresh' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_Ndrop-MTThresh'+ lwcstr + '_split.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_Ndrop' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_Ndrop'+ lwcstr + '_split.png'

    plt.savefig(fileout, dpi=300)
    plt.close()
    # plt.show()
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_lwcProfiles_split(obs_data,lwcvar,lwcstr, thresholding, plots_out_dir,dates,prof_time, **args): #, lon, lat):

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

    ylims=[0,2]
    yticks=np.arange(0,2e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting LWC mean profiles split times:' + lwcvar)
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    if pmonc==True:
        lwc_zvar=[]
        lwc_tvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]
            lwc_tvar+=[monc_data[m]['tvar']['lwc_tot_mean']]

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
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    ### define axis instance
    ####LWC
    plt.figure(figsize=(14,7))
    plt.subplots_adjust(top = 0.8, bottom = 0.1, right = 0.98, left = 0.08)

    for pt in range(0,len(prof_time)):
        plt.subplot(1,len(prof_time),pt+1)
        ax1 = plt.gca()
        sstr=datenum2date(prof_time[pt][0])
        estr=datenum2date(prof_time[pt][1])
        plt.title(sstr.strftime("%H") +'-' + estr.strftime("%H") + ' UTC')
        obsid= np.squeeze(np.argwhere((obs_data['time']>=prof_time[pt][0]) & (obs_data['time']<prof_time[pt][1])))
        plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid'  + lwcstr, zorder = obs_zorder)
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3,
            np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.xlim([0,0.2])
        plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id=  np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                ax1.fill_betweenx(np.nanmean(um_data[m]['height'],0),np.nanmean(um_data[m]['model_lwc'][id,:],0)*1e3 - np.nanstd(um_data[m]['model_lwc'][id,:]*1e3,0),
                    np.nanmean(um_data[m]['model_lwc'][id,:],0)*1e3 + np.nanstd(um_data[m]['model_lwc'][id,:],0)*1e3, color = fcols[m], alpha = 0.05)
                plt.plot(np.nanmean(um_data[m]['model_lwc'][id,:],0)*1e3 - np.nanstd(um_data[m]['model_lwc'][id,:],0)*1e3, np.nanmean(um_data[m]['height'],0),
                    '--', color =lcols[m], linewidth = 0.5)
                plt.plot(np.nanmean(um_data[m]['model_lwc'][id,:],0)*1e3 + np.nanstd(um_data[m]['model_lwc'][id,:],0)*1e3, np.nanmean(um_data[m]['height'],0),
                    '--', color = lcols[m], linewidth = 0.5)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][lwc_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][lwc_tvar[m]]<prof_time[pt][1])))
                ax1.fill_betweenx(monc_data[m][lwc_zvar[m]],np.nanmean(monc_data[m]['model_lwc'][id,:],0)*1e3 - np.nanstd(monc_data[m]['model_lwc'][id,:]*1e3,0),
                    np.nanmean(monc_data[m]['model_lwc'][id,:],0)*1e3 + np.nanstd(monc_data[m]['model_lwc'][id,:],0)*1e3, color = fcolsmonc[m], alpha = 0.05)
                plt.plot(np.nanmean(monc_data[m]['model_lwc'][id,:],0)*1e3 - np.nanstd(monc_data[m]['model_lwc'][id,:],0)*1e3, monc_data[m][lwc_zvar[m]],
                    '--', color =lcolsmonc[m], linewidth = 0.5)
                plt.plot(np.nanmean(monc_data[m]['model_lwc'][id,:],0)*1e3 + np.nanstd(monc_data[m]['model_lwc'][id,:],0)*1e3, monc_data[m][lwc_zvar[m]],
                    '--', color = lcolsmonc[m], linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id= np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                plt.plot(np.nanmean(um_data[m]['model_lwc'][id,:],0)*1e3,np.nanmean(um_data[m]['height'],0), color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][lwc_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][lwc_tvar[m]]<prof_time[pt][1])))
                plt.plot(np.nanmean(monc_data[m]['model_lwc'][id,:],0)*1e3,monc_data[m][lwc_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)
        if pt == 1:
            plt.legend(bbox_to_anchor=(1.5, 1.05), loc=4, ncol=2)


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
        plt.xticks(np.arange(0,0.45,0.1))
        #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
        ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))


    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_LWC-MTThresh' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWC-MTThresh'+ lwcstr + '_split.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_LWC' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWC'+ lwcstr + '_split.png'

    plt.savefig(fileout, dpi=300)
    plt.close()
    # plt.show()
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_iwcProfiles_split(obs_data,twcvar,twcstr, thresholding, plots_out_dir,dates,prof_time, **args): #, lon, lat):

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

    ylims=[0,2]
    yticks=np.arange(0,2e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting IWC mean profiles split times:' + twcvar)
    print ('')

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    if pmonc==True:
        twc_zvar=[]
        twc_tvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            monc_data[m]['model_iwc'] = monc_data[m]['iwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            twc_zvar+=[monc_data[m]['zvar']['iwc_tot_mean']]
            twc_tvar+=[monc_data[m]['tvar']['iwc_tot_mean']]

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
                            um_data[m]['model_iwc'][t,k] = np.nan
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
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    ### define axis instance
    ####IWC
    plt.figure(figsize=(14,7))
    plt.subplots_adjust(top = 0.8, bottom = 0.1, right = 0.98, left = 0.08)

    for pt in range(0,len(prof_time)):
        plt.subplot(1,len(prof_time),pt+1)
        ax1 = plt.gca()
        sstr=datenum2date(prof_time[pt][0])
        estr=datenum2date(prof_time[pt][1])
        plt.title(sstr.strftime("%H") +'-' + estr.strftime("%H") + ' UTC')
        obsid= np.squeeze(np.argwhere((obs_data['time']>=prof_time[pt][0]) & (obs_data['time']<prof_time[pt][1])))
        plt.plot(np.nanmean(obs_data['iwc'][obsid,:],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid'  + twcstr, zorder = obs_zorder)
        ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data['iwc'][obsid,:],0)*1e3 - np.nanstd(obs_data['iwc'][obsid,:],0)*1e3,
            np.nanmean(obs_data['iwc'][obsid,:],0)*1e3 + np.nanstd(obs_data['iwc'][obsid,:],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.xlim([0,0.2])
        plt.plot(np.nanmean(obs_data['iwc'][obsid,:],0)*1e3 - np.nanstd(obs_data['iwc'][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs_data['iwc'][obsid,:],0)*1e3 + np.nanstd(obs_data['iwc'][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
            '--', color = 'k', linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id=  np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                ax1.fill_betweenx(np.nanmean(um_data[m]['height'],0),np.nanmean(um_data[m]['model_iwc'][id,:],0)*1e3 - np.nanstd(um_data[m]['model_iwc'][id,:]*1e3,0),
                    np.nanmean(um_data[m]['model_iwc'][id,:],0)*1e3 + np.nanstd(um_data[m]['model_iwc'][id,:],0)*1e3, color = fcols[m], alpha = 0.05)
                plt.plot(np.nanmean(um_data[m]['model_iwc'][id,:],0)*1e3 - np.nanstd(um_data[m]['model_iwc'][id,:],0)*1e3, np.nanmean(um_data[m]['height'],0),
                    '--', color =lcols[m], linewidth = 0.5)
                plt.plot(np.nanmean(um_data[m]['model_iwc'][id,:],0)*1e3 + np.nanstd(um_data[m]['model_iwc'][id,:],0)*1e3, np.nanmean(um_data[m]['height'],0),
                    '--', color = lcols[m], linewidth = 0.5)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][twc_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][twc_tvar[m]]<prof_time[pt][1])))
                ax1.fill_betweenx(monc_data[m][twc_zvar[m]],np.nanmean(monc_data[m]['model_iwc'][id,:],0)*1e3 - np.nanstd(monc_data[m]['model_iwc'][id,:]*1e3,0),
                    np.nanmean(monc_data[m]['model_iwc'][id,:],0)*1e3 + np.nanstd(monc_data[m]['model_iwc'][id,:],0)*1e3, color = fcolsmonc[m], alpha = 0.05)
                plt.plot(np.nanmean(monc_data[m]['model_iwc'][id,:],0)*1e3 - np.nanstd(monc_data[m]['model_iwc'][id,:],0)*1e3, monc_data[m][twc_zvar[m]],
                    '--', color =lcolsmonc[m], linewidth = 0.5)
                plt.plot(np.nanmean(monc_data[m]['model_iwc'][id,:],0)*1e3 + np.nanstd(monc_data[m]['model_iwc'][id,:],0)*1e3, monc_data[m][twc_zvar[m]],
                    '--', color = lcolsmonc[m], linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id= np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                plt.plot(np.nanmean(um_data[m]['model_iwc'][id,:],0)*1e3,np.nanmean(um_data[m]['height'],0), color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][twc_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][twc_tvar[m]]<prof_time[pt][1])))
                plt.plot(np.nanmean(monc_data[m]['model_iwc'][id,:],0)*1e3,monc_data[m][twc_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)
        if pt == 1:
            plt.legend(bbox_to_anchor=(1.5, 1.05), loc=4, ncol=2)


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


    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_IWC-MTThresh' + twcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWC-MTThresh'+ twcstr + '_split.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_IWC' + twcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_IWC'+ twcstr + '_split.png'

    plt.savefig(fileout, dpi=300)
    plt.close()
    # plt.show()
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_NisgProfiles_split(obs_data, lwcvar, lwcstr, thresholding, plots_out_dir, dates, prof_time, **args): #, lon, lat):

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
            elif  list(args.keys())[n] == 'raw_data':
                raw_data=args[list(args.keys())[n]]
                obs_zorder += len(raw_data)
                praw = True
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    ylims=[0,2]
    yticks=np.arange(0,2e3,0.5e3)
    ytlabels=yticks/1e3


    print ('******')
    print ('')
    print ('Plotting Nice (UM) and Nisg (MONC) mean profiles split times:')
    print ('')
    # print (raw_data[0].keys())

    ###----------------------------------------------------------------
    ###         Calculate total water content
    ###----------------------------------------------------------------
    obs_data['twc'] = obs_data['lwc'] + obs_data['iwc']
    obs_data['twc_ad'] = obs_data['lwc_adiabatic'] + obs_data['iwc']
    obs_data['twc_ad_nolwp'] = obs_data['lwc_adiabatic_inc_nolwp'] + obs_data['iwc']
    if pum==True:
        for m in range(0,len(um_data)):
            um_data[m]['model_twc'] = um_data[m]['model_lwc'] + um_data[m]['model_iwc']
    if pmonc==True:
        lwc_zvar=[]
        lwc_tvar=[]
        nisg_zvar=[]
        nisg_tvar=[]
        for m in range(0,len(monc_data)):
            #monc_data[m]['model_iwc']= (monc_data[m]['ice_mmr_mean']+monc_data[m]['graupel_mmr_mean']+monc_data[m]['snow_mmr_mean'])*monc_data[m]['rho']
            #monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']*monc_data[m]['rho']
            #monc_data[m]['model_twc'] = monc_data[m]['model_lwc'] +monc_data[m]['model_iwc']
            monc_data[m]['model_twc'] = monc_data[m]['twc_tot_mean']
            monc_data[m]['model_lwc'] = monc_data[m]['lwc_tot_mean']
            #lwc_zvar=monc_data[m]['zvar']['liquid_mmr_mean']
            lwc_zvar+=[monc_data[m]['zvar']['lwc_tot_mean']]
            lwc_tvar+=[monc_data[m]['tvar']['lwc_tot_mean']]
            nisg_zvar+=[monc_data[m]['zvar']['nisg_tot_mean']]
            nisg_tvar+=[monc_data[m]['tvar']['nisg_tot_mean']]
    if praw==True:
        for m in range(0,len(raw_data)):
            raw_data[m]['qnice'][raw_data[m]['qnice'] < 0] = 0.0

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
                        if um_data[m]['model_twc'][t,k] < twc_thresh_um[k]:
                            um_data[m]['model_twc'][t,k] = np.nan
                            um_data[m]['model_lwc'][t,k] = np.nan
                            if praw == True:
                                raw_data[m]['qice'][t,k] = np.nan       ## apply nan mask to qice based on TWC
                                raw_data[m]['qnice'][t,k] = np.nan       ## apply nan mask to qnice based on TWC
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
                            monc_data[m]['model_lwc'][t,k] = np.nan
                            monc_data[m]['nisg_tot_mean'][t,k] = np.nan

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
    plt.rc('legend',fontsize=SMALL_SIZE)
    # plt.subplots_adjust(top = 0.95, bottom = 0.12, right = 0.95, left = 0.15,
    #         hspace = 0.4, wspace = 0.1)
    ###define colors
    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    ### define axis instance
    ####LWC
    plt.figure(figsize=(14,7))
    plt.subplots_adjust(top = 0.8, bottom = 0.1, right = 0.98, left = 0.08)

    for pt in range(0,len(prof_time)):
        plt.subplot(1,len(prof_time),pt+1)
        ax1 = plt.gca()
        sstr=datenum2date(prof_time[pt][0])
        estr=datenum2date(prof_time[pt][1])
        plt.title(sstr.strftime("%H") +'-' + estr.strftime("%H") + ' UTC')
        # obsid= np.squeeze(np.argwhere((raw_data['time']>=prof_time[pt][0]) & (raw_data['time']<prof_time[pt][1])))
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3,np.nanmean(obs_data['height'],0), color = 'k', linewidth = 3, label = 'Obs_UMgrid'  + lwcstr, zorder = obs_zorder)
        # ax1.fill_betweenx(np.nanmean(obs_data['height'],0),np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3,
        #     np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, color = 'lightgrey', alpha = 0.5)
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 - np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)
        # plt.plot(np.nanmean(obs_data[lwcvar][obsid,:],0)*1e3 + np.nanstd(obs_data[lwcvar][obsid,:],0)*1e3, np.nanmean(obs_data['height'],0),
        #     '--', color = 'k', linewidth = 0.5)
        # if praw==True:
        #     for m in range(0,len(raw_data)):
        #         id=  np.squeeze(np.argwhere((raw_data[m]['time']>=prof_time[pt][0]) & (raw_data[m]['time']<prof_time[pt][1])))
        #         ax1.fill_betweenx(raw_data[m]['height'],np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 - np.nanstd(raw_data[m]['qnliq'][id,:]/1e6,0),
        #             np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 + np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, color = fcols[m], alpha = 0.05)
        #         plt.plot(np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 - np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, raw_data[m]['height'],
        #             '--', color =lcols[m], linewidth = 0.5)
        #         plt.plot(np.nanmean(raw_data[m]['qnliq'][id,:],0)/1e6 + np.nanstd(raw_data[m]['qnliq'][id,:],0)/1e6, raw_data[m]['height'],
        #             '--', color = lcols[m], linewidth = 0.5)
        # if pmonc==True:
        #     for m in range(0,len(monc_data)):
        #         id= np.squeeze(np.argwhere((monc_data[m][ndrop_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][ndrop_tvar[m]]<prof_time[pt][1])))
        #         ax1.fill_betweenx(monc_data[m][ndrop_zvar[m]],np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 - np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:]/1e6,0),
        #             np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 + np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, color = fcolsmonc[m], alpha = 0.05)
        #         plt.plot(np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 - np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, monc_data[m][ndrop_zvar[m]],
        #             '--', color =lcolsmonc[m], linewidth = 0.5)
        #         plt.plot(np.nanmean(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6 + np.nanstd(monc_data[m]['ndrop_tot_mean'][id,:],0)/1e6, monc_data[m][ndrop_zvar[m]],
        #             '--', color = lcolsmonc[m], linewidth = 0.5)
        if praw==True:
            for m in range(0,len(raw_data)):
                id= np.squeeze(np.argwhere((raw_data[m]['time']>=prof_time[pt][0]) & (raw_data[m]['time']<prof_time[pt][1])))
                plt.plot(np.nanmean(raw_data[m]['qnice'][id,:],0)/1e6,raw_data[m]['height'], color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][nisg_tvar[m]]>=prof_time[pt][0]) & (monc_data[m][nisg_tvar[m]]<prof_time[pt][1])))
                plt.plot(np.nanmean(monc_data[m]['nisg_tot_mean'][id,:],0)/1e6,monc_data[m][nisg_zvar[m]], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)
        if pt == 1:
            plt.legend(bbox_to_anchor=(1.5, 1.05), loc=4, ncol=2)


        plt.xlabel('N$_{isg}$ [L$^{-1}$]')
        plt.ylabel('Z [km]')
        # plt.ylim([0,5e3])
        # plt.yticks(np.arange(0,5.01e3,0.5e3))
        # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.set_yticklabels(ytlabels)
        plt.xlim([0,0.005])
        # plt.xticks(np.arange(0,0.001,0.006))
        # #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
        # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.025))


    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if thresholding == True:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_Nisg-MTThresh' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_Nisg-MTThresh'+ lwcstr + '_split.png'
    else:
        if pmonc==True:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_Nisg' + lwcstr + '_split.png'
        else:
            fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_Nisg'+ lwcstr + '_split.png'

    plt.savefig(fileout, dpi=300)
    plt.close()
    # plt.show()
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')

def plot_lwp(obs_data, plots_out_dir, dates,**args ):

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


    print ('******')
    print ('')
    print ('Plotting  timeseries of lwp:')
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
    plt.subplots_adjust(top = 0.95, bottom = 0.05, right = 0.95, left = 0.05,
            hspace = 0.4, wspace = 0.13)

    lcols=['lightseagreen','steelblue','royalblue','darkblue']
    fcols=['lightcyan','lightblue','skyblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    #from IPython import embed;
    fig = plt.figure(figsize=(18,10 ))
    #ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    ax  = fig.add_axes([0.07,0.7,0.7,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(obs_data['time'], obs_data['lwp'][:,0], color = 'black', label = 'ice_station')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    if pum==True:
        for m in range(0,len(um_data)):
            plt.plot(um_data[m]['time'], um_data[m]['model_lwp'], color = lcols[m], label = label[m])
    if pmonc==True:
            for m in range(0,len(monc_data)):
                plt.plot(monc_data[m][monc_data[m]['tvar']['LWP_mean']], monc_data[m]['LWP_mean'], color = lcolsmonc[m], label = mlabel[m])
    plt.ylabel('LWP [kg m$^2$]')

    plt.legend(bbox_to_anchor=(-0.08, 0.77, 1., .102), loc=4, ncol=4)
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))

    #plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim([dates[0], dates[1]])
    plt.xlabel('Time [UTC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    dstr=datenum2date(dates[0])
#    from IPython import embed; embed()
    dstr=datenum2date(dates[1])
    if pmonc==True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_LWP.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_LWP.png'
    #fileout = os.path.join(plots_out_dir,dstr.strftime('%Y%m%d') + '_lwp_ts.png')
    plt.savefig(fileout)


def plot_monc_comparison(obs_data,lwcvar,lwcstr, plots_out_dir, dates, **args): #, lon, lat):

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

    ylims=[0,2.5]
    yticks=np.arange(0,2.5e3,0.5e3)
    ylims=[0,6]
    yticks=np.arange(0,6e3,0.5e3)
    ytlabels=yticks/1e3


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
            monc_data[m]['model_lwc']= monc_data[m]['liquid_mmr_mean']
            monc_data[m]['model_lwc'][monc_data[m]['model_lwc'] <= 0.0] = np.nan
            monc_data[m]['man_lwc1']= monc_data[m]['q_cloud_liquid_mass_mean']
            monc_data[m]['man_lwc2']= monc_data[m]['q_cloud_liquid_mass_mean']*monc_data[m]['rho_mean']
            monc_data[m]['man_lwc3']= monc_data[m]['lwc_tot_mean']
            lwc_tvar+=[monc_data[m]['tvar']['liquid_mmr_mean']]
            lwc_zvar+=[monc_data[m]['zvar']['liquid_mmr_mean']]
    print ('******')
    print ('')
    print ('Plotting LWC timeseries for whole drift period:')
    print ('')

    cmax=0.3
    clev=np.arange(-0.1,0.01,0.005)

    #####PlotLwc###############################################
    yheight=3
    fig = plt.figure(figsize=(9.5,yheight*numsp+1))
    plt.subplots_adjust(top = 0.92, bottom = 0.06, right = 0.92, left = 0.08,
            hspace = 0.4, wspace = 0.2)


    plt.subplot(3,1,1)
    ax = plt.gca()
    img=plt.contourf(monc_data[m][lwc_tvar[m]]/60/60, np.squeeze(monc_data[m][lwc_zvar[m]][:]),
     np.transpose(monc_data[m]['model_lwc'])*1e3- np.transpose(monc_data[m]['man_lwc1'])*1e3,
    levels=clev,cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    ax.set_yticklabels(ytlabels)
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.xlabel('Time (UTC)')
    plt.title('liquid_mmr_mean - q_cloud_liquid_mass_mean')
    ax2 = ax.twinx()
    ax2.set_ylabel('model_lwc', rotation = 270, labelpad = 27)
    ax2.set_yticks([])

    plt.subplot(3,1,2)
    ax = plt.gca()
    plt.contourf(monc_data[m][lwc_tvar[m]]/60/60, np.squeeze(monc_data[m][lwc_zvar[m]][:]),
     np.transpose(monc_data[m]['man_lwc2'])*1e3- np.transpose(monc_data[m]['man_lwc3'])*1e3,
    levels=clev,cmap = newcmp)
    plt.ylabel('Z [km]')
    plt.ylim(ylims)
    plt.yticks(yticks)
    plt.title('q_cloud_liquid_mass_mean*rho_mean - lwc_tot_mean')
    ax.set_yticklabels(ytlabels)
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.xlabel('Time (UTC)')
    ax2 = ax.twinx()
    ax2.set_ylabel('q_liquid_mean', rotation = 270, labelpad = 27)
    ax2.set_yticks([])
    cbaxes = fig.add_axes([0.225, 0.25, 0.6, 0.015])
    cb = plt.colorbar(img, cax = cbaxes, orientation = 'horizontal')

    dstr=datenum2date(dates[1])
    if pmonc == True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_MONC_LWCTimeseries_comparison'+ lwcstr + '.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_MONC_LWCTimeseries_comparison'+ lwcstr + '.png'
    plt.savefig(fileout)
    plt.close()

    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')


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
    print ('')

    ### build nanmask
    nanmask = np.zeros([np.size(obs_data['time']), np.size(obs_data['Cv'],1)])
    nanindex = np.zeros([np.size(obs_data['time'])])
    lwcindex = np.zeros([np.size(obs_data['time'])])
    iwcindex = np.zeros([np.size(obs_data['time'])])
    wc0index = np.zeros([np.size(obs_data['time'])])
    lwpindex = np.zeros([np.size(obs_data['time'])])
    for i in range(len(obs_data['time'])):
        if np.isnan(np.nanmean(obs_data['Cv'][i,:], 0)):       ## if there are only nans in the time profile
            nanmask[i,:] = 1.0
            #nanmask[i-1,:] = 1.0   ##why setting -1 and +1 to nan???
            #nanmask[i+1,:] = 1.0
            nanindex[i] = 1
        # if np.logical_and(np.isnan(np.nanmean(obs_data['lwc_adiabatic'][i,:], 0)), np.isnan(np.nanmean(obs_data['iwc'][i,:], 0))):       ## if both wc profiles contain only nans
            # wcindex[i] = 1
        if np.isnan(np.nanmean(obs_data['lwc'][i,:], 0)):
            lwcindex[i] = 1
        if np.isnan(np.nanmean(obs_data['iwc'][i,:], 0)):
            iwcindex[i] = 1
        elif np.logical_and(np.nanmean(obs_data['lwc'][i,:], 0) == 0, np.nanmean(obs_data['iwc'][i,:], 0) == 0):       ## if both wc profiles contain only zeros
            wc0index[i] = 1
        elif np.isnan(obs_data['lwp'][i,0]):       ## if there are only nans in the lwp timeseries
            lwpindex[i] = 1
    nanmask[nanmask == 0.0] = np.nan
    nanind = np.where(nanindex == 1)
    lwcind = np.where(lwcindex == 1)
    iwcind = np.where(iwcindex == 1)
    wc0ind = np.where(wc0index == 1)
    lwpind = np.where(lwpindex == 1)

    return nanind, nanmask, lwcind, iwcind, wc0ind, lwpind
    print('done')
    print('**********')

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
    print ('done')
    print ('*******')
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

def CaseStudySelection(obs_data, um_data, monc_data, raw_data, dates):
    print('')
    print('************')
    print('Shortening data to case study time frame:')

    obs_vars=list(obs_data.keys())
    time_ind = np.argwhere((obs_data['time']<dates[0]) | (obs_data['time']>dates[1]+1/60/24 ))
    tinit=len(obs_data['time'])
    for j in range(0,len(obs_vars)):
        tmp2=np.argwhere(np.array(obs_data[obs_vars[j]].shape) == tinit)
        if tmp2 == 0:
            obs_data[obs_vars[j]]=np.delete(obs_data[obs_vars[j]],time_ind,0)
        elif tmp2 == 1:
            obs_data[obs_vars[j]]=np.delete(obs_data[obs_vars[j]],time_ind,1)
        elif tmp2 == 2:
            obs_data[obs_vars[j]]=np.delete(obs_data[obs_vars[j]],time_ind,2)
        elif tmp2 == 3:
            obs_data[obs_vars[j]]=np.delete(obs_data[obs_vars[j]],time_ind,3)

    for m in range(0,len(um_data)):
        um_vars=list(um_data[m].keys())
        time_ind = np.argwhere((um_data[m]['time']<dates[0]) | (um_data[m]['time']>dates[1]+1/60/24 ))
        tinit=len(um_data[m]['time'])
        for j in range(0,len(um_vars)):
            tmp2=np.argwhere(np.array(um_data[m][um_vars[j]].shape) == tinit)
            if tmp2 == 0:
                um_data[m][um_vars[j]]=np.delete(um_data[m][um_vars[j]],time_ind,0)
            elif tmp2 == 1:
                um_data[m][um_vars[j]]=np.delete(um_data[m][um_vars[j]],time_ind,1)
            elif tmp2 == 2:
                um_data[m][um_vars[j]]=np.delete(um_data[m][um_vars[j]],time_ind,2)
            elif tmp2 == 3:
                um_data[m][um_vars[j]]=np.delete(um_data[m][um_vars[j]],time_ind,3)

    for m in range(0,len(raw_data)):
        raw_vars=list(raw_data[m].keys())
        time_ind = np.argwhere((raw_data[m]['time']<dates[0]) | (raw_data[m]['time']>dates[1]+1/60/24 ))
        tinit=len(raw_data[m]['time'])
        for j in range(0,len(raw_vars)):
            tmp2=np.argwhere(np.array(raw_data[m][raw_vars[j]].shape) == tinit)
            if tmp2 == 0:
                raw_data[m][raw_vars[j]]=np.delete(raw_data[m][raw_vars[j]],time_ind,0)
            elif tmp2 == 1:
                raw_data[m][raw_vars[j]]=np.delete(raw_data[m][raw_vars[j]],time_ind,1)
            elif tmp2 == 2:
                raw_data[m][raw_vars[j]]=np.delete(raw_data[m][raw_vars[j]],time_ind,2)
            elif tmp2 == 3:
                raw_data[m][raw_vars[j]]=np.delete(raw_data[m][raw_vars[j]],time_ind,3)

    for m in range(0,len(monc_data)):
        monc_vars=list(monc_data[m].keys())
        monc_vars.remove('tvar')
        monc_vars.remove('zvar')
        mt1=monc_data[m]['time1']
        mt2=monc_data[m]['time2']

        # mt1=monc_data[m]['time1']/60/60/24 + dates[0]
        # mt2=monc_data[m]['time2']/60/60/24 + dates[0]
        # mt3=monc_data[m]['time3']/60/60/24 + dates[0]
        time1_ind = np.argwhere((mt1>dates[1] ))
        t1init=len(mt1)
        time2_ind = np.argwhere((mt2>dates[1] ))
        t2init=len(mt2)
        if 'time3' in monc_data:
            mt3=monc_data[m]['time3']
            time3_ind = np.argwhere((mt3>dates[1] ))
            t3init=len(mt3)
        for j in range(0,len(monc_vars)):
            if any(np.array(monc_data[m][monc_vars[j]].shape) == t1init):
                monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time1_ind,0)
            elif any(np.array(monc_data[m][monc_vars[j]].shape) == t2init):
                tmp2=np.argwhere(np.array(monc_data[m][monc_vars[j]].shape) == t2init)
                if tmp2 == 0:
                    monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time2_ind,0)
                elif tmp2 == 1:
                    monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time2_ind,1)
                elif tmp2 == 2:
                    monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time2_ind,2)
                elif tmp2 == 3:
                    monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time2_ind,3)
            elif 'time3' in monc_data:
                if any(np.array(monc_data[m][monc_vars[j]].shape) == t3init):
                    tmp2=np.argwhere(np.array(monc_data[m][monc_vars[j]].shape) == t2init)
                    if tmp2 == 0:
                        monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time3_ind,0)
                    elif tmp2 == 1:
                        monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time3_ind,1)
                    elif tmp2 == 2:
                        monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time3_ind,2)
                    elif tmp2 == 3:
                        monc_data[m][monc_vars[j]]=np.delete(monc_data[m][monc_vars[j]],time3_ind,3)

    return obs_data, um_data, monc_data, raw_data
    print('done')
    print('************')

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
        cn_obs_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/Cloudnet/Observations/measurements_V7/'
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

    ### -----------------------------------------------------------------
    ### CHOOSE UM RUNS - MODEL DATA
    #out_dir = ['23_u-cc278_RA1M_CASIM/',
    #           '24_u-cc324_RA2T_CON/',
    #           '25_u-cc568_RA2M_CON/']
    out_dir = ['23_u-cc278_RA1M_CASIM/',
               '30_u-cg179_RA1M_CASIM/',
               '26_u-cd847_RA1M_CASIM/',
               '27_u-ce112_RA1M_CASIM/',
               ]

    # out_dir = ['23_u-cc278_RA1M_CASIM/']
    ### CHOOSE MONC RUNS
    # m_out_dir = ['5_control_20180913T0000Z_Wsub-1.5_Fletcher/',
    #              '6_control_20180913T0000Z_Wsub-1.5-1km/',
    #              '8_control_20180913T0000Z_Wsub-1.0-1km/',
    #              '9_control_20180913T0000Z_Wsub-0.5-1km/']
    m_out_dir =[#'22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/',
               # '23_20180913T0000Z_6hSpin-up_12h0600-0000thTend_20h1200-0600thTend/',
               # '24_20180913T0000Z_12h0600-0000_20h1200-0600thTend_0.1Cooper/',
               # '25_20180913T0000Z_20h0600-0000thTend/'
               # '26A_20180913T0000Z_6hSpinUp_12h0600-0000thTend_20h1200-0600thTend_6-20h0.1Cooper/',
               # '26B_20180913T0000Z_6hSpinUp_12h0600-0000thTend_20h1200-0600thTend_6-20h0-Cooper/',
               # '27A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24hCooper/',
               # '27B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.5Cooper/',
               '27C_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper/',
               '27D_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd25/',
               '27E_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd10/',
               '27F_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd5/',
               '30A_20180913T0000Z_8hSpinUp_8-14hUVRelax0600_14-24hUVRelax1200_8-24h0.1Cooper_FixedNd10/']
               # '28A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-100-ARG/',
               # '28B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-100-Twomey/'
               # '29A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-20-ARG/',
               # '29B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-20-allAct/',
            #'4_control_20180913T0000Z_Wsub-1.5/',
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
        elif out_dir[m][:10] == '30_u-cg179':
            label.append('UM_CASIM-100-PasProc')
            outstr.append('CASIM100-PasProc')
        else:
            label.append('undefined_label')
            outstr.append('')
            print(label)
    mlabel=[]
    moutstr=[]
    for m in range(0, len(m_out_dir)):
        # if m_out_dir[m][:1] == '3':
        #     mlabel.append('MONC nosub')
        #     moutstr.append('Mnowsub')
        # elif m_out_dir[m][:1] == '4':
        #     mlabel.append('MONC Wsub1.5')
        #     moutstr.append('Mwsub')
        # elif m_out_dir[m][:1] == '5':
        #     mlabel.append('MONC Wsub1.5 \n Fletcher')
        #     moutstr.append('Mwsubfle')
        # elif m_out_dir[m][:1] == '6':
        #     mlabel.append('MONC Wsub1.5-1km')
        #     moutstr.append('Mwsub1.5-1km')
        # elif m_out_dir[m][:1] == '7':
        #     mlabel.append('MONC Wsub1.5-1km \n solACC-100')
        #     moutstr.append('Mwsub1kmsolACC100')
        # elif m_out_dir[m][:1] == '8':
        #     mlabel.append('MONC Wsub1.0-1km')
        #     moutstr.append('Mwsub1.0-1km')
        # elif m_out_dir[m][:1] == '9':
        #     mlabel.append('MONC Wsub0.5-1km')
        #     moutstr.append('Mwsub0.5-1km')
        if m_out_dir[m][:2] == '20':
            mlabel.append('MONC qinit2 800m \n thqvTend noice')
            moutstr.append('qin2-thqvTend-noice')
        elif m_out_dir[m][:2] == '22':
            mlabel.append('MONC qinit2 800m \n thForcing-0000-0600')
            moutstr.append('MONC-22')
        elif m_out_dir[m][:2] == '23':
            mlabel.append('MONC thForcing-0600-0000')
            moutstr.append('MONC-23')
        elif m_out_dir[m][:2] == '24':
            mlabel.append('MONC thForcing-12h0600-0000-20h1200-0600 0.1*Cooper')
            moutstr.append('MONC-24')
        elif m_out_dir[m][:2] == '25':
            mlabel.append('MONC thForcing-20h0600-0000')
            moutstr.append('MONC-25')
        elif m_out_dir[m][:3] == '26A':
            mlabel.append('MONC 6hSpinUp thForcing-0-12h0600-0000 6-20h-0.1*Cooper')
            moutstr.append('MONC-26A')
        elif m_out_dir[m][:3] == '26B':
            mlabel.append('MONC 6hSpinUp thForcing-0-12h0600-0000 6-20h-0*Cooper')
            moutstr.append('MONC-26B')
        elif m_out_dir[m][:3] == '27A':
            mlabel.append('MONC_Cooper_FixedNd50')
            moutstr.append('MONC-27A')
        elif m_out_dir[m][:3] == '27B':
            mlabel.append('MONC_0.5Cooper_FixedNd50')
            moutstr.append('MONC-27B')
        elif m_out_dir[m][:3] == '27C':
            mlabel.append('MONC_0.1Cooper_FixedNd50')
            moutstr.append('MONC-27C')
        elif m_out_dir[m][:3] == '27D':
            mlabel.append('MONC_0.1Cooper_FixedNd25')
            moutstr.append('MONC-27D')
        elif m_out_dir[m][:3] == '27E':
            mlabel.append('MONC_0.1Cooper_FixedNd10')
            moutstr.append('MONC-27E')
        elif m_out_dir[m][:3] == '27F':
            mlabel.append('MONC_0.1Cooper_FixedNd5')
            moutstr.append('MONC-27F')
        elif m_out_dir[m][:3] == '28A':
            mlabel.append('MONC_0.1Cooper_CASIM-100-ARG')
            moutstr.append('MONC-28A')
        elif m_out_dir[m][:3] == '28B':
            mlabel.append('MONC_0.1Cooper_CASIM-100-Twomey')
            moutstr.append('MONC-28B')
        elif m_out_dir[m][:3] == '29A':
            mlabel.append('MONC_0.1Cooper_CASIM-20-ARG')
            moutstr.append('MONC-29A')
        elif m_out_dir[m][:3] == '29B':
            mlabel.append('MONC_0.1Cooper_CASIM-20-allAct')
            moutstr.append('MONC-29B')
        elif m_out_dir[m][:3] == '30A':
            mlabel.append('MONC_0.1Cooper \n UVRelax FixedNd10')
            moutstr.append('MONC-30A')
        else:
            label.append('undefined_label')
            moutstr.append('')

    #---- MONC SPIN UP TIME
    spin6 = ['26']
    spin8 = ['27','28','29','30']

    if m_out_dir[0][:2] in spin6:
        monc_spin = 6 *60 *60
    elif m_out_dir[0][:2] in spin8:
        monc_spin = 8 *60 *60

    ### -----------------------------------------------------------------
    ### create monc filenames
    monc_filename=[]
    for m in range(0, len(m_out_dir)):
        fname=glob.glob(monc_root_dir + m_out_dir[m] +'*dg*.nc')
        monc_filename.append(fname)

    monc_3d_filename=[]
    for m in range(0, len(m_out_dir)):
        fname=glob.glob(monc_avg_dir + m_out_dir[m] +'3d*npy')
        monc_3d_filename.append(fname)

    ### -----------------------------------------------------------------
    ### CHOSEN RUN - CLOUDNET DATA
    cn_um_out_dir = ['cloud-fraction-metum-grid/2018/',
                     'lwc-scaled-metum-grid/2018/',
                     'iwc-Z-T-metum-grid/2018/']
    cn_obs_out_dir = ['cloud-fraction-metum-grid/2018/qual70/',
                        'lwc-adiabatic-metum-grid/2018/',
                        'iwc-Z-T-metum-grid/2018/']

    um_out_dir = 'OUT_R0/'

    # # -------------------------------------------------------------
    # # Load cloudnet data
    # # -------------------------------------------------------------
    print ('******')
    print ('')
    print ('Begin nc read in: ')
    print (' ')

    for i in range(0,len(names)):
        dstr=names[i][0:8]
        datenum = date2datenum(dtime.datetime.strptime(dstr,'%Y%m%d'))
        cn_filename_um=[]
        filename_um=[]
        for m in range(0,len(out_dir)):
            cn_filename_um.append([cn_um_dir + out_dir[m] + cn_um_out_dir[0] + names[i] + cn_um_out_dir[0][-31:-6] + '.nc',
                        cn_um_dir + out_dir[m] + cn_um_out_dir[1] + names[i] + cn_um_out_dir[1][-27:-6] + '.nc',
                        cn_um_dir + out_dir[m] + cn_um_out_dir[2] + names[i] + cn_um_out_dir[2][-24:-6] + '.nc'])
            filename_um.append([um_root_dir + out_dir[m] + um_out_dir + names[i] + 'metum.nc'])
        cn_filename_obs = [cn_obs_dir + cn_obs_out_dir[0] + names[i] + cn_obs_out_dir[0][:-13] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[1] + names[i] + cn_obs_out_dir[1][:-6] + '.nc',
                        cn_obs_dir + cn_obs_out_dir[2] + names[i] + cn_obs_out_dir[2][:-6] + '.nc']

    ### --------------------------------------------------------------------
    ###     READ IN ALL CLOUDNET FILES: reinitialise diagnostic dictionaries
    ### --------------------------------------------------------------------
        print ('** LOAD DAY:', i)
        print ('   Loading multiple diagnostics from cloudnet results:')
        cn_nc0 = {} # observations
        for c in range(0,3):
            cn_nc0[c] = Dataset(cn_filename_obs[c],'r')
        cn_nc1={}
        nc1 = {}
        for m in range(0,len(out_dir)): #UM model data
            cn_nc1[m]={}
            for c in range(0,3):
                print(cn_filename_um[m][c])
                cn_nc1[m][c] = Dataset(cn_filename_um[m][c],'r')
            nc1[m] = {}
            print (filename_um[m])
            nc1[m] = Dataset(filename_um[m][0],'r')

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

        print ('OBS data loaded!')

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

        print ('UM Cloudnet data loaded!')
        print ('')

        ### --------------------------------------------------------------------
        ###     LOAD RAW UM MODEL DATA INTO DICTIONARY
        ### --------------------------------------------------------------------

        var_list = ['height','temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
                    'air_temperature_at_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','uwind','vwind','wwind',
                    'cloud_fraction','radr_refl','qnliq','qnice','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation', 'latent_heat_flux',
                    'toa_outgoing_longwave_flux','toa_incoming_shortwave_flux','toa_outgoing_shortwave_flux','seaice_albedo_agg','qsnow','qice','qicecrystals']

        if i == 0:
            raw_data={}
            time_raw={}

        for m in range(0,len(out_dir)): #UM model data
            if i == 0:
                raw_data[m]={}
                time_raw[m]={}
                ### create time arrays for all model data
                time_raw[m] = datenum + np.float64((nc1[m].variables['forecast_time'][:])/24.0)
                ### load in initial UM data
                for j in range(0,len(var_list)):
                    if np.ndim(nc1[m].variables[var_list[j]]) == 1:  # 1d timeseries only
                        raw_data[m][var_list[j]] = nc1[m].variables[var_list[j]][:]
                    else:                                   # 2d column um_data
                        raw_data[m][var_list[j]] = nc1[m].variables[var_list[j]][:]
            ### fill arrays with remaining data
            else:
                time_raw[m] = np.append(time_raw[m], datenum + np.float64((nc1[m].variables['forecast_time'][:])/24.0))
                ### append rest of UM data
                for j in range(0,len(var_list)):
                    if np.ndim(nc1[m].variables[var_list[j]]) == 1:
                        raw_data[m][var_list[j]] = np.append(raw_data[m][var_list[j]],nc1[m].variables[var_list[j]][:])
                    else:
                        raw_data[m][var_list[j]] = np.append(raw_data[m][var_list[j]],nc1[m].variables[var_list[j]][:],0)

        for m in range(0,len(out_dir)): nc1[m].close()

        print ('UM raw model data loaded!')
        print ('')

    #################################################################
    ## save time to dictionaries now we're not looping over all diags anymore
    #################################################################
    for m in range(0,len(out_dir)):
        um_data[m]['time'] = time_um[m]
        raw_data[m]['time'] = time_raw[m]
    obs_data['time'] = time_obs

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


    monc_var_3d_list =['T_mean','p_mean','th_mean','rho_mean','q_vapour_mean',
                        'q_cloud_liquid_mass_mean','q_ice_mass_mean','q_snow_mass_mean','q_graupel_mass_mean',
                        'q_cloud_liquid_number_mean','q_ice_number_mean','q_snow_number_mean','q_graupel_number_mean',
                        'twc_tot_mean','iwc_tot_mean','nisg_tot_mean','ndrop_tot_mean','lwc_tot_mean']
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
    monc_data=removeSpinUp(monc_data,monc_spin)

    ## -------------------------------------------------------------
    ## convert monc time to datenum
    ## -------------------------------------------------------------
    for m in range(0,len(monc_data)):
        monc_data[m]['time2']=dates[0] + monc_data[m]['time2']/60/60/24
        monc_data[m]['time1']=dates[0] + monc_data[m]['time1']/60/60/24
        if 'time3' in monc_data:
            monc_data[m]['time3']=dates[0] + monc_data[m]['time3']/60/60/24


    ## -------------------------------------------------------------
    ## shorten obs/model data to case study time period
    ## -------------------------------------------------------------

    obs_data, um_data, monc_data, raw_data = CaseStudySelection(obs_data, um_data, monc_data, raw_data, dates)

    ## -------------------------------------------------------------
    ## maximise obs data available and build mask for available data
    ## -------------------------------------------------------------
    print ('*******************')
    print('setting missing data to nan and interpolate missing obs')
    obs_data, um_data = setFlags(obs_data, um_data, obs_var_list, um_var_list)
    obs_data = interpCloudnet(obs_data)
    nanind, nanmask, lwcind, iwcind, wc0ind, lwpind = buildNaNMask(obs_data)
    #nanindadv, nanmaskadv, wcindadv, wc0indadv, lwpindadv = buildNaNMaskadv(obs_data)

    varlist_obs = ['Cv', 'lwc','lwc_adiabatic_inc_nolwp','lwc_adiabatic', 'iwc', 'lwp']
    varlist_um = ['model_Cv_filtered', 'model_lwc', 'model_iwc_filtered', 'model_lwp']
    # print(um_data.keys())

    ## remove missing Cv obs timesteps (remove from all)
    print(' remove missing Cv obs timesteps')
    for c in range(0, 1):
        for m in range(0, len(out_dir)):
            um_data[m][varlist_um[c]][nanind, :] = np.nan

    ### remove missing water content obs timestep (only remove from water contents)
    for c in range(1, 2):
        for m in range(0, len(out_dir)):
            um_data[m][varlist_um[c]][lwcind, :] = np.nan
    for c in range(2, 3):
        for m in range(0, len(out_dir)):
            um_data[m][varlist_um[c]][iwcind, :] = np.nan

    ## lwp only 1d
    for m in range(0, len(out_dir)):
        um_data[m]['model_lwp'][lwpind] = np.nan


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

    # figure = plot_iwcProfiles(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates, um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    # figure = plot_twcProfiles(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_NdropProfiles_split(obs_data, lwcvar,lwcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,raw_data=raw_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_lwcProfiles_split(obs_data, lwcvar,lwcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_NisgProfiles_split(obs_data, lwcvar,lwcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,raw_data=raw_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_iwcProfiles_split(obs_data, twcvar,twcstr,thresholding, plots_out_dir,dates, prof_times,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)


    # -------------------------------------------------------------
    # Cloudnet plot: Plot contour timeseries
    # -------------------------------------------------------------
    # figure = plot_CvTimeseries(obs_data,plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_LWCTimeseries(obs_data,  lwcvar,lwcstr, plots_out_dir, dates, um_data=um_data,label=label, outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_TWCTimeseries( obs_data, twcvar,twcstr, plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr,monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_IWCTimeseries(obs_data, plots_out_dir, dates,  um_data=um_data,label=label,outstr=outstr,monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
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
