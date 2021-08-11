### SCRIPT TO COMPARE MODEL TO MEASURMENTS
from  __future__ import print_function
import time
import datetime as dtime
from IPython import embed
import numpy as np
#import pandas as pd
from netCDF4 import Dataset
#import diags_MOCCHA as diags
#import diags_varnames as varnames
#import cartopy.crs as ccrs
#import iris
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as ticker

#import matplotlib.cm as mpl_cm
import os
import glob
#import seaborn as sns

### import python functions
import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from manipFuncts import intersect_mtlb
#from physFuncts import calcThetaE, calcThetaVL
#from pyFixes import py3_FixNPLoad


def plot_surfaceVariables(obs, plot_out_dir, dates,**args  ):

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
    print ('Plotting  timeseries of surface variables:')
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

    lcols=['mediumseagreen','steelblue','darkblue']
    fcols=['mediumaquamarine','lightblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    #from IPython import embed; embed()
    fig = plt.figure(figsize=(18,10 ))
    #ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    ax  = fig.add_axes([0.07,0.7,0.7,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(obs['metalley']['mday'], obs['metalley']['t'], color = 'black', label = 'ice_station')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    if pum==True:
        for m in range(0,len(um_data)):
            plt.plot(um_data[m]['time'], um_data[m]['air_temperature_at_1.5m']-273.15, color = lcols[m], label = label[m])
    if pmonc==True:
        for m in range(0,len(monc_data)):
            plt.plot(monc_data[m]['time'], monc_data[m]['air_temperature_at_1.5m']-273.15, color = lcolsmonc[m], label = mlabel[m])
    plt.ylabel('T [$^\circ$C]')
    plt.legend(bbox_to_anchor=(-0.08, 0.77, 1., .102), loc=4, ncol=4)
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.ylim([-10,0])

    ax  = fig.add_axes([0.07,0.4,0.7,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(data1['time'], data1['rh_1.5m'], color = 'darkblue', label = label1)
    plt.plot(data3['time'], data3['rh_1.5m'], color = 'steelblue', label = label3[:-4])
    plt.plot(data2['time'], data2['rh_1.5m'], color = 'mediumseagreen', label = label2)
    plt.plot(obs['metalley']['mday'], obs['metalley']['rh'], color = 'black', label = 'Obs')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    plt.ylabel('RH [%]')
    #plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim(dates[0], dates[1])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.ylim([80,110])
    plt.xlabel('Time [UTC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    date=datenum2date(dates[0])
#    from IPython import embed; embed()
    fileout = os.path.join(plot_out_dir,date.strftime('%Y%m%d') + '_surfaceVariables_ts.png')
    plt.savefig(fileout)

def plot_lwp(obs_data, plot_out_dir, dates,**args ):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
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

    lcols=['mediumseagreen','steelblue','darkblue']
    fcols=['mediumaquamarine','lightblue','blue']
    lcolsmonc=['gold','darkgoldenrod','darkorange','orangered','firebrick']
    fcolsmonc=['navajowhite','goldenrod','moccasin','lightsalmon','lightcoral']
    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    #from IPython import embed; embed()
    fig = plt.figure(figsize=(18,10 ))
    #ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    ax  = fig.add_axes([0.07,0.7,0.7,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(obs_data['hatpro']['mday'], obs_data['hatpro']['lwp']/1e3, color = 'black', label = 'ice_station')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    for m in range(0,len(um_data)):
        plt.plot(um_data[m]['time'], um_data[m]['LWP']-273.15,'^', color = lcols[m], label = label[m])
    for m in range(0,len(monc_data)):
        plt.plot(monc_data[m][monc_data[m]['tvar']['LWP_mean']], monc_data[m]['LWP_mean']-273.15,'o', color = lcolsmonc[m], label = mlabel[m])
    plt.ylabel('LWP [g m$^2$]')

    plt.legend(bbox_to_anchor=(-0.08, 0.77, 1., .102), loc=4, ncol=4)
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))

    #plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim([datenum, edatenum])
    plt.xlabel('Time [UTC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    date=datenum2date(datenum)
#    from IPython import embed; embed()
    fileout = os.path.join(plot_out_dir,date.strftime('%Y%m%d') + '_lwp_ts.png')
    plt.savefig(fileout)


def plot_BLDepth_SMLDepth(obs_data, plot_out_dir, dates,**args ):

    numsp=1
    if bool(args):
        for n in range(0,len(args)):
            if  list(args.keys())[n] == 'monc_data':
                monc_data=args[list(args.keys())[n]]
                numsp += len(monc_data)
            elif list(args.keys())[n] == 'mlabel':
                mlabel = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'moutstr':
                moutstr= args[list(args.keys())[n]]
            elif  list(args.keys())[n] == 'um_data':
                um_data=args[list(args.keys())[n]]
                numsp += len(um_data)
            elif list(args.keys())[n] == 'label':
                label = args[list(args.keys())[n]]
            elif list(args.keys())[n] == 'outstr':
                outstr= args[list(args.keys())[n]]

    print ('******')
    print ('')
    print ('Plotting  timeseries of BLDepth and decoupling:')
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
    #from IPython import embed; embed()
    fig = plt.figure(figsize=(10,12 ))
    ax  = fig.add_axes([0.1,0.7,0.8,0.2])   # left, bottom, width, height
    yB = [-10, 120]
    lines=[]
    ax.fill_between(np.squeeze(obs_data['dec']['mday']),np.squeeze(obs_data['dec']['cbase_sandeep']),np.squeeze(obs_data['dec']['ct']), color = 'skyblue', alpha = 0.3)
    lines +=plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['invbase'], color = 'k', label = 'OBS inv')
    lines +=plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['decbase'],'-', color = 'gray', label = 'OBS sml')
    #legend1=plt.legend(loc='best')
    for m in range(0,len(um_data)):
        #plt.plot(um_data[m]['time'], um_data[m]['bl_depth'], 'o',color = lcols[m], label = 'model bl')
        lines +=plt.plot(um_data[m]['inv']['mday'], um_data[m]['inv']['invbase'], color = lcols[m], label = 'invbase')
    for m in range(0,len(monc_data)):
        lines +=plt.plot(monc_data[m]['inv']['mday'], monc_data[m]['inv']['invbase'], color = lcolsmonc[m], label = 'invbase')
    plt.ylabel('Height [m]')
    llabels= ['Obs inv','OBS sml']
    llabels+=label
    llabels+=mlabel
    fig.legend(lines,llabels, bbox_to_anchor=(0.1, 0.93, 0.8, .05),loc=3, ncol=4)
    #plt.gca().add_artist(legend1)
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('Main inversion')
    ax  = fig.add_axes([0.1,0.4,0.8,0.2])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    lines=[]
    ax.fill_between(np.squeeze(obs_data['dec']['mday']),np.squeeze(obs_data['dec']['cbase_sandeep']),np.squeeze(obs_data['dec']['ct']), color = 'skyblue', alpha = 0.3)
    lines+=plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['invbase'], color = 'k', label = 'inv')
    lines+= plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['decbase'],'-', color = 'gray', label = 'sml')
    #legend1=plt.legend(loc='best')
    for m in range(0,len(um_data)):
        #plt.plot(um_data[m]['time'], um_data[m]['bl_depth'], 'o',color = lcols[m], label = 'model bl')
        lines+=plt.plot(um_data[m]['inv']['mday'], um_data[m]['inv']['decbase'],'-', color = lcols[m], label = 'sml')
        #if m ==0: legend1=plt.legend(loc='best')
    for m in range(0,len(monc_data)):
        lines+=plt.plot(monc_data[m]['inv']['mday'], monc_data[m]['inv']['decbase'],'-', color = lcolsmonc[m], label = 'sml')
    plt.ylabel('Height [m]')
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('sml height')
    #ax.set_xlim([datenum, edatenum])

    ax  = fig.add_axes([0.1,0.1,0.8,0.2])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    lines=[]
    ax.fill_between(np.squeeze(obs_data['dec']['mday']),np.squeeze(obs_data['dec']['cbase_sandeep']),np.squeeze(obs_data['dec']['ct']), color = 'skyblue', alpha = 0.3)
    plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['invbase'], color = 'k', label = 'inv')
    lines+= plt.plot(obs_data['hatpro_temp']['mday'], obs_data['hatpro_temp']['decbase'],'-', color = 'gray', label = 'sml')
    for m in range(0,len(um_data)):
        lines+=plt.plot(um_data[m]['time'], um_data[m]['bl_depth'], '-',color = lcols[m], label = 'model bl')
        #if m ==0: legend1=plt.legend(loc='best')
    plt.ylabel('Height [m]')
    ax.set_xlim([dates[0], dates[1]])
    plt.grid(which='both')
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.title('Model BL height')
    plt.xlabel('Time [UTC]')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    date=datenum2date(dates[0])
#    from IPython import embed; embed()
    fileout = os.path.join(plot_out_dir + dstr.strftime('%Y%m%d') + '_Obs' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_BLdepth-SML.png.png')
    plt.savefig(fileout)



def plot_T_profiles_split(obs, plots_out_dir,dates,prof_time, **args): #, lon, lat):
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
    print ('Plotting T mean profiles split times:')
    print ('')

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
    plt.figure(figsize=(18,8))
    plt.subplots_adjust(top = 0.8, bottom = 0.1, right = 0.92, left = 0.08)
    for pt in range(0,len(prof_time)):
        plt.subplot(1,len(prof_time),pt+1)
        ax1 = plt.gca()
        sstr=datenum2date(prof_time[pt][0])
        estr=datenum2date(prof_time[pt][1])
        plt.title(sstr.strftime("%H") +'-' + estr.strftime("%H") + ' UTC')
        obsid= np.squeeze(np.argwhere((obs['hatpro_temp']['mday']>=prof_time[pt][0]) & (obs['hatpro_temp']['mday']<prof_time[pt][1])))
        plt.plot(np.nanmean(obs['hatpro_temp']['temperature'][:,obsid],1),obs['hatpro_temp']['Z'], color = 'k', linewidth = 3, label = 'HATPRO', zorder = obs_zorder)
        ax1.fill_betweenx(obs['hatpro_temp']['Z'],np.nanmean(obs['hatpro_temp']['temperature'][:,obsid],1) - np.nanstd(obs['hatpro_temp']['temperature'][:,obsid],1),
            np.nanmean(obs['hatpro_temp']['temperature'][:,obsid],1) + np.nanstd(obs['hatpro_temp']['temperature'][:,obsid],1), color = 'lightgrey', alpha = 0.5)
        # plt.xlim([0,0.2])
        plt.plot(np.nanmean(obs['hatpro_temp']['temperature'][:,obsid],1) - np.nanstd(obs['hatpro_temp']['temperature'][:,obsid],1),obs['hatpro_temp']['Z'],
            '--', color = 'k', linewidth = 0.5)
        plt.plot(np.nanmean(obs['hatpro_temp']['temperature'][:,obsid],1) + np.nanstd(obs['hatpro_temp']['temperature'][:,obsid],1), obs['hatpro_temp']['Z'],
            '--', color = 'k', linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id=  np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                ax1.fill_betweenx(um_data[m]['height'],np.nanmean(um_data[m]['temperature'][id,:],0) - np.nanstd(um_data[m]['temperature'][id,:],0),
                    np.nanmean(um_data[m]['temperature'][id,:],0) + np.nanstd(um_data[m]['temperature'][id,:],0), color = fcols[m], alpha = 0.05)
                plt.plot(np.nanmean(um_data[m]['temperature'][id,:],0) - np.nanstd(um_data[m]['temperature'][id,:],0), um_data[m]['height'],
                    '--', color =lcols[m], linewidth = 0.5)
                plt.plot(np.nanmean(um_data[m]['temperature'][id,:],0) + np.nanstd(um_data[m]['temperature'][id,:],0),um_data[m]['height'],
                    '--', color = lcols[m], linewidth = 0.5)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                tvar=monc_data[m]['tvar']['T_mean']
                zvar=monc_data[m]['zvar']['T_mean']
                id= np.squeeze(np.argwhere((monc_data[m][tvar]>=prof_time[pt][0]) & (monc_data[m][tvar]<prof_time[pt][1])))
                ax1.fill_betweenx(monc_data[m][zvar],np.nanmean(monc_data[m]['T_mean'][id,:],0) - np.nanstd(monc_data[m]['T_mean'][id,:],0),
                    np.nanmean(monc_data[m]['T_mean'][id,:],0) + np.nanstd(monc_data[m]['T_mean'][id,:],0), color = fcolsmonc[m], alpha = 0.05)
                plt.plot(np.nanmean(monc_data[m]['T_mean'][id,:],0) - np.nanstd(monc_data[m]['T_mean'][id,:],0), monc_data[m][zvar],
                    '--', color =lcolsmonc[m], linewidth = 0.5)
                plt.plot(np.nanmean(monc_data[m]['T_mean'][id,:],0) + np.nanstd(monc_data[m]['T_mean'][id,:],0), monc_data[m][zvar],
                    '--', color = lcolsmonc[m], linewidth = 0.5)
        if pum==True:
            for m in range(0,len(um_data)):
                id= np.squeeze(np.argwhere((um_data[m]['time']>=prof_time[pt][0]) & (um_data[m]['time']<prof_time[pt][1])))
                plt.plot(np.nanmean(um_data[m]['temperature'][id,:],0),um_data[m]['height'], color = lcols[m], linewidth = 3, label = label[m], zorder = 1)
        if pmonc==True:
            for m in range(0,len(monc_data)):
                id= np.squeeze(np.argwhere((monc_data[m][tvar]>=prof_time[pt][0]) & (monc_data[m][tvar]<prof_time[pt][1])))
                plt.plot(np.nanmean(monc_data[m]['T_mean'][id,:],0),monc_data[m][zvar], color = lcolsmonc[m], linewidth = 3, label = mlabel[m], zorder = 1)
        if pt == 1:
            plt.legend(bbox_to_anchor=(1.5, 1.05), loc=4, ncol=4)


        plt.xlabel('Temperture [K]')
        plt.ylabel('Z [km]')
        plt.xlim([260,271])
        # plt.yticks(np.arange(0,5.01e3,0.5e3))
        # ax1.set_yticklabels([0,' ',1,' ',2,' ',3,' ',4,' ',5])
        plt.ylim(ylims)
        plt.yticks(yticks)
        ax1.yaxis.set_minor_locator(ticker.MultipleLocator(100))
        ax1.set_yticklabels(ytlabels)
        # plt.xlim([0,0.05])
        # plt.xticks(np.arange(0,0.051,0.015))
        #ax1.set_xticklabels([0,' ',0.015,' ',0.03,' ',0.045,' ',0.06])
        # ax1.xaxis.set_minor_locator(ticker.MultipleLocator(0.0075))
    dstr=datenum2date(dates[1])
    # plt.grid('on')
    if pmonc==True:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) + '_' +'_'.join(moutstr) + '_Tprofile'  + '_split.png'
    else:
        fileout = plots_out_dir + dstr.strftime('%Y%m%d') + '_Obs-UMGrid_' + '_'.join(outstr) +'_Tprofile'  + '_split.png'

    plt.savefig(fileout)
    print ('')
    print ('Finished plotting! :)')
    print ('')
    print ('******')


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

def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    machine = 'JASMIN'

    if machine =='JASMIN':
        plots_out_dir='/gws/nopw/j04/ncas_radar_vol1/jutta/PLOTS/CaseStudy/'
        if not os.path.exists(plots_out_dir):
            os.makedirs(plots_out_dir)
        um_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/UM/'
        obs_met_dir=  '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/';
        obs_acas_dir= '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/ACAS/ACAS_AO2018_v2_May2019/';
        obs_rs_dir=   '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/radiosondes/';
        obs_hatpro_dir='/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/HATPRO/';
        obs_albedo_dir='/nfs/a96/MOCCHA/working/data/'
        obs_rad_dir='/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/radiation/'
        obs_dec_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/OBS/HATPRO/'
        inv_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/Inversions/'
        monc_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
        #monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/MONC/output/'
        monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
        ### SET OUTPUT DIRECTORY FOR PLOTS
        plot_out_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/PLOTS/CaseStudy/'

    elif machine =='LEEDS':
    ### INPUT FOLDERS
        um_root_dir = '/nfs/a96/MOCCHA/working/gillian/UM/DATA/'
        obs_met_dir=  '/nfs/a96/MOCCHA/working/jutta/final_data/met_alley/concatenated/';
        obs_acas_dir= '/nfs/a96/MOCCHA/data/ACAS/ACAS_AO2018_v2_May2019/';
        obs_rs_dir=   '/nfs/a96/MOCCHA/working/jutta/final_data/radiosondes/V3/';
        obs_hatpro_dir='/nfs/a96/MOCCHA/working/jutta/final_data/HATPRO/';
        obs_albedo_dir='/nfs/a96/MOCCHA/working/data/'
        obs_rad_dir='/nfs/a96/MOCCHA/working/jutta/requests/Gillian/'
        obs_dec_dir = '/nfs/a96/MOCCHA/working/jutta/requests/Sandeep/'
        monc_root_dir = '/nfs/a96/MOCCHA/working/gillian/MONC_CASES/MOCCHA/output/'
        ### SET OUTPUT DIRECTORY FOR PLOTS
        plot_out_dir = '/nfs/a96/MOCCHA/working/jutta/plots/CaseStudies/ModelComparison/'


    # ### CHOSEN RUN
    # out_dir1 = '25_u-cc568_RA2M_CON/'
    # out_dir2 = '23_u-cc278_RA1M_CASIM/'
    # out_dir3 = '24_u-cc324_RA2T_CON/'

    out_dir = ['23_u-cc278_RA1M_CASIM/',
               '30_u-cg179_RA1M_CASIM/',
              '26_u-cd847_RA1M_CASIM/',
              '27_u-ce112_RA1M_CASIM/']
#    out_dir = [  '30_u-cg179_RA1M_CASIM/' ]
    ### CHOOSE MONC RUNS
    m_out_dir =['22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/']

    # m_out_dir = ['5_control_20180913T0000Z_Wsub-1.5_Fletcher/',
    #             '6_control_20180913T0000Z_Wsub-1.5-1km/',
    #             '7_20180913T0000Z_Wsub-1.5-1km_solAccum-100_inuc-0_iact-3/']



    um_sub_dir = 'OUT_R0/'
    ### CHOOSE DATES TO PLOT
    DATE = 20180913
    strdate=str(DATE)

    sdate = dtime.datetime.strptime('2018091300','%Y%m%d%H')
    edate = dtime.datetime.strptime('2018091313','%Y%m%d%H')
    dates = [date2datenum(sdate),date2datenum(edate)]

    #---- MONC SPIN UP TIME
    monc_spin = 6 *60 *60

    #---- SPLIT PROFILES IN TIME JUNKS
    prof_time=[[dates[0], dates[0]+4/24],
                [dates[0]+4/24, dates[0]+8/24],
                [dates[0]+8/24, dates[0]+14/24]]


    if not os.path.exists(plot_out_dir):
        os.mkdir(plot_out_dir)


    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')

    ### -------------------------------------------------------------------------
    ### define UM input filename
    ### -------------------------------------------------------------------------
    filename_um=[]
    nc={}
    for m in range(0, len(out_dir)):
        filename_um = um_root_dir + out_dir[m] + um_sub_dir + strdate + '_oden_metum.nc'
        nc[m] = Dataset(filename_um,'r')
    # -------------------------------------------------------------
    print ('')
    #### LOAD IN SPECIFIC DIAGNOSTICS
    ### BASE UM RUNS (UM_RA2M/UM_RA2T)

    for var in nc[0].variables: print(var)
    var_list1 = ['u_10m','v_10m', 'air_temperature_at_1.5m','q_1.5m','rh_1.5m','visibility','dew_point_temperature_at_1.5m','LWP','IWP',
                'surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'sensible_heat_flux','latent_heat_flux', 'bl_depth','bl_type','temperature','theta']
                #PLOT FROM CLOUDNET:
                #'temperature','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                #'cloud_fraction','radr_refl','rainfall_flux','snowfall_flux',]#

    um_data = {}
    for m in range(0,len(out_dir)):
        um_data[m]={}
        datenum = date2datenum(dtime.datetime.strptime(strdate,'%Y%m%d'))
        um_data[m]['time'] = datenum + (nc[m].variables['forecast_time'][:]/24.0)
        ### define height arrays explicitly
        um_data[m]['height'] = nc[m].variables['height'][:]
        for var in nc[m].variables: print(var)

        print ('Starting on t=0 RA2M data:')
        for j in range(0,len(var_list1)):
            if np.ndim(nc[m].variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                continue
            elif np.ndim(nc[m].variables[var_list1[j]]) >= 1:
                um_data[m][var_list1[j]] = nc[m].variables[var_list1[j]][:]
        nc[m].close()

    #---- load UM Inversions
    print ('**************************')
    print ('Load UM INVERSION DATA')
    for m in range(0,len(out_dir)):
        filename=glob.glob(inv_dir + out_dir[m][0:-1] + '*20m.mat')
        tmp  = readMatlabStruct(filename[0])
        um_data[m]['inv'] = tmp['dec']
        print (um_data[m]['inv'].keys())
        for var in um_data[m]['inv'].keys():
            um_data[m]['inv'][var]=np.squeeze(um_data[m]['inv'][var])

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

    ### -----------------------------------------------------------------
    ###     READ IN MONC DATA
    print ('Loading MONC data:')
    print ('')
    ###1d variables, 2d variables (time,height), 3d variables (time,x,y), 4d variables(time,x,y,z)
    # monc_var_list =[['time_series_2_60','time_series_20_600' ,'time_series_2_600','z','rho', 'LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean'],
    #                 ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction',
    #                 'vapour_mmr_mean','liquid_mmr_mean','rain_mmr_mean','ice_mmr_mean','snow_mmr_mean',
    #                 'graupel_mmr_mean']]
                #    ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp'],
                #    ['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]
    monc_var_list =[['z', 'zn','LWP_mean','IWP_mean','SWP_mean','TOT_IWP_mean','GWP_mean']]
                #    ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction'],
                #    ['liquid_mmr_mean','ice_mmr_mean','graupel_mmr_mean','snow_mmr_mean']]
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
                    time_var_list=time_var_list+[var]
            full_var_list=monc_var_list
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
                        monc_data[m][var] = np.append(monc_data[m][var],ncm.variables[var][:])
    #loading 3d variables
    for m in range(0, len(m_out_dir)):
        for n in range(0,len(monc_3d_filename)):
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
                    monc_data[m][var] =np.addpend(monc_data[m][var],pyd[var])


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

    #################################################################################################################################
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

    #---- load MONC Inversions
    print ('**************************')
    print ('Load MONC INVERSION DATA')
    for m in range(0,len(m_out_dir)):
        filename=glob.glob(inv_dir + m_out_dir[m][0:-1] + '*20m.mat')
        tmp  = readMatlabStruct(filename[0])
        tmp = tmp['dec']
        print (tmp.keys())
        monc_data[m]['inv']={}
        for var in tmp.keys():
            c, ia, ib = intersect_mtlb(monc_data[m]['time1'],tmp['mday'])
        #    print(ia,ib)
            tmp2=np.argwhere(np.array(tmp[var].shape) == len(np.squeeze(tmp['mday'])))
            if tmp2 ==0:
                if np.squeeze(len(tmp[var]))==1:
                    monc_data[m]['inv'][var]=np.squeeze(tmp[var][ib])
                else:
                    monc_data[m]['inv'][var]=np.squeeze(tmp[var][ib,:])
            elif tmp2 ==1:
                monc_data[m]['inv'][var]=np.squeeze(tmp[var][:,ib])

# -------------------------------------------------------------
# Load observations
# -------------------------------------------------------------
    print ('Loading observations:')
    print('')
            # -------------------------------------------------------------
            # Which file does what?
            # -------------------------------------------------------------
            #### ice station: net LW / net SW
                    #### obs['ice_station_fluxes']/mast_radiation_30min_v2.3.mat
            #### obs['foremast']:
                    #### obs['foremast']/ACAS_AO2018_obs['foremast']_30min_v2_0.nc
            #### 7th deck: temperature, surface temperature, RH, downwelling SW, downwelling LW
                    #### 7thDeck/ACAS_AO2018_WX_30min_v2_0.nc
    obs={}
    print ('**************************')
    print ('Load ice station data from Jutta...')
    filename = 'AO2018_metalley_01min_v3.0.mat'
    obs['metalley'] = readMatlabStruct(obs_met_dir + filename)
    print(obs['metalley'].keys())

    print ('**************************')
    print ('Load decoupling height data from Jutta...')
    filename = '2018091300-2018091314_smc_decoupling_sandeep_Scb_V3.mat'
    obs['dec'] = readMatlabStruct(obs_dec_dir + filename)
    print(obs['dec'].keys())
    for var in obs['dec'].keys():
        obs['dec'][var]=np.squeeze(obs['dec'][var])


    print ('**************************')
    print ('Load HATPRO data used by Cloudnet...')
    filename='HATPRO_LWP_IWV_30s_V3_userready.mat'
    obs['hatpro'] = readMatlabStruct(obs_hatpro_dir + filename)
    print (obs['hatpro'].keys())
    for var in obs['hatpro'].keys():
        obs['hatpro'][var]=np.squeeze(obs['hatpro'][var])

    filename='HATPRO_T_corrected_inversionheights_thetaE_V1.mat'
    tmp= readMatlabStruct(obs_hatpro_dir + filename)
    obs['hatpro_temp'] = tmp['HTdecepot']
    obs['hatpro_temp']['Z'] = tmp['HTcompepot']['Z']
    print (obs['hatpro_temp'].keys())
    for var in obs['hatpro_temp'].keys():
        obs['hatpro_temp'][var]=np.squeeze(obs['hatpro_temp'][var])

    #print ('Load albedo estimates from Michael...')
    #obs['albedo'] = readMatlabStruct(obs_albedo_dir + 'MOCCHA_Albedo_estimates_Michael.mat')
    print ('**************************')
    print ('Load cleaned and pre processed radiation data ..')
    obs['ship_rad']={}
    obs['ice_rad']={}

    nc4 =Dataset(obs_rad_dir + 'MetData_Gillian_V3_30minres.nc','r')

    var_list_srad=(['SWdship', 'LWdship', 'SWnetship', 'LWnetship', 'SWuship'])
    var_list_irad=(['SWuice', 'LWdice', 'LWuice', 'SWdice'])

    obs['ship_rad']['time']=nc4.variables['time2']
    for j in range(0,len(var_list_srad)):
        obs['ship_rad'][var_list_srad[j]] = nc4.variables[var_list_srad[j]][:]

    obs['ice_rad']['time']=nc4.variables['time3']
    for j in range(0,len(var_list_irad)):
        obs['ice_rad'][var_list_irad[j]] = nc4.variables[var_list_irad[j]][:]
    nc4.close

    print('')
    print(obs['ice_rad'].keys())
    print(obs['ship_rad'].keys())
    print('')
    print ('Load radiosonde data from Jutta...')
    obs['sondes'] = readMatlabStruct(obs_rs_dir + '/SondeData_h10int_V03.mat')

    print ('**************************')
    print ('Load RS observations inversion height data from Jutta...')
    obs['inversions'] = readMatlabStruct(obs_rs_dir + '/InversionHeights_RSh05int_final_V03.mat')

    #print ('Load foremast data from John...')
    #obs['foremast'] = Dataset(obs_acas_dir + '/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    #print ('Load 7th deck weather station data from John...')
    #obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print ('**************************')
    print ('Load weather sensor data from John...')
    obs['pwd'] = readMatlabStruct(obs_acas_dir + 'ACAS_AO2018_PWD_30min_v1_0.mat')

    print ('...')

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
            moutstr.append('qin2-thqvTend-noice')
        elif m_out_dir[m][:2] == '23':
            mlabel.append('MONC thForcing-0600-0000')
            moutstr.append('thForcing-0600-0000')
        else:
            label.append('undefined_label')
            moutstr.append('')


    # -------------------------------------------------------------
    # Plot paper figures
    # -------------------------------------------------------------
    #figure = plot_surfaceVariables(obs,plot_out_dir, dates, um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_lwp(obs,plot_out_dir, dates, um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    #figure = plot_BLDepth_SMLDepth(obs,plot_out_dir, dates, um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_T_profiles_split(obs,plots_out_dir,dates, prof_time,um_data=um_data,label=label,outstr=outstr,  monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)
    figure = plot_BLDepth_SMLDepth(obs,plot_out_dir, dates, um_data=um_data,label=label,outstr=outstr, monc_data=monc_data,mlabel=mlabel,moutstr=moutstr)

    embed()
    #figure = plot_radiation(obs,plot_out_dir, dates,plot_out_dir, um_data=um_data,label=label,outsr=outsr, monc_data=monc_data,mlabel=mlabel,moutsr=moutsr)
    # figure = plot_paperFluxes(data1, data2, data3, month_flag, missing_files, out_dir1, out_dir2, out_dir3, obs, doy, label1, label2, label3)
    # figure = plot_paperRadiation(data1, data2, data3, out_dir1, out_dir2, out_dir3,datenum,label1,label2,label3,plot_out_dir)
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
    # save out working data for debugging purposes
    # -------------------------------------------------------------
    np.save('working_data1', data1)
    np.save('working_data2', data2)
    np.save('working_data3', data3)

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
