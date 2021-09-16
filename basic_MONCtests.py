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


def plot_basicTests( monc_data, monc_spin, plots_out_dir, moutstr, mlabel, m_out_dir, data  ):


    print ('******')
    print ('')
    print ('Plotting test figures to check MONC output:')
    print ('')


    checkpoint = []

    if monc_spin == 21600.:
        checkpoint2 = 12 # checkpoint restart at 12h
    elif monc_spin == 28800.:
        checkpoint1 = 8 # checkpoint restart at 8h
        checkpoint2 = 14 # checkpoint restart at 14h

    st_id = int(checkpoint1*4) - 1
    st_ts = np.float(checkpoint1)*3600.
    cp_id = int(checkpoint2*4) - 1
    cp_ts = np.float(checkpoint2)*3600.

    SMALL_SIZE = 12
    MED_SIZE = 14
    LARGE_SIZE = 16

    plt.rc('font',size=MED_SIZE)
    plt.rc('axes',titlesize=MED_SIZE)
    plt.rc('axes',labelsize=MED_SIZE)
    plt.rc('xtick',labelsize=MED_SIZE)
    plt.rc('ytick',labelsize=MED_SIZE)
    plt.rc('legend',fontsize=MED_SIZE)

    ### theta
    fig = plt.figure(figsize=(6,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.9, left = 0.15,
            hspace = 0.3, wspace = 0.1)
    plt.plot(monc_data[0]['th_mean'][0,:],monc_data[0]['zn'],label = 'start')
    if np.size(monc_data[0]['th_mean'],0) >= cp_id:
        plt.plot(monc_data[0]['th_mean'][cp_id,:],monc_data[0]['zn'],label = 'checkpoint restart')
        plt.plot(monc_data[0]['th_mean'][int(cp_id)+1,:],monc_data[0]['zn'],label = 'checkpoint restart+1')
    plt.plot(monc_data[0]['th_mean'][-1,:],monc_data[0]['zn'],label = 'end')
    plt.xlabel('$\Theta$ [K]')
    plt.ylabel('Z [m]')
    plt.legend()
    fileout = plots_out_dir + moutstr[0] + '_meanTH_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()

    ### calculate change in theta (K/hr)
    # theta1 = ((monc_data[0]['th_mean'][cp_id,:] - monc_data[0]['th_mean'][0,:]) / 12) * 24
    # theta2 = ((monc_data[0]['th_mean'][-1,:] - monc_data[0]['th_mean'][int(cp_id)+1,:]) / 8) * 24
    # theta1a = ((monc_data[0]['th_mean'][23,:] - monc_data[0]['th_mean'][0,:]) / 6) * 24
    # theta1b = ((monc_data[0]['th_mean'][cp_id,:] - monc_data[0]['th_mean'][23,:]) / 6) * 24
    #
    # fig = plt.figure(figsize=(6,5))
    # plt.subplots_adjust(top = 0.9, bottom = 0.12, right = 0.9, left = 0.15,
    #         hspace = 0.3, wspace = 0.1)
    # plt.plot([0,0],[0,2.5e3],'--',color='lightgrey')
    # plt.plot(theta1a,monc_data[0]['zn'],label = '0-6h')
    # plt.plot(theta1b,monc_data[0]['zn'],label = '6-12h')
    # plt.plot(theta1,monc_data[0]['zn'],label = '0-12h')
    # plt.plot(theta2,monc_data[0]['zn'],label = '12-20h')
    # plt.xlabel('$\Delta \Theta$ [K/day]')
    # plt.ylabel('Z [m]')
    # plt.legend()
    # fileout = plots_out_dir + moutstr[0] + '_DeltaTH_' + mlabel[0] + '.png'
    # plt.savefig(fileout)
    # plt.close()

    ### liquid mass
    fig = plt.figure(figsize=(9,4))
    plt.subplots_adjust(top = 0.93, bottom = 0.14, right = 0.98, left = 0.12,
            hspace = 0.3, wspace = 0.26)
    plt.subplot(121)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_cloud_liquid_mass_mean'])*1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.ylabel('Z [m]')
    plt.title('LWMR [g/kg]')
    plt.colorbar()
    plt.subplot(122)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['ndrop_tot_mean'])/1e6,
        vmin=0., vmax=10.,
        )
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.title('N$_{d}$ [cm$^{-3}$]')
    plt.colorbar()
    fileout = plots_out_dir + moutstr[0] + '_LWC-NdTimeseries_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()

    ### ice mass
    fig = plt.figure(figsize=(8,7))
    plt.subplots_adjust(top = 0.93, bottom = 0.11, right = 0.98, left = 0.12,
            hspace = 0.3, wspace = 0.26)
    plt.subplot(221)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose((monc_data[0]['q_ice_mass_mean']+monc_data[0]['q_snow_mass_mean']+monc_data[0]['q_graupel_mass_mean']))*1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.ylabel('Z [m]')
    plt.title('QISG_MR [g/kg]')
    plt.colorbar()
    plt.subplot(222)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_ice_mass_mean'])*1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.title('QI_MR [g/kg]')
    plt.colorbar()
    plt.subplot(223)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_snow_mass_mean'])*1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.ylabel('Z [m]')
    plt.title('QS_MR [g/kg]')
    plt.colorbar()
    plt.subplot(224)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_graupel_mass_mean'])*1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.title('QG_MR [g/kg]')
    plt.colorbar()
    fileout = plots_out_dir + moutstr[0] + '_IWCTimeseries_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()

    ### ice number
    fig = plt.figure(figsize=(8,7))
    plt.subplots_adjust(top = 0.93, bottom = 0.11, right = 0.98, left = 0.12,
            hspace = 0.3, wspace = 0.26)
    plt.subplot(221)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['nisg_tot_mean'])/1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.ylabel('Z [m]')
    plt.title('N$_{isg}$ [/L]')
    plt.colorbar()
    plt.subplot(222)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_ice_number_mean'])/1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.title('N$_{i}$ [/L]')
    plt.colorbar()
    plt.subplot(223)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_snow_number_mean'])/1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.ylabel('Z [m]')
    plt.title('N$_{s}$ [/L]')
    plt.colorbar()
    plt.subplot(224)
    plt.pcolor(monc_data[0]['time2'],monc_data[0]['zn'],
        np.transpose(monc_data[0]['q_graupel_number_mean'])/1e3)
    plt.plot([cp_ts,cp_ts],[0,2.5e3],'--',color='white',zorder=3)
    plt.plot([monc_spin,monc_spin],[0,2.5e3],'--',color='red',zorder=3)
    plt.xlabel('Time [s]')
    plt.title('N$_{g}$ [/L]')
    plt.colorbar()
    fileout = plots_out_dir + moutstr[0] + '_ICNCTimeseries_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()


    ### u profiles
    fig = plt.figure(figsize=(14,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.98, left = 0.08,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(141)
    plt.plot(data['sonde']['u'], data['sonde']['Z'], 'k--')
    plt.plot(monc_data[0]['u_wind_mean'][0,:],monc_data[0]['zn'],label = 't=0h')
    plt.ylabel('Z [m]')
    plt.xlabel('U [m/s]')
    plt.ylim([0,2.5e3])
    plt.title('t=0h')
    plt.subplot(142)
    plt.plot(data['sonde']['u'], data['sonde']['Z'], 'k--')
    plt.plot(monc_data[0]['u_wind_mean'][st_id,:],monc_data[0]['zn'],label = 't=8h')
    plt.title('t=' + str(checkpoint1) + 'h')
    plt.ylim([0,2.5e3])
    plt.xlabel('U [m/s]')
    plt.subplot(143)
    if np.size(monc_data[0]['u_wind_mean'],0) >= cp_id:
        plt.title('t=' + str(checkpoint2) + 'h')
        plt.plot(data['sonde+1']['u'], data['sonde+1']['Z'], 'k--')
        plt.plot(monc_data[0]['u_wind_mean'][cp_id,:],monc_data[0]['zn'],label = 't=14h')
        plt.ylim([0,2.5e3])
        plt.xlabel('U [m/s]')
    #     plt.plot(monc_data[0]['u_wind_mean'][int(cp_id)+1,:],monc_data[0]['zn'],label = 'checkpoint restart+1')
    plt.subplot(144)
    plt.plot(data['sonde+2']['u'], data['sonde+2']['Z'], 'k--')
    plt.plot(monc_data[0]['u_wind_mean'][-1,:],monc_data[0]['zn'],label = 't=24h')
    plt.title('t=24h')
    plt.ylim([0,2.5e3])
    plt.xlabel('U [m/s]')
    # plt.legend()
    fileout = plots_out_dir + moutstr[0] + '_Uprofiles_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()

    ### v profiles
    fig = plt.figure(figsize=(14,5))
    plt.subplots_adjust(top = 0.9, bottom = 0.1, right = 0.98, left = 0.08,
            hspace = 0.3, wspace = 0.3)
    plt.subplot(141)
    plt.plot(data['sonde']['v'], data['sonde']['Z'], 'k--')
    plt.plot(monc_data[0]['v_wind_mean'][0,:],monc_data[0]['zn'],label = 't=0h')
    plt.ylabel('Z [m]')
    plt.ylim([0,2.5e3])
    plt.title('t=0h')
    plt.xlabel('V [m/s]')
    plt.subplot(142)
    plt.plot(data['sonde']['v'], data['sonde']['Z'], 'k--')
    plt.plot(monc_data[0]['v_wind_mean'][st_id,:],monc_data[0]['zn'],label = 't=8h')
    plt.title('t=' + str(checkpoint1) + 'h')
    plt.ylim([0,2.5e3])
    plt.xlabel('V [m/s]')
    plt.subplot(143)
    if np.size(monc_data[0]['v_wind_mean'],0) >= cp_id:
        plt.title('t=' + str(checkpoint2) + 'h')
        plt.plot(data['sonde+1']['v'], data['sonde+1']['Z'], 'k--')
        plt.plot(monc_data[0]['v_wind_mean'][cp_id,:],monc_data[0]['zn'],label = 't=14h')
        plt.ylim([0,2.5e3])
        plt.xlabel('V [m/s]')
    #     plt.plot(monc_data[0]['u_wind_mean'][int(cp_id)+1,:],monc_data[0]['zn'],label = 'checkpoint restart+1')
    plt.subplot(144)
    plt.plot(data['sonde+2']['v'], data['sonde+2']['Z'], 'k--')
    plt.plot(monc_data[0]['v_wind_mean'][-1,:],monc_data[0]['zn'],label = 't=24h')
    plt.title('t=24h')
    plt.ylim([0,2.5e3])
    plt.xlabel('V [m/s]')
    # plt.legend()
    fileout = plots_out_dir + moutstr[0] + '_Vprofiles_' + mlabel[0] + '.png'
    plt.savefig(fileout)
    plt.close()

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

    #---- SPLIT PROFILES IN TIME JUNKS
    prof_times=[[dates[0], dates[0]+4/24],
                [dates[0]+4/24, dates[0]+8/24],
                [dates[0]+8/24, dates[0]+14/24]]

    #which server
    machine = 'JASMIN'
    ### Choose observations vertical gridding used in Cloudnet processing (UM/IFS/RADAR)
    obs_switch = 'UM' #RADAR

    if machine =='JASMIN':
        plots_out_dir='/gws/nopw/j04/ncas_radar_vol1/gillian/PLOTS/CaseStudy/'
        if not os.path.exists(plots_out_dir):
            os.makedirs(plots_out_dir)

        um_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/UM/DATA/'
        cn_obs_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/DATA/Cloudnet/Observations/measurements_V7/'
        cn_um_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/Cloudnet/UM/'
        monc_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
        #monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/jutta/MONC/output/'
        monc_avg_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/MONC/output/'
        obs_root_dir = '/gws/nopw/j04/ncas_radar_vol1/gillian/Obs/'

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
               # '23_20180913T0000Z_6hSpin-up_12h0600-0000thTend_20h1200-0600thTend/'
               # '24_20180913T0000Z_12h0600-0000_20h1200-0600thTend_0.1Cooper/'
               # '25_20180913T0000Z_20h0600-0000thTend/'
               # '26A_20180913T0000Z_6hSpinUp_12h0600-0000thTend_20h1200-0600thTend_6-20h0.1Cooper/',
               # '26B_20180913T0000Z_6hSpinUp_12h0600-0000thTend_20h1200-0600thTend_6-20h0-Cooper/'
               # '27A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24hCooper/',
               # '27B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.5Cooper/',
               # '27C_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper/',
               # '27D_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd25/',
               # '27E_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd10/',
               # '27E_CASIMvn0.3.4-MONCr8166-test/',
               # '27F_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd5/',
               # '27G_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_FixedNd10_5KDecouple14-24h/',
               # '28A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-100-ARG/',
               # '28B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-100-Twomey/',
               # '29A_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-20-ARG/',
               # '29B_20180913T0000Z_8hSpinUp_14h0600-0000thTend_24h1200-0600thTend_8-24h0.1Cooper_AccumSolAero-CASIM-20-allAct/',
               '30A_20180913T0000Z_8hSpinUp_8-14hUVRelax0600_14-24hUVRelax1200_8-24h0.1Cooper_FixedNd10/',
               ]
            #'4_control_20180913T0000Z_Wsub-1.5/',
    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    mlabel=[]
    moutstr=[]
    for m in range(0, len(m_out_dir)):
        if m_out_dir[m][:2] == '3_':
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
            mlabel.append('MONC thForcing-12h0600-0000-20h1200-0600')
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
            if m_out_dir[m][-5:] == 'test/':
                mlabel.append('MONC_0.1Cooper_FixedNd10_vnTest')
                moutstr.append('MONC-27C-test')
            else:
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
        elif m_out_dir[m][:3] == '27G':
            mlabel.append('MONC_0.1Cooper_FixedNd10_5KDecouple')
            moutstr.append('MONC-27G')
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
            mlabel.append('MONC_0.1Cooper_FixedNd10_uvRelax')
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
        fname=glob.glob(monc_root_dir + m_out_dir[m] +'*_dg_*.nc')
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
                    ['u_wind_mean'],
                    ['theta_mean','total_cloud_fraction', 'liquid_cloud_fraction','ice_cloud_fraction'],
                    ['liquid_mmr_mean','ice_mmr_mean','graupel_mmr_mean','snow_mmr_mean'],
                #    ['vwp','lwp','rwp','iwp','swp','gwp','tot_iwp'],
                #    ['q_vapour','q_cloud_liquid_mass','q_rain_mass','q_ice_mass','q_snow_mass','q_graupel_mass']]
                    ]


    monc_var_3d_list =['T_mean','p_mean','th_mean','rho_mean','q_vapour_mean',
                        'q_cloud_liquid_mass_mean','q_ice_mass_mean','q_snow_mass_mean','q_graupel_mass_mean',
                        'q_cloud_liquid_number_mean','q_ice_number_mean','q_snow_number_mean','q_graupel_number_mean',
                        'twc_tot_mean','iwc_tot_mean','lwc_tot_mean','nisg_tot_mean','ndrop_tot_mean',
                        #'z', 'zn', 'u_mean', 'v_mean', 'w_mean', 'q_vapour_mean',
                        # 'time1', 'time2', 'p_mean', 'T_mean', 'th_mean', 'rho_mean',
                        # 'q_cloud_liquid_mass_mean', 'q_ice_mass_mean', 'q_snow_mass_mean',
                        # 'q_graupel_mass_mean', 'iwc_tot_mean', 'lwc_tot_mean', 'twc_tot_mean', 'zvar', 'tvar',
                        ]
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

    #################################################################
    ###     READ IN OBS DATA
    #################################################################
    sondes = readMatlabStruct(obs_root_dir + 'radiosondes/SondeData_h10int_V02.mat')

    print ('')
    print (sondes.keys())

    ## -------------------------------------------------------------
    ## Choose sonde for initialisation:
    ## -------------------------------------------------------------
    data = {}
    data['sonde_option'] = '20180913T0000' # '20180912T1800' #'20180913T0000'#

    if data['sonde_option'] == '20180912T1800':
        numindex = 0
    elif data['sonde_option'] == '20180913T0000':
        numindex = 1
    elif data['sonde_option'] == '20180913T0600':
        numindex = 2

    ## -------------------------------------------------------------
    ## Load radiosonde (relative to 20180912 1200UTC
    ## -------------------------------------------------------------
    index256 = np.where(np.logical_or(np.round(sondes['doy'][:,:]) == 256., np.round(sondes['doy'][:,:]) == 257.))

    print (sondes['doy'][:,index256[1][numindex]])
    data['sonde'] = {}
    for k in sondes.keys():
        if k == 'Z': continue
        data['sonde'][k] = sondes[k][:,index256[1][numindex]]
    data['sonde']['Z'] = sondes['Z']
    data['sonde']['u'][data['sonde']['u'] > 1e3] = np.nan
    data['sonde']['v'][data['sonde']['v'] > 1e3] = np.nan

    ### load subsequent sondes
    for i in np.arange(0,3):
        data['sonde+' + str(i+1)] = {}
        print ('sonde+' + str(i+1))
        print (sondes['doy'][:,index256[1][i+1+numindex]])
        for k in sondes.keys():
            if k == 'Z': continue
            data['sonde+' + str(i+1)][k] = sondes[k][:,index256[1][i+1+numindex]]
        data['sonde+' + str(i+1)]['Z'] = sondes['Z']
        data['sonde+' + str(i+1)]['u'][data['sonde+' + str(i+1)]['u'] > 1e3] = np.nan
        data['sonde+' + str(i+1)]['v'][data['sonde+' + str(i+1)]['v'] > 1e3] = np.nan

    ##################################################################################################################################
    ## -------------------------------------------------------------
    ## remove spin up time from monc data
    ## -------------------------------------------------------------
    # monc_data=removeSpinUp(monc_data,monc_spin)

    ## -------------------------------------------------------------
    ## convert monc time to datenum
    ## -------------------------------------------------------------
    # for m in range(0,len(monc_data)):
    #     monc_data[m]['time2']=dates[0] + monc_data[m]['time2']/60/60/24
    #     monc_data[m]['time1']=dates[0] + monc_data[m]['time1']/60/60/24
    #     if 'time3' in monc_data:
    #         monc_data[m]['time3']=dates[0] + monc_data[m]['time3']/60/60/24

###################################################################################################################
###################################################################################################################
################################################ FIGURES ##########################################################
###################################################################################################################
###################################################################################################################

    # -------------------------------------------------------------
    # Plot some basic data to check monc run worked successfully
    # -------------------------------------------------------------
    # print (monc_data.keys())
    # print (monc_data[0].keys())
    # print (monc_data[0]['time2'].shape)
    # print (monc_data[0]['time2'])

    figure = plot_basicTests( monc_data, monc_spin, plots_out_dir, moutstr, mlabel, m_out_dir, data )

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
