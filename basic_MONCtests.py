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


def plot_basicTests( monc_data ):


    print ('******')
    print ('')
    print ('Plotting test figures to check MONC output:')
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

    checkpoint = 47 ### checkpoint at 12h (48th timestep)

    plt.plot(monc_data[0]['th_mean'][0,:],monc_data[0]['zn'],label = 'start')
    plt.plot(monc_data[0]['th_mean'][47,:],monc_data[0]['zn'],label = 'checkpoint restart')
    plt.plot(monc_data[0]['th_mean'][-1,:],monc_data[0]['zn'],label = 'end')

    plt.show()


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
    print (monc_data[0].keys())
    print (monc_data[0]['time2'].shape)
    print (monc_data[0]['time2'])

    figure = plot_basicTests( monc_data )

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
