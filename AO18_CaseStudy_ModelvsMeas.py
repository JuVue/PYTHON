### SCRIPT TO COMPARE MODEL TO MEASURMENTS
from  __future__ import print_function
import time
import datetime
import numpy as np
#import pandas as pd
from netCDF4 import Dataset
#import diags_MOCCHA as diags
#import diags_varnames as varnames
#import cartopy.crs as ccrs
#import iris
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
#import matplotlib.cm as mpl_cm
import os
#import seaborn as sns

### import python functions
import sys
sys.path.insert(1, '/nfs/a26/lecimb/jutta/GITHUB/PYTHON/py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
#from physFuncts import calcThetaE, calcThetaVL
#from pyFixes import py3_FixNPLoad


def plot_surfaceVariables(data1, data2, data3, obs, out_dir1, out_dir2, out_dir3, datenum, label1, label2, label3,plot_out_dir):

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


    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    #from IPython import embed; embed()
    fig = plt.figure(figsize=(16,12 ))
    #ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    plt.subplot(2,1,1)
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(data1['time'], data1['air_temperature_at_1.5m']-273.15, color = 'darkblue', label = label1)
    plt.plot(data3['time'], data3['air_temperature_at_1.5m']-273.15, color = 'steelblue', label = label3[:-4])
    plt.plot(data2['time'], data2['air_temperature_at_1.5m']-273.15, color = 'mediumseagreen', label = label2[:-7])
    plt.plot(obs['metalley']['mday'], obs['metalley']['t'], color = 'black', label = 'Obs')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    plt.ylabel('T [$^\circ$C]')
    plt.legend(bbox_to_anchor=(-0.08, 0.77, 1., .102), loc=4, ncol=4)
    ax.set_xlim([datenum, datenum+1])
    plt.grid()
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.ylim([-10,0])

    plt.subplot(2,1,2)
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(data1['time'], data1['rh_1.5m'], color = 'darkblue', label = label1)
    plt.plot(data3['time'], data3['rh_1.5m'], color = 'steelblue', label = label3[:-4])
    plt.plot(data2['time'], data2['rh_1.5m'], color = 'mediumseagreen', label = label2[:-7])
    plt.plot(obs['metalley']['mday'], obs['metalley']['rh'], color = 'black', label = 'Obs')#plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    plt.ylabel('RH [%]')
    #plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim([datenum, datenum+1])
    plt.grid()
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H%M'))
    plt.ylim([80,110])

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')

    date=datenum2date(datenum)
    from IPython import embed; embed()
    fileout = os.path.join(plot_out_dir,date.strftime('%Y%m%d') + '_testplot_line.png')
    plt.savefig(fileout)
    #plt.ion()


def plot_radiation(data1, data2, data3, obs, out_dir1, out_dir2, out_dir3, datenum, label1, label2, label3,plot_out_dir):
    print ('******')
    print ('')
    print ('Plotting  timeseries  of radiation terms:')
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
    from IPython import embed; embed()

    ### for reference in figures
    zeros = np.zeros(len(data2['time']))

    # #### add override for data2 to allow 24h data to be used for testing purposes
    # if out_dir2[-4:] == '24h/':
    #     data2['surface_net_LW_radiation'][data2['surface_net_LW_radiation'] == 0] = np.nan
    #     data2['surface_net_SW_radiation'][data2['surface_net_SW_radiation'] == 0] = np.nan

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(12,10))

    ax  = fig.add_axes([0.07,0.7,0.55,0.22])   # left, bottom, width, height
    obs['ice_rad']['netLW'] = obs['ice_rad']['LWdice'][:] - obs['obs_temp'].variables['LWuice'][:]
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




def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### INPUT FOLDERS
    um_root_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/data/calibrated/metum/'
    obs_met_dir=  '/nfs/a96/MOCCHA/working/jutta/final_data/met_alley/concatenated/';
    obs_acas_dir= '/nfs/a96/MOCCHA/data/ACAS/ACAS_AO2018_v2_May2019/';
    obs_rs_dir=   '/nfs/a96/MOCCHA/working/jutta/final_data/radiosondes/V3/';
    obs_hatpro_dir='/nfs/a96/MOCCHA/working/jutta/final_data/HATPRO/';
    obs_albedo_dir='/nfs/a96/MOCCHA/working/data/'
    obs_rad_dir='/nfs/a96/MOCCHA/working/jutta/requests/Gillian/'
    ### CHOSEN RUN
    out_dir1 = '25_u-cc568_RA2M_CON/'
    out_dir2 = '23_u-cc278_RA1M_CASIM/'
    out_dir3 = '24_u-cc324_RA2T_CON/'

    ### CHOOSE DATES TO PLOT
    DATES = 20180913

    ### SET OUTPUT DIRECTORY FOR PLOTS
    plot_out_dir = '/nfs/a96/MOCCHA/working/jutta/plots/CaseStudies/ModelComparison/'

    if not os.path.exists(plot_out_dir):
        os.mkdir(plot_out_dir)


    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')

    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------

    strdate = str(DATES)
    datenum = date2datenum(datetime.datetime.strptime(strdate,'%Y%m%d'))

    filename_um1 = um_root_dir + out_dir1 + strdate + '_oden_metum.nc'
    filename_um2 = um_root_dir + out_dir2 + strdate + '_oden_metum.nc'
    filename_um3 = um_root_dir + out_dir3 + strdate + '_oden_metum.nc'

    print (filename_um1)
    print (filename_um2)
    print (filename_um3)
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
    # -------------------------------------------------------------
    print ('')
    ##  for var in nc.variables: print (var)
    #### LOAD IN SPECIFIC DIAGNOSTICS
    ### BASE UM RUNS (UM_RA2M/UM_RA2T)
    var_list1 = ['u_10m','v_10m', 'air_temperature_at_1.5m','q_1.5m','rh_1.5m','visibility','dew_point_temperature_at_1.5m','LWP','IWP',
                'surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'sensible_heat_flux','latent_heat_flux']
                #PLOT FROM CLOUDNET:
                #'temperature','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                #'cloud_fraction','radr_refl','rainfall_flux','snowfall_flux',]#

    #### CASIM var_list1 only different for  qnice and qndrop,qice / qsnow
    var_list2 = var_list1

    data1 = {}
    data2 = {}
    data3 = {}
    data1['time'] = datenum + (nc1.variables['forecast_time'][:]/24.0)
    data2['time'] = datenum + (nc2.variables['forecast_time'][:]/24.0)
    data3['time'] = datenum + (nc3.variables['forecast_time'][:]/24.0)

    ### define height arrays explicitly
    data1['height'] = nc1.variables['height'][:]
    data2['height'] = nc2.variables['height'][:]
    data3['height'] = nc3.variables['height'][:]

    ## um1
    print ('Starting on t=0 RA2M data:')
    for j in range(0,len(var_list1)):
        if np.ndim(nc1.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
            continue
        elif np.ndim(nc1.variables[var_list1[j]]) >= 1:
            data1[var_list1[j]] = nc1.variables[var_list1[j]][:]
    nc1.close()

    #### um2
    print ('Starting on t=0 CASIM data:')
    for j in range(0,len(var_list2)):
        if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
            continue
        elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
            data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
    nc2.close()

    print ('Starting on t=0 RA2T data:')
    for j in range(0,len(var_list1)):
        if var_list1[j] in nc3.variables:
            if np.ndim(nc3.variables[var_list1[j]]) == 0:     # ignore horizontal_resolution
                continue
            elif np.ndim(nc3.variables[var_list1[j]]) >= 1:
                data3[var_list1[j]] = nc3.variables[var_list1[j]][:]
    nc3.close()
    print('')
    print (var_list2[j])
    print('')
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
    
    print ('Load ice station data from Jutta...')
    filename = 'AO2018_metalley_01min_v3.0.mat'
    obs['metalley'] = readMatlabStruct(obs_met_dir + filename)
    print(obs['metalley'].keys())

    # print ('Load HATPRO data used by Cloudnet...')
    # filename='HATPRO_LWP_IWV_30s_V3_userready.mat'
    # obs['hatpro'] = readMatlabStruct(obs_hatpro_dir + filename)
    # print (obs['hatpro'].keys())
    #
    # obs['hatpro']['iwv'] = np.squeeze(obs['hatpro']['iwv'])
    # obs['hatpro']['mday'] = np.squeeze(obs['hatpro']['mday'])
    # obs['hatpro']['lwp'] = np.squeeze(obs['hatpro']['lwp'])
    # obs['hatpro']['lwpflag'] = np.squeeze(obs['hatpro']['lwp_corflag'])
    # obs['hatpro']['iwvflag'] = np.squeeze(obs['hatpro']['iwv_corflag'])
    # obs['hatpro']['rainflag'] = np.squeeze(obs['hatpro']['rainflag'])
    # obs['hatpro']['doy'] = calcTime_Mat2DOY(obs['hatpro']['mday'])


    #print ('Load albedo estimates from Michael...')
    #obs['albedo'] = readMatlabStruct(obs_albedo_dir + 'MOCCHA_Albedo_estimates_Michael.mat')

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

    print ('Load observations inversion height data from Jutta...')
    obs['inversions'] = readMatlabStruct(obs_rs_dir + '/InversionHeights_RSh05int_final_V03.mat')

    #print ('Load foremast data from John...')
    #obs['foremast'] = Dataset(obs_acas_dir + '/ACAS_AO2018_foremast_30min_v2_0.nc','r')

    #print ('Load 7th deck weather station data from John...')
    #obs['deck7th'] = Dataset(obs_root_dir + '7thDeck/ACAS_AO2018_WX_30min_v2_0.nc','r')

    print ('Load weather sensor data from John...')
    obs['pwd'] = readMatlabStruct(obs_acas_dir + 'ACAS_AO2018_PWD_30min_v1_0.mat')

    print ('...')

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
    # Plot paper figures
    # -------------------------------------------------------------
    #figure = plot_surfaceVariables(data1, data2, data3, obs, out_dir1, out_dir2, out_dir3,datenum,label1,label2,label3,plot_out_dir)
    figure = plot_radiation(data1, data2, data3, obs, out_dir1, out_dir2, out_dir3,datenum,label1,label2,label3,plot_out_dir)
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
