### SCRIPT TO COMPARE MODEL TO MEASURMENTS
from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset
import os

#### import python functions
import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcThetaE, calcThetaVL
from pyFixes import py3_FixNPLoad


def plot_paperRadiation(data1, data2, data3, out_dir1, out_dir2, out_dir3, DATES, label1, label2, label3):

    import iris.plot as iplt
    import iris.quickplot as qplt
    import iris.analysis.cartography
    import cartopy.crs as ccrs
    import cartopy
    import matplotlib.cm as mpl_cm
    from time_functions import calcTime_Mat2DOY


    print ('******')
    print ('')
    print ('Plotting combined timeseries and PDFs of radiation terms:')
    print ('')

    ##################################################
    #### 	CARTOPY
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

    #################################################################
    ## create figure and axes instances
    #################################################################
    ### -------------------------------
    ### Build figure (timeseries)
    ### -------------------------------
    fig = plt.figure(figsize=(18,12 ))

    ax  = fig.add_axes([0.07,0.7,0.53,0.22])   # left, bottom, width, height
    ax = plt.gca()
    yB = [-10, 120]
    plt.plot(data1['time'], data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    plt.plot(data3['time'], data3['surface_net_SW_radiation'].data, color = 'steelblue', label = label3[:-4])
    plt.plot(data2['time'], data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    plt.ylabel('SW$_{net}$ [W m$^{-2}$]')
    plt.legend(bbox_to_anchor=(-0.08, 0.67, 1., .102), loc=4, ncol=3)
    ax.set_xlim([datenum datenum+1])

    #plt.xticks([230,235,240,245,250,255])
    #ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    plt.ylim([-3,120])
    #
    # ax  = fig.add_axes([0.07,0.4,0.53,0.22])   # left, bottom, width, height
    # ax = plt.gca()
    # yC = [-90, 10]
    # plt.plot([240.0,240.0],[yC[0],yC[-1]],'--', color='grey')
    # plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # plt.plot(data1['time'], data1['surface_net_LW_radiation'].data, color = 'darkblue')
    # plt.plot(data4['time'], data3['surface_net_LW_radiation'].data, color = 'steelblue')
    # plt.plot(data2['time'], data2['surface_net_LW_radiation'].data, color = 'mediumseagreen')
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_net_lw'].data, color = 'gold')
    # else:
    #     plt.plot(data3['time'], data3['surface_net_LW_radiation'].data, color = 'gold')
    # plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['LWnet_ice'], color = 'grey', linewidth = 3)
    # plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['LWnet_ship'], color = 'k')
    # plt.ylabel('LW$_{net}$ [W m$^{-2}$]')
    # ax.set_xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    # plt.ylim([-90,5])
    #
    # ax  = fig.add_axes([0.07,0.1,0.53,0.22])   # left, bottom, width, height
    # ax = plt.gca()
    # yA = [-65, 85]
    # #### plot periods
    # ##  melt/freeze threshold
    # plt.plot([240.0,240.0],[yA[0],yA[-1]],'--', color='grey')
    # ## p3
    # plt.plot([data1['time_hrly'][p3[0][0]], data1['time_hrly'][p3[0][-1]]], [-50, -50], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p3[0][0]], -50, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p3[0][-1]], -50, '.', color = 'r')
    # plt.annotate('P3', xy=(227.5,-45), xytext=(227.5,-45), fontsize = 12, color = 'r')
    # plt.plot([data1['time_hrly'][p3[0][-1]], data1['time_hrly'][p3[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    # ## p4
    # plt.plot([data1['time_hrly'][p4[0][0]], data1['time_hrly'][p4[0][-1]]], [-50, -50], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p4[0][0]], -50, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p4[0][-1]], -50, '.', color = 'r')
    # plt.annotate('P4', xy=(234.5,-45), xytext=(234.5,-45), fontsize = 12, color = 'r')
    # # plt.plot([data1['time_6hrly'][p4[0][-1]], data1['time_6hrly'][p4[0][-1]]], [yA[0],yA[-1]], '--', color = 'r', linewidth = 1)
    # ## p5
    # plt.plot([data1['time_hrly'][p5[0][0]], data1['time_hrly'][p5[0][-1]]], [70, 70], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p5[0][0]], 70, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p5[0][-1]], 70, '.', color = 'r')
    # plt.annotate('P5', xy=(243,59), xytext=(243,59), fontsize = 12, color = 'r')
    # plt.plot([data1['time_hrly'][p5[0][-1]], data1['time_hrly'][p5[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    # ## p6
    # plt.plot([data1['time_hrly'][p6[0][0]], data1['time_hrly'][p6[0][-1]]], [70, 70], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p6[0][0]], 70, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p6[0][-1]], 70, '.', color = 'r')
    # plt.annotate('P6', xy=(248.75,59), xytext=(248.75,59), fontsize = 12, color = 'r')
    # plt.plot([data1['time_hrly'][p6[0][-1]], data1['time_hrly'][p6[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    # ## p7
    # plt.plot([data1['time_hrly'][p7[0][0]], data1['time_hrly'][p7[0][-1]]], [70, 70], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p7[0][0]], 70, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p7[0][-1]], 70, '.', color = 'r')
    # plt.annotate('P7', xy=(253,59), xytext=(253,59), fontsize = 12, color = 'r')
    # plt.plot([data1['time_hrly'][p7[0][-1]], data1['time_hrly'][p7[0][-1]]], [yA[0],yA[-1]], '-.', color = 'r', linewidth = 1)
    # ## p8
    # plt.plot([data1['time_hrly'][p8[0][0]], data1['time_hrly'][p8[0][-1]]], [70, 70], '-.', color = 'r')
    # plt.plot(data1['time_hrly'][p8[0][0]], 70, '.', color = 'r')
    # plt.plot(data1['time_hrly'][p8[0][-1]], 70, '.', color = 'r')
    # plt.annotate('P8', xy=(256.5,59), xytext=(256.5,59), fontsize = 12, color = 'r')
    # ##
    # plt.plot(data2['time'], zeros,'--', color='lightgrey')
    # plt.plot(data1['time'], data1['surface_net_LW_radiation'].data + data1['surface_net_SW_radiation'].data, color = 'darkblue', label = label1)
    # plt.plot(data4['time'], data4['surface_net_LW_radiation'].data + data4['surface_net_SW_radiation'].data, color = 'steelblue', label = label4[:-4])
    # plt.plot(data2['time'], data2['surface_net_LW_radiation'].data + data2['surface_net_SW_radiation'].data, color = 'mediumseagreen', label = label2)
    # if ifs_flag == True:
    #     plt.plot(data3['time'], data3['sfc_net_lw'].data + data3['sfc_net_sw'].data, color = 'gold', label = label3)
    # else:
    #     plt.plot(data3['time'], data3['surface_net_LW_radiation'].data + data3['surface_net_SW_radiation'].data, color = 'gold', label = label3)
    # plt.plot(obs['fixed_radiation']['time_ice'], obs['fixed_radiation']['LWnet_ice'] + obs['fixed_radiation']['SWnet_ice'], color = 'grey',  linewidth = 3, label = 'Ice_station')
    # plt.plot(obs['fixed_radiation']['time_ship'], obs['fixed_radiation']['LWnet_ship'] + obs['fixed_radiation']['SWnet_ship'], color = 'k', label = 'Ship')
    # plt.ylabel('Net Radiation [W m$^{-2}$]')
    # ax.set_xlim([doy[0],doy[-1]])
    # plt.xticks([230,235,240,245,250,255])
    # ax.set_xticklabels(['18 Aug','23 Aug','28 Aug','2 Sep','7 Sep','12 Sep'])
    # plt.ylim([-60,80])
    # plt.xlabel('Date')
    #
    # ### -------------------------------
    # ### Build figure (PDFs)
    # ### -------------------------------
    # # f, axes = plt.subplots(2, 1, figsize=(7, 7))#, sharex=True)
    # # fig = plt.figure(figsize=(7,9))
    # # plt.subplots_adjust(top = 0.95, bottom = 0.1, right = 0.95, left = 0.1,
    # #         hspace = 0.3, wspace = 0.15)
    # # plt.subplot(211)
    #
    # #### only compare with observations where we have data:
    # ####        all model data share a timestamp
    # melt = np.where(np.logical_and(data1['time_hrly'][:-3] >= obs['fixed_radiation']['time_ice'][0], data1['time_hrly'][:-3] < 240.0))
    # freeze = np.where(data1['time_hrly'][:-3] >= 240.0)
    #
    # icemelt = np.where(obs['fixed_radiation']['time_ice'] < 240.0)
    # icefreeze = np.where(obs['fixed_radiation']['time_ice'] >= 240.0)
    # shipmelt = np.where(obs['fixed_radiation']['time_ship'] < 240.0)
    # shipfreeze = np.where(obs['fixed_radiation']['time_ship'] >= 240.0)
    #
    # # sw1 = data1['surface_net_SW_radiation'][data1['hrly_flag']]
    # # lw1 = data1['surface_net_LW_radiation'][data1['hrly_flag']]
    # # sw2 = data2['surface_net_SW_radiation'][data2['hrly_flag']]
    # # lw2 = data2['surface_net_LW_radiation'][data2['hrly_flag']]
    # # sw3 = data3['sfc_net_sw'][data3['hrly_flag']]
    # # lw3 = data3['sfc_net_lw'][data3['hrly_flag']]
    # # sw4 = data4['surface_net_SW_radiation'][data4['hrly_flag']]
    # # lw4 = data4['surface_net_LW_radiation'][data4['hrly_flag']]
    # sw1 = data1['fixed_radiation']['SWnet']
    # lw1 = data1['fixed_radiation']['LWnet']
    # sw2 = data2['fixed_radiation']['SWnet']
    # lw2 = data2['fixed_radiation']['LWnet']
    # sw3 = data3['fixed_radiation']['SWnet']
    # lw3 = data3['fixed_radiation']['LWnet']
    # sw4 = data4['fixed_radiation']['SWnet']
    # lw4 = data4['fixed_radiation']['LWnet']
    #
    # ax  = fig.add_axes([0.64,0.7,0.15,0.22])   # left, bottom, width, height
    # yEmax = 0.1
    # plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    # f = sns.distplot(sw1[melt], hist=False, color="darkblue", kde_kws={"shade": True})
    # f = sns.distplot(sw4[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    # f = sns.distplot(sw2[melt], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # f = sns.distplot(sw3[melt], hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipmelt], hist=False, color="black")
    # plt.annotate('Melt', xy=(87,0.087), xytext=(87,0.087), fontsize = 14)
    # plt.xlim([-10,110])
    # plt.ylim([0,yEmax])
    # plt.xlabel('SW$_{net}$ [W m$^{-2}$]')
    # # ## Obs_ship peak SWnet = 18.22
    # # ## Obs_ice peak SWnet = 13.30
    # # ## ECMWF_IFS peak SWnet = 21.61
    # # ## UM_CASIM-100 peak SWnet = 26.92
    # # ## UM_RA2T peak SWnet = 41.85
    # # ## UM_RA2M peak SWnet = 40.41
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('UM_RA2M peak SWnet = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')
    #
    # ax  = fig.add_axes([0.64,0.4,0.15,0.22])   # left, bottom, width, height
    # yFmax = 0.16
    # plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    # f = sns.distplot(lw1[melt], hist=False, color="darkblue", kde_kws={"shade": True})
    # f = sns.distplot(lw4[melt], hist=False, color="steelblue", kde_kws={"shade": True})
    # f = sns.distplot(lw2[melt], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # f = sns.distplot(lw3[melt], hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['LWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['LWnet_ship'][shipmelt], hist=False, color="black")
    # plt.annotate('Melt', xy=(0,0.14), xytext=(0,0.14), fontsize = 14)
    # plt.xlim([-80,20])
    # plt.ylim([0,yFmax])
    # plt.xlabel('LW$_{net}$ [W m$^{-2}$]')
    # # # UM_RA2M peak LWnet = -3.40
    # # # UM_RA2T peak LWnet = -4.01
    # # # UM_CASIM-100 peak LWnet = -2.05
    # # # ECMWF_IFS peak LWnet = -3.68
    # # # Obs_ice peak LWnet = -3.71
    # # # Obs_ship peak LWnet = -2.12
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('Obs_ship peak LWnet = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')
    #
    # ax  = fig.add_axes([0.64,0.1,0.15,0.22])   # left, bottom, width, height
    # yDmax = 0.08
    # plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    # crf1 = sw1[melt] + lw1[melt]
    # f = sns.distplot(crf1, hist=False, color="darkblue", kde_kws={"shade": True})
    # crf4 = sw4[melt] + lw4[melt]
    # f = sns.distplot(crf4, hist=False, color="steelblue", kde_kws={"shade": True})
    # crf2 = sw2[melt] + lw2[melt]
    # f = sns.distplot(crf2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # crf3 = sw3[melt] + lw3[melt]
    # f = sns.distplot(crf3, hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icemelt] + obs['fixed_radiation']['LWnet_ice'][icemelt], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipmelt] + obs['fixed_radiation']['LWnet_ship'][shipmelt], hist=False, color="black")
    # plt.annotate('Melt', xy=(47,0.07), xytext=(47,0.07), fontsize = 14)
    # plt.xlabel('Net Radiation [W m$^{-2}$]')
    # plt.xlim([-80,80])
    # plt.ylim([0,yDmax])
    # # # UM_RA2M peak netRad = 33.57
    # # # UM_RA2T peak netRad = 37.46
    # # # UM_CASIM-100 peak netRad = 22.47
    # # # ECMWF_IFS peak netRad = 15.56
    # # # Obs_ice peak netRad = 8.31
    # # # Obs_ship peak netRad = 16.94
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('UM_RA2T peak netRad = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')
    #
    # ax  = fig.add_axes([0.83,0.7,0.15,0.22])   # left, bottom, width, height
    # yEmax = 0.1
    # plt.plot([0,0],[0,yEmax],'--', color='lightgrey')
    # f = sns.distplot(sw1[freeze], hist=False, color="darkblue", kde_kws={"shade": True})
    # f = sns.distplot(sw4[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    # f = sns.distplot(sw2[freeze], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # f = sns.distplot(sw3[freeze], hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipfreeze], hist=False, color="black")
    # plt.annotate('Freeze', xy=(77,0.087), xytext=(77,0.087), fontsize = 14)
    # plt.xlim([-10,110])
    # plt.ylim([0,yEmax])
    # plt.xlabel('SW$_{net}$ [W m$^{-2}$]')
    # # # UM_RA2M peak SWnet = 25.25
    # # # UM_RA2T peak SWnet = 26.57
    # # # UM_CASIM-100 peak SWnet = 10.65
    # # # ECMWF_IFS peak SWnet = 10.00
    # # # Obs_ice peak SWnet = 8.67
    # # # Obs_ship peak SWnet = 7.87
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('UM_RA2T peak SWnet = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')
    #
    # ax  = fig.add_axes([0.83,0.4,0.15,0.22])   # left, bottom, width, height
    # yFmax = 0.16
    # plt.plot([0,0],[0,yFmax],'--', color='lightgrey')
    # f = sns.distplot(lw1[freeze], hist=False, color="darkblue", kde_kws={"shade": True})
    # f = sns.distplot(lw4[freeze], hist=False, color="steelblue", kde_kws={"shade": True})
    # f = sns.distplot(lw2[freeze], hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # f = sns.distplot(lw3[freeze], hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['LWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['LWnet_ship'][shipfreeze], hist=False, color="black")
    # plt.annotate('Freeze', xy=(-8,0.14), xytext=(-8,0.14), fontsize = 14)
    # plt.xlim([-80,20])
    # plt.ylim([0,yFmax])
    # plt.xlabel('LW$_{net}$ [W m$^{-2}$]')
    # # # UM_RA2M peak LWnet = -9.31
    # # # UM_RA2T peak LWnet = -11.49
    # # # UM_CASIM-100 peak LWnet = -5.65
    # # # ECMWF_IFS peak LWnet = -6.83
    # # # Obs_ice peak LWnet = -7.57
    # # # Obs_ship peak LWnet = -6.50
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('UM_RA2T peak LWnet = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')
    #
    # # plt.subplot(212)
    # ax  = fig.add_axes([0.83,0.1,0.15,0.22])   # left, bottom, width, height
    # yDmax = 0.08
    # plt.plot([0,0],[0,yDmax],'--', color='lightgrey')
    # crf1 = sw1[freeze] + lw1[freeze]
    # f = sns.distplot(crf1, hist=False, color="darkblue", kde_kws={"shade": True})
    # crf4 = sw4[freeze] + lw4[freeze]
    # f = sns.distplot(crf4, hist=False, color="steelblue", kde_kws={"shade": True})
    # crf2 = sw2[freeze] + lw2[freeze]
    # f = sns.distplot(crf2, hist=False, color="mediumseagreen", kde_kws={"shade": True})
    # crf3 = sw3[freeze] + lw3[freeze]
    # f = sns.distplot(crf3, hist=False, color="gold", kde_kws={"shade": True})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ice'][icefreeze] + obs['fixed_radiation']['LWnet_ice'][icefreeze], hist=False, color="grey", kde_kws={"linewidth": 3})
    # f = sns.distplot(obs['fixed_radiation']['SWnet_ship'][shipfreeze] + obs['fixed_radiation']['LWnet_ship'][shipfreeze], hist=False, color="black")
    # plt.annotate('Freeze', xy=(35,0.07), xytext=(35,0.07), fontsize = 14)
    # plt.xlim([-80,80])
    # plt.ylim([0,yDmax])
    # plt.xlabel('Net Radiation [W m$^{-2}$]')
    # # # UM_RA2M peak netRad = 14.91
    # # # UM_RA2T peak netRad = 18.18
    # # # UM_CASIM-100 peak netRad = 6.10
    # # # ECMWF_IFS peak netRad = 5.62
    # # # Obs_ice peak netRad = 0.98
    # # # Obs_ship peak netRad = 1.73
    # # y = f.lines[0].get_ydata()
    # # x = f.lines[0].get_xdata()
    # # maxid = np.argmax(y)
    # # print ('UM_RA2T peak netRad = ' + ('%.2f' % x[maxid]))
    # # plt.plot(x[maxid],y[maxid],'o')

    print ('******')
    print ('')
    print ('Finished plotting! :)')
    print ('')


    fileout = os.path.join(plot_out_dir, strdate '_netSW_netLW_netRad_line.png')
    # plt.savefig(fileout)
    plt.show()





def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### INPUT FOLDERS
    um_root_dir = '/nfs/a96/MOCCHA/working/jutta/CloudNet/Oden/data/calibrated/metum/'
    obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        #ship_filename = '/nfs/a96/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'

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
    datenum = date2datenum(datetime.datetime.fromtimestamp(DATES))

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
    #### LOAD IN SPECIFIC DIAGNOSTICS
    ### BASE UM RUNS (UM_RA2M/UM_RA2T)
    if out_dir1[:21] == '25_u-cc568_RA2M_CON':
        var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','surface_downwelling_LW_radiation','surface_downwelling_SW_radiation',
                'sensible_heat_flux','latent_heat_flux','rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                'cloud_fraction','radr_refl']#,'temp_1.5m']
    else:
            var_list1 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
                'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','qice','uwind','vwind','wwind',
                'cloud_fraction','radr_refl']#,'temp_1.5m','latent_heat_flux']
            var_list3 = var_list1

    ### CASIM RUNS
    var_list2 = ['temperature','surface_net_SW_radiation','surface_net_LW_radiation','sensible_heat_flux',
            'air_temperature_at_1.5m', 'rainfall_flux','snowfall_flux','q','pressure','bl_depth','bl_type','qliq','uwind','vwind','wwind',
            'cloud_fraction','radr_refl']#, 'tke']#,'latent_heat_flux']#,'tke'] #'qnliq','qnice','mixing_length_for_momentum','qice',
            #, 'latent_heat_flux']



    data1 = {}
    data2 = {}
    data3 = {}
    data1['time'] = datenum + (nc1.variables['forecast_time'][:]/24.0)
    data2['time'] = datenum + (nc2.variables['forecast_time'][:]/24.0)
    data3['time'] = datenum + (nc3.variables['forecast_time'][:]/24.0)
    # else:
    #     time_um1 = float(filename_um1[-16:-14]) + (nc1.variables['forecast_time'][:]/24.0)
    #     time_um2 = float(filename_um2[-16:-14]) + (nc2.variables['forecast_time'][:]/24.0)
    #     if ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['time'][:]/24.0)
    #     if not ifs_flag: time_um3 = float(filename_um3[-16:-14]) + (nc3.variables['forecast_time'][:]/24.0)
    #     time_um4 = float(filename_um4[-16:-14]) + (nc4.variables['forecast_time'][:]/24.0)
    #     time_um5 = float(filename_um5[-16:-14]) + (nc5.variables['forecast_time'][:]/24.0)

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
        print (var_list2[j])
        if np.ndim(nc2.variables[var_list2[j]]) == 0:     # ignore horizontal_resolution
            continue
        elif np.ndim(nc2.variables[var_list2[j]]) >= 1:
            data2[var_list2[j]] = nc2.variables[var_list2[j]][:]
    nc2.close()
    print ('Starting on t=0 RA2T data, initialising:')
    for j in range(0,len(var_list1)):
    # print (var_list4[j])
        if var_list3[j] in nc3.variables:
            if np.ndim(nc3.variables[var_list3[j]]) == 0:     # ignore horizontal_resolution
                continue
            elif np.ndim(nc3.variables[var_list3[j]]) >= 1:
                data3[var_list3[j]] = nc3.variables[var_list3[j]][:]
    nc3.close()

    #################################################################
    ## create labels for figure legends - done here so only needs to be done once!
    #################################################################
    label1 = 'undefined_label'
    if out_dir1[:10] == '25_u-cc568': label1 = 'UM_RA2M'
    if out_dir1[:10] == '24_u-cc324': label1 = 'UM_RA2T_' + out_dir1[-4:-1]
    if out_dir1[:10] == '23_u-cc278': label1 = 'UM_CASIM-100_GA6alb'

    label2 = 'undefined_label'
    if out_dir2[:10] == '25_u-cc568': label2 = 'UM_RA2M'
    if out_dir2[:10] == '24_u-cc324': label2 = 'UM_RA2T_' + out_dir2[-4:-1]
    if out_dir2[:10] == '23_u-cc278': label2 = 'UM_CASIM-100_GA6alb'

    label3 = 'undefined_label'
    if out_dir3 == 'OUT_25H/': label3 = 'ECMWF_IFS'
    if out_dir3[:10] == '25_u-cc568': label3 = 'UM_RA2M'
    if out_dir3[:10] == '24_u-cc324': label3 = 'UM_RA2T_' + out_dir3[-4:-1]
    if out_dir3[:10] == '23_u-cc278': label3 = 'UM_CASIM-100_GA6alb'

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
