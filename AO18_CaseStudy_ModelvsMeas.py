### SCRIPT TO COMPARE MODEL TO MEASURMENTS
from __future__ import print_function
import time
import datetime
import numpy as np
from netCDF4 import Dataset


#### import python functions
import sys
sys.path.insert(1, './py_functions/')
from time_functions import datenum2date, date2datenum, calcTime_Mat2DOY, calcTime_Date2DOY
from readMAT import readMatlabStruct
from physFuncts import calcThetaE, calcThetaVL
from pyFixes import py3_FixNPLoad
def main():

    START_TIME = time.time()
    print ('******')
    print ('')
    print ('Start: ' + time.strftime("%c"))
    print ('')

    ### INPUT FOLDERS
    um_root_dir = '/nfs/a96/MOCCHA/working/juta/CloudNet/Oden/data/calibrated/metum/'
    obs_root_dir = '/home/gillian/MOCCHA/ODEN/DATA/'
        #ship_filename = '/nfs/a96/MOCCHA/ODEN/DATA/2018_shipposition_1hour.txt'

    ### CHOSEN RUN
    out_dir1 = '25_u-cc568_RA2M_CON/'
    out_dir2 = '23_u-cc278_RA1M_CASIM/'
    out_dir3 = '24_u-cc324_RA2T_CON/'

    ### CHOOSE DATES TO PLOT

    DATES = 20180913

    print ('******')
    print ('')
    print ('Identifying .nc file: ')
    print ('')


    ### -------------------------------------------------------------------------
    ### define input filename
    ### -------------------------------------------------------------------------

    strdate = str(DATES)
    datenum = date2datenum(datetime.datetime(DATES))

    filename_um1 = um_root_dir + out_dir1 + strdate + 'metum.nc'
    filename_um2 = um_root_dir + out_dir2 + strdate + 'metum.nc'
    filename_um3 = um_root_dir + out_dir3 + strdate + 'metum.nc'

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
                data3[var_list3[j]] = nc4.variables[var_list3[j]][:]
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

    label4 = 'undefined_label'
    if out_dir4[:10] == '25_u-cc568': label4 = 'UM_RA2M'
    if out_dir4[:10] == '24_u-cc324': label4 = 'UM_RA2T_' + out_dir4[-4:-1]
    if out_dir4[:10] == '23_u-cc278': label4 = 'UM_CASIM-100_GA6alb'

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
