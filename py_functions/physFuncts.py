"""
Functions to calculate physical properties
==============================

"""

import numpy as np
from IPython import embed

# from __future__ import print_function

def calcAirDensity(temperature, pressure):

    """
    Function to calculate air density from temperature and pressure
    ==============================

    """

        #### EXAMPLE OF USE:
        #### data = calcAirDensity(data['temperature'][:], data['pressure'][:])
        ####        temperature = K
        ####        pressure = hPa

    R = 2.8704  #### hPa kg-1 K-1

    print('Calculating air density profile:')
    print('')
    #converting to np array
    if not isinstance(pressure,np.ndarray):
        pressure= np.array(pressure)
    if not isinstance(temperature, np.ndarray):
        temperature= np.array(temperature)

    #converting to hpa and K
    temperature[temperature<200] = temperature[temperature<200]+273.15
    pressure[pressure>1100] = pressure[pressure>1100]/100

    rho = pressure / (R*temperature)
    # ### if temperature is 1D
    # if np.ndim(temperature) == 1:
    #     rho = np.zeros([np.size(temperature)])
    #     for k in range(0,np.size(temperature)):
    #         rho[k] = pressure[k] / (R * temperature[k])
    # ### if temperature is 2D
    # elif np.ndim(temperature) == 2:
    #     rho = np.zeros([np.size(temperature,0), np.size(temperature,1)])
    #     for k in range(0,np.size(temperature, 1)):
    #         rho[:,k] = pressure[:,k] / (R * temperature[:,k])

    # print rho

    return rho

def calcThetaE(temperature, pressure, q):

    """
    Function to calculate equivalent potential temperature
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = q / qvs
                #### gives RH as a fraction

    print('Calculating theta:')
    theta = temperature * np.power(1e5 / pressure, (Rd/cp))
    print('...')

    print('Calculating theta of dry air:')
    thetad = temperature * np.power(1e5 / (pressure - evs), kd)
    print('...')

    print('Calculating theta_e:')
    tempvar = (-1.0 * kd * q) / eps
    thetaE = thetad * np.power( rh, tempvar ) * np.exp(L_vap * q / (temperature * cpd) )         ###Bryan 2008

    print('...')
    print('Done!')

    return theta, thetaE

def calcThetaVL(temperature, pressure, q, ql, qi, tim, height):

    """
    Function to calculate virtual liquid potential temperature
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg
    liquid mass mixing ratio = kg/kg
    ice mass mixing ratio = kg/kg

    Note that theta_l is based on 'liquid/frozen water static energy' (= cp.T + g.z - L.ql - Ls.qi)
        rather than potential temperature.
    Theta_vl is a conserved variable that is equal to virtual potential temperature (theta_v) in
        cloud-free air and so is used as a simplified measure of buoyancy

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    g = 9.806       # m/s2

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### total water mixing ratio
    qt = q + ql + qi

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = q / qvs
                #### gives RH as a fraction

    print('Calculating theta:')
    theta = temperature * np.power(1e5 / pressure, (Rd/cp))
    print('...')

    print('Calculating theta_l:')
    theta_l = temperature - ((L_vap * ql)/cp) - ((L_sub * qi)/cp) + ((g * height)/cp)
    print('...')

    print('Calculating theta_vl:')
    cv = (1/eps) - 1
    theta_vl = theta_l + theta_l * cv * qt
    print('...')

    print('...')
    print('Done!')

    return theta, theta_l, theta_vl

def calcsvp(T):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature in K or C
    output:
    SVP in hPa

    """
    #converting K
    tempC = T
    tempC[tempC>200] = tempC[tempC>200] - 273.15

    satvappres = 6.112 * np.exp( 17.67*tempC / (tempC + 243.5) )

    return satvappres

def calcvp(T):

    """
    Function to calculate  vapour pressure
    ==============================
    inputs:
    temperature = K or C
    if T=ambient then vp=saturation vapour pressure at T
    if T=DP then vp=ambient vapour pressure at ambient air temp

    """
    #converting K
    T[T<100] =T[T<100] + 273.15

    vappres =np.exp((-6867.7/T)-(5.2952*np.log(T))+56.658)

    return vappres

def polysvp(t,type):

    """
    Function to calculate saturation vapour pressure
    ==============================
    inputs:
    temperature = K

    POLYSVP RETURNED IN UNITS OF PA.
    TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)
    COPIED FROM PARCEL MODEL
    """

    ### ! ICE
    if type == 1:
        dt = np.maximum(-80,t-273.16)
        svp = 6.11147274 + dt * (0.503160820 + dt * (0.188439774e-1 + dt * (0.420895665e-3 + dt *
            (0.615021634e-5 + dt * (0.602588177e-7 + dt * (0.385852041e-9 + dt * (0.146898966e-11 +
            0.252751365e-14 * dt)))))))
        svp = svp*100

    ### ! LIQUID
    elif type == 0:
        dt = np.maximum(-80,t-273.16)
        svp = 6.11239921 + dt * (0.443987641 + dt * (0.142986287e-1 + dt * (0.264847430e-3 + dt *
            (0.302950461e-5 + dt * (0.206739458e-7 + dt * (0.640689451e-10 + dt * (-0.952447341e-13 +
            -0.976195544e-15 * dt)))))))
        svp = svp*100


    return svp

def calcRH(temperature, pressure, q):

    """
    Function to calculate RH from given water vapour mixing ratio
    also works with specific humidity
    ==============================
    inputs:
    pressure = Pa
    temperature = K
    water vapour mixing ratio = kg/kg

    """

    L_vap = 2.555e6    # J/kg
    L_sub = 2.836e6  # J/kg
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K

    Rd = 287.04   # dry air J kg^-1 K^-1
    Rv = 461.50

    kd = Rd/cpd     # k dry air

    eps = Rd/Rv

    ### saturation vapour pressue
    evs = polysvp(temperature, 0)

    ### saturation mixing ratio and relative humidity
    qvs = (eps * evs) / (pressure - evs)
    rh = (q / qvs) * 100
                #### gives RH as a percentage

    return rh

def calcT(theta,pressure):
    cp = 1004.6      # J/kg.K
    cpd = 1005.7     # J/kg.K
    Rd = 287.04   # dry air J kg^-1 K^-1
    p0 = 1e3   # hPa

    #converting to np array
    if not isinstance(pressure,np.ndarray):
        pressure= np.array(pressure)
    if not isinstance(theta, np.ndarray):
        theta= np.array(theta)

    #converting to hpa and K
    pressure[pressure>1100] = pressure[pressure>1100]/100
    theta[theta<200] = theta[theta<200]+273.15

    temperature = theta / np.power(p0 / pressure,Rd/cp)
    return temperature

def calcDewPoint(mr,p ):

    """
    Function to calculate  dewpoint
    equation & constants from Roger & Yau, 'Short course in cloud physics', Eqn 2.28
    ==============================
    inputs:
    pressure in hPa
    vapour mixing ratio (mr)  in g/kg or kg/kg
    """

    A=2.53e9
    B=5.42e3
    epsilon=0.622

    # convert to hPa if necessary
    if (p>9000).any():
        p=p/100
    #convert Q to kg/kg
    if (mr>1).any():
        mr=mr/1000
    # check Q for -ve values (yes it has happened!!)
    mr[mr<0]=mr[mr<0]*0+1e-6

    dp=B/np.log(A*epsilon/(mr*p))

    return dp

def calcSH_mr(mr,p):

    """
    Function to calculate  specific humidity
    ==============================
    inputs:
    vapour mixing ratio in in g/kg or kg/kg
    pressure = hpa
    """
    Td = calcDewPoint(mr,p)
    wvp = calcvp(Td)  # vp at Td is actual vapour pressure
    sh=0.622*wvp/(p-0.378*wvp)*1000

    return sh

def calcSH_wvp(wvp,p):

    """
    Function to calculate  specific humidity
    ==============================
    inputs:
    water vapour pressure hPa
    pressure = hpa
    """
    #convert vp to hpa
    if (wvp>100).any():
        wvp=wvp/100

    # convert to hPa if necessary
    if (p>9000).any():
        p=p/100

    sh=0.622*wvp/(p-0.378*wvp)*1000

    return sh

def calcP(T,Theta):

    cpd = 1005.7     # J/kg.K
    Rd = 287.04   # dry air J kg^-1 K^-1

    kd = Rd/cpd     # k dry air

    if (T<100).any():
        T+=273.15

    if (Theta<100).any():
        Theta+=273.15

    p0=1000
    p=p0/np.power(Theta/T,1/kd)

    return p

def windcomp2windvec(u_ms,v_ms):
    # windcomp2windvec.m will convert wind komponents into horizontal wind
    # velocitiy and direction
    # input: u,v
    # output: wsp, wd

    #Check to make sure sizes of the vectors are the same
    #if (uc.shape != vc.shape):
    #    print('u vector and v vector must be the same size')

    #calculate wind velocity
    wsp = np.sqrt(u_ms*u_ms +v_ms*v_ms)
    wind_dir_trig = np.arctan2(u_ms/wsp, v_ms/wsp)  # arctan needs normalized wind components
    wind_dir_trig_to_degrees = wind_dir_trig * 180/np.pi ## convert back to degrees from radians
    wind_dir_trig_from_degrees = wind_dir_trig_to_degrees + 180 ## convert to met. wind directions
    wd = wind_dir_trig_from_degrees
    wd[wd==0]=360.
    return wsp,wd
