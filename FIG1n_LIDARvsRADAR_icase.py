#Case studies of sub-cloud and cloud regions....with theta profiles..
#input dir is
from __future__ import division
import matplotlib
matplotlib.use('agg')
import pygtk
pygtk.require('2.0')
import gtk
import numpy as np
import os, glob, sys
import scipy.io
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import datetime, time
from datetime import date
import time
from scipy.stats import skew
from scipy.stats import kurtosis
import math
from netCDF4 import Dataset
from datetime import date, timedelta
import scipy.io as sio
from scipy import stats

def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print "Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds."
    else:
        print "Toc: start time not set"

tic()
idate = sys.argv[1]
print idate
tDay = datetime.datetime.strptime(idate, "%Y%m%d")
#jdate = (tDay + datetime.timedelta(days=-1)).strftime('%Y%m%d')
#print jdate
print "idate and len idate:", idate, len(idate)
tt1 = sys.argv[2]
tt2 = sys.argv[3]
print "time lengths of scloud", idate, 't1',int(tt1), 't2', int(tt2)
if idate=='20180831':
  san=1
else:
  san=0

print "cloud base height correction:", san, idate
#PART-0....retrieving the radiosonde data.........now..of idate
#directories...
SONDEdir = '/nfs/a96/MOCCHA/working/jutta/UserReadyData/radiosondes/'
outdir = os.getcwd() + '/'

#PART1 Reading..the mean wind speed profiles from RADIOSONDE data....
#sonde = loadmat('SondeData_hint_V02.mat')
sonde = loadmat(SONDEdir+'SondeData_hint_V02.mat', squeeze_me=True, struct_as_record=False)
rwspd = sonde['RSintall'].ws   ##wind speed in m/s
rmday = sonde['RSintall'].mday ##day in matlab format.
rht   = sonde['RSintall'].Z    ##height in meter
rrh   = sonde['RSintall'].RH    ##relative humidity in %
rsh   = sonde['RSintall'].sphum    ##specific humidity in g/kg
rmr   = sonde['RSintall'].mr    ##mixing ratio in g/kg
rmrs  = sonde['RSintall'].mrs    ##saturation mixing ratio in g/kg
rvp   = sonde['RSintall'].vp    ##vapour pressure in
rsvp  = sonde['RSintall'].svp    ##saturation vapour pressure
rpt   = sonde['RSintall'].pottemp    ##potential temperature in degC
rept  = sonde['RSintall'].epottemp    ##equivalent potential temperature in degC
rdpt  = sonde['RSintall'].dewp    ##dewpoint temperature in degC
rtemp = sonde['RSintall'].pottemp    ##temperature in degC
rpres = sonde['RSintall'].pressure    ##pressure in hPa

print "RS: shapes:", rwspd.shape, rmday.shape, rht.shape
print "temp, ept, sh, rh:", rtemp.shape, rept.shape, rsh.shape, rrh.shape

rmdays = np.zeros((len(rmday),8))
for m in range(0,len(rmday)):
  rmdays[m,0] = int(date.fromordinal(int(rmday[m]-366)).year)
  rmdays[m,1] = int(date.fromordinal(int(rmday[m]-366)).month)
  rmdays[m,2] = int(date.fromordinal(int(rmday[m]-366)).day)
  rmdays[m,3] = int(rmday[m]%1*24-rmday[m]%1*24%1)
  rmdays[m,4] = int(rmday[m]%1*24%1*60-rmday[m]%1*24%1*60%1)
  rmdays[m,5] = int(np.round(rmday[m]%1*24%1*60%1*60))
  rmdays[m,6] = float(rmday[m]%1*24)
  rmdays[m,7] = float(rmdays[m,1] + rmdays[m,2]/100)

print rmdays.shape
##Selecting input array date data only...from total set of radiosonde data....
idm = int(idate[4:6]) + int(idate[6:8])/100
inx = np.where(rmdays[:,7] == idm)
if len(inx[0])>=1:
 RSpres =  rpres[:,inx[0]] #pressure in hPa
 RSwspd =  rwspd[:,inx[0]]
 RStemp =  rtemp[:,inx[0]]
 RSept  =  rept[:,inx[0]]
 RSpt  =  rpt[:,inx[0]]
 RSsh   =  rsh[:,inx[0]]
 RSrh   =  rrh[:,inx[0]]
 RStim  =  rmdays[inx[0],6]
 RSht   =  rht/1000
 RSinx  = 1
else:
 RSinx = 0
 print "No Radiosonde data is found ", idate

#del idm; del inx
print "running date shapes:", len(RStim), RSht.shape, RSpres.shape, RStemp.shape, RSept.shape, RSsh.shape, RSrh.shape, RSpt.shape
#sys.exit()

##Reading inversion height file...
SONDEdir = '/nfs/a96/MOCCHA/working/araveti/radisonde/'
#inversion = loadmat(SONDEdir + 'InversionHeights_RSh05int_final_V03.mat', squeeze_me=True, struct_as_record=False)
inversion = loadmat(SONDEdir + 'InversionHeights_RSh05int_final_v2.0.mat', squeeze_me=True, struct_as_record=False)
imday = inversion['RSdec'].mday ##day in matlab format.
invbase = inversion['RSdec'].invbase   ##main inversion base height (m)
invdepth = inversion['RSdec'].invdepth   ##depth of inversion (m)
invstrength = inversion['RSdec'].invstrength   ##main inversion strength (K)
decbase = inversion['RSdec'].decbase   ##decoupling height (m)
sfmlbase = inversion['RSdec'].sfmlbase   ##surface mixed layer height (m)
bltype = inversion['RSdec'].bltype   ##boundary layer type:, 0-coupled, 1-surface inv, 2-decoupling by temp inv, 3-decoupling by stable layer below temp inv, 4-no temp inv
qinvbase = inversion['RSdec'].q_invbase   #q#main inversion base height (m)
qinvdepth = inversion['RSdec'].q_invdepth   #q#depth of inversion (m)
qinvstrength = inversion['RSdec'].q_invstrength   #q#main inversion strength (K)
qdecbase = inversion['RSdec'].q_decbase   #q#decoupling height (m)
print "RS: inversion shapes,", bltype.shape, invbase.shape, decbase.shape, sfmlbase.shape
print "RS: q-inversion shapes,", qinvbase.shape, qdecbase.shape

##Selecting input array date data only...from total set of radiosonde data....
if len(inx[0])>=1:
 RSbltype =  bltype[inx[0]]
 RSinvbase =  invbase[inx[0]]
 RSdecbase  =  decbase[inx[0]]
 RSqinvbase   =  qinvbase[inx[0]]
 RSqdecbase   =  qdecbase[inx[0]]
 RSsfmlbase   =  sfmlbase[inx[0]]
 RSinx  = 1
else:
 RSinx = 0
 print "No Radiosonde data is found ", idate

del idm; del inx
print "running date lengths of inversion heights:RSdec:", len(RStim), len(RSbltype), len(RSinvbase), len(RSdecbase), len(RSqinvbase), len(RSqdecbase), len(RSsfmlbase)
del rmdays; del sonde; del rwspd; del rtemp; del rept; del rsh; del rrh; del rht; 
del rvp; del rsvp; del rpt; del rmr; del rmrs;
del inversion; del imday; del invbase; del invdepth; del invstrength; del decbase; del sfmlbase; del bltype; 
del qinvbase; del qinvdepth; del qinvstrength; del qdecbase; 
#sys.exit()


#ORIGINAL.....
#directories...
indir = '/nfs/a96/MOCCHA/working/araveti/lidarvsradar/tvariance/'
outdir = os.getcwd() + '/'
print "DATA directory......Running directory.......Script directory....all:", indir, outdir
#Part-A selecting the well mixed day first.................
VARepsL=[];  VARepsR=[]; VARerrL=[];  VARerrR=[];  FlagL=[];  FlagR=[]
VARsnrL=[];  VARsnrR=[]; VARsigL=[];  VARsigR=[];  VARdopL=[];  VARdopR=[];
VARskwL=[];  VARskwR=[]; VARbetL=[];  VARspwR=[];
VARlsL1=[];  VARlsLs=[];  VARlsR1=[];  VARlsRs=[];
for day in range(1):
  jdate = idate
  print "-------------------------------------------"
  print "processing the run date: CASE study", jdate
  print "-------------------------------------------"
  vfile = indir + 'LIDARattributes_'+jdate+'.mat'
  vfile2 = indir + 'RADARattributes_'+jdate+'.mat'
  if (os.path.exists(vfile) and os.path.exists(vfile2)):
    #lidar
    vdata = sio.loadmat(vfile)
    varepsL = vdata['epsilon_L']
    vardopL = vdata['LDOP']
    varsnrL = vdata['LSNR']
    varbetaL = vdata['LBETA']
    varsigdopL = vdata['sigma_dopL']
    varskewdopL = vdata['skew_dopL']
    varrngL = vdata['Lranges']
    vartime = vdata['atime']
    varerrL = vdata['epsilon_Lepserr']
    varerrL = np.abs(varerrL);
    #length scales
    varlsL1 = vdata['L1_lid']
    varlsLs = vdata['Ls_lid']
    print "LIDAR: attributes + epsilon shapes:", varepsL.shape, varrngL.shape, vartime.shape, varerrL.shape
    print "LIDAR: length scales:", varlsL1.shape, varlsLs.shape
    #radar
    vdata2 = sio.loadmat(vfile2)
    varepsR = vdata2['epsilon_R']
    vardopR = vdata2['RDOP']
    varsnrR = vdata2['RSNR']
    varspwR = vdata2['RSPW']
    varsigdopR = vdata2['sigma_dopR']
    varskewdopR = vdata2['skew_dopR']
    varrngR = vdata2['Rranges']
    vartime = vdata2['atime']
    varerrR = vdata2['epsilon_Repserr']
    varerrR = np.abs(varerrR);
    print "RADAR: attributes + epsilon shapes:", varepsR.shape, varrngR.shape, vartime.shape, varerrR.shape
    #length scales
    varlsR1 = vdata2['L1_rad']
    varlsRs = vdata2['Ls_rad']
    print "RADAR: length scales:", varlsR1.shape, varlsRs.shape
    ##
    for i in range(len(vartime)):
      for j in range(len(varrngL)):
        a = varepsL[i,j];  b = varerrL[i,j]
        if ((a>=0) or (a<=-8) or (b>=350)):
         varepsL[i,j]=np.nan; varerrL[i,j]=np.nan
         varsnrL[i,j]=np.nan; vardopL[i,j]=np.nan
         varbetaL[i,j]=np.nan; varsigdopL[i,j]=np.nan
         varskewdopL[i,j]=np.nan;
        del a; del b;

      for j in range(len(varrngR)):
        c = varepsR[i,j];  d = varerrR[i,j]
        if ((c>=0) or (c<=-8) or (d>=350)):
         varepsR[i,j]=np.nan; varerrR[i,j]=np.nan
         varsnrR[i,j]=np.nan; vardopR[i,j]=np.nan
         varspwR[i,j]=np.nan; varsigdopR[i,j]=np.nan
         varskewdopR[i,j]=np.nan;
        del c; del d;

    del i; del j;
    ##NAN quality-lidar
    for i in range(len(vartime)):
      for j in range(len(varrngL)):
        if np.isnan(varepsL[i,j])==True:
         varerrL[i,j]=np.nan;  varsnrL[i,j]=np.nan;  vardopL[i,j]=np.nan;  varbetaL[i,j]=np.nan
         varsigdopL[i,j]=np.nan; varskewdopL[i,j]=np.nan;

    del i; del j;

    for i in range(len(vartime)):
      for j in range(len(varrngR)):
        if np.isnan(varepsR[i,j])==True:
         varerrR[i,j]=np.nan;  varsnrR[i,j]=np.nan;  vardopR[i,j]=np.nan;  varspwR[i,j]=np.nan
         varsigdopR[i,j]=np.nan; varskewdopR[i,j]=np.nan;

    del i; del j;
    #converting lidar data into radar data.................
    linx=np.where(varrngL>=2.9)[0][0];    rinx=np.where(varrngR>=2.9)[0][0]
    leps = varepsL[:,:linx];     lerr = varerrL[:,:linx]
    ldop = vardopL[:,:linx];     lsig = varsigdopL[:,:linx];    lskw = varskewdopL[:,:linx]
    lrng = varrngL[:linx];     rrng = varrngR[:rinx]
    varepsLnew = np.asarray(np.zeros((len(vartime), len(rrng))));  varepsLnew[:] = np.nan;
    varerrLnew = np.asarray(np.zeros((len(vartime), len(rrng))));  varerrLnew[:] = np.nan;
    vardopLnew = np.asarray(np.zeros((len(vartime), len(rrng))));  vardopLnew[:] = np.nan;
    varsigLnew = np.asarray(np.zeros((len(vartime), len(rrng))));  varsigLnew[:] = np.nan;
    varskwLnew = np.asarray(np.zeros((len(vartime), len(rrng))));  varskwLnew[:] = np.nan;
    for r in range(len(rrng)):
      rng = rrng[r][0]
      inx = np.where((lrng>=rng-0.01) & (lrng<=rng+0.01))[0]
      for t in range(len(vartime)):
        eps = leps[t,inx];  err = lerr[t,inx];  dop = ldop[t,inx];  sig = lsig[t,inx];  skw = lskw[t,inx]
        varepsLnew[t,r] = eps[~np.isnan(eps)].mean()
        varerrLnew[t,r] = err[~np.isnan(err)].mean()
        vardopLnew[t,r] = dop[~np.isnan(dop)].mean()
        varsigLnew[t,r] = sig[~np.isnan(sig)].mean()
        varskwLnew[t,r] = skw[~np.isnan(skw)].mean()
        del eps; del err; del dop; del sig; del skw;
      del rng; del inx;
    print "samples lidar epsilon shape:", idate,'--', varepsLnew.shape, varerrLnew.shape, vardopLnew.shape, varsigLnew.shape, varskwLnew.shape
    varrng =  rrng;
    del linx; del rinx; del leps; del lerr; del rrng; del ldop; del lsig; del lskw;

    #store1-from file1...
    VARepsL.append(varepsL);  VARsnrL.append(varsnrL);  VARerrL.append(varerrL);  VARdopL.append(vardopL)
    VARbetL.append(varbetaL);  VARsigL.append(varsigdopL);  VARskwL.append(varskewdopL)
    VARlsL1.append(varlsL1);   VARlsLs.append(varlsLs);
    #store2--from file2...
    VARepsR.append(varepsR);  VARsnrR.append(varsnrR);  VARerrR.append(varerrR);  VARdopR.append(vardopR)
    VARspwR.append(varspwR);  VARsigR.append(varsigdopR);  VARskwR.append(varskewdopR)
    VARlsR1.append(varlsR1);   VARlsRs.append(varlsRs);
    del vdata; del vdata2;

print "EPSILON--length of total files:", len(VARepsL), len(VARepsR)
print "EPS-ERROR--length of total files:", len(VARerrL), len(VARerrR)
print "MEAN-SNR--length of total files:", len(VARsnrL), len(VARsnrR)
print "MEAN-DOP--length of total files:", len(VARdopL), len(VARdopR)
print "length scales LIDAR vs RADAR:", len(VARlsL1), len(VARlsLs), len(VARlsR1), len(VARlsRs)
toc()
print "PART-A is done:  all dates post-process is done...IMP variables are extracted"
print "converted LIDAR data shapes:", varepsLnew.shape, varerrLnew.shape, vardopLnew.shape, varsigLnew.shape, varskwLnew.shape
sys.exit()


#PART-B....#Reading the micorwave radiometer data...
#directories...
MWRdir = '/nfs/a96/MOCCHA/working/araveti/DISSIPATION/'
outdir = os.getcwd() + '/'

#PART1 Reading..the mean temperature profiles from MWR data....
mwr = loadmat(MWRdir + 'HATPRO_Tprofiles_adjusted_inversionheights.mat', squeeze_me=True, struct_as_record=False)
mmday = mwr['HTcomp'].time   #matlab time
altitude = mwr['HTcomp'].altitude  #altitude in m
temp = mwr['HTcomp'].temperature  #temperature in Kelvin
#
invbase = mwr['HTcomp'].invbase   ##main inversion base height (m)
invbasetemp = mwr['HTcomp'].invbasetemp   ##main inversion base height temperature (K)
invtop = mwr['HTcomp'].invtop   ##main inversion top height (m)
invdepth = mwr['HTcomp'].invdepth   ##depth of inversion (m)
invstrength = mwr['HTcomp'].invstrength   ##main inversion strength (K)
decbase = mwr['HTcomp'].decbase   ##decoupling height (m)
decoupled = mwr['HTcomp'].decoupled   ##decoupling number yes or no...1 is yes.

print "MWR: shapes:", mmday.shape, altitude.shape, temp.shape
print "MWR: shapes stats:", invbase.shape, invtop.shape, invdepth.shape, invbasetemp.shape, invstrength.shape, decbase.shape, decoupled.shape

#microwave radiometer matlab time conversion. dates...
mmdays = np.zeros((len(mmday),10))
for m in range(0,len(mmday)):
  mmdays[m,0] = int(date.fromordinal(int(mmday[m]-366)).year)
  mmdays[m,1] = int(date.fromordinal(int(mmday[m]-366)).month)
  mmdays[m,2] = int(date.fromordinal(int(mmday[m]-366)).day)
  mmdays[m,3] = int(mmday[m]%1*24-mmday[m]%1*24%1)
  mmdays[m,4] = int(mmday[m]%1*24%1*60-mmday[m]%1*24%1*60%1)
  mmdays[m,5] = int(np.round(mmday[m]%1*24%1*60%1*60))
  mmdays[m,6] = float(mmday[m]%1*24)
  mmdays[m,7] = float(mmdays[m,1] + mmdays[m,2]/100)
  mmdays[m,8] = int(date.fromordinal(int(mmday[m]-366)).timetuple().tm_yday)
  mmdays[m,9] = mmdays[m,8] + (mmdays[m,3]*60 + mmdays[m,4])/1440

print mmdays.shape
mwrtime = mmdays[:,0:5]
##considering only from 20180812 to 20180919............
#start search # mwrtime[np.where(mwrtime[:,2]==12)[0],:][0]
#end search   # mwrtime[np.where(mwrtime[:,2]==19)[0],:][-1]
for i in range(len(mwrtime)):
  if ((mwrtime[i,0]==2018) and (mwrtime[i,1]==8) and (mwrtime[i,2]==12) and (np.round(mwrtime[i,3],2)==0) and (np.round(mwrtime[i,4],2)==0)):
     t1=i
     print "MWR: start time index is :", i,

  if ((mwrtime[i,0]==2018) and (mwrtime[i,1]==9) and (mwrtime[i,2]==19) and (np.round(mwrtime[i,3],2)==23) and (np.round(mwrtime[i,4],2)==59)):
     t2=i
     print "MWR: end time index is :", i,

print "MWR: start and end time indexes are:", t1, t2

##Sub set of total microwave radiometer temperature data....
mtime = mmdays[t1:t2+1,:]
mtemp = temp[:,t1:t2+1]
minvbase = invbase[t1:t2+1]
minvtop = invtop[t1:t2+1]
minvdepth = invdepth[t1:t2+1]
minvstrength = invstrength[t1:t2+1]
minvbasetemp = invbasetemp[t1:t2+1]
mdecbase = decbase[t1:t2+1]
mdecoupled = decoupled[t1:t2+1]
print "MWR: Subset Study period shapes stats:", minvbase.shape, minvtop.shape, minvdepth.shape, minvbasetemp.shape, minvstrength.shape, mdecbase.shape, mdecoupled.shape, 'time-alt-temp', mtime.shape, altitude.shape, mtemp.shape

#overlay the invbase and invtop.....height regions............
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

IDMstore=[]; PROFtemp=[]; PROFinvbase=[]; PROFinvtop=[]; PROFdecoupled=[]; s=0;
for day in range(1):
  jdate = idate
  print "-------------------------------"
  print "processing the run date-MWR:", jdate
  print "-------------------------------"
  idm = int(jdate[4:6]) + int(jdate[6:8])/100;
  inx = np.where(mtime[:,7] == idm)[0]
  if len(inx)>1:
    IDMstore.append(float(idm))
    s=s+1;
    MWRtim = mtime[inx, 6]
    MWRtemp = mtemp[:,inx]
    MWRinvbase  = minvbase[inx]
    MWRinvtop  = minvtop[inx]
    MWRdecoupled = mdecoupled[inx]
    MWRdecbase = mdecbase[inx];
    MWRpot = np.zeros((len(altitude), len(MWRtim))); MWRpot[:]=np.nan;
    #original temp convert to potential temp.....MWRtemp....
    for k in range(len(MWRtim)):
      #pot-temp calc....
      nval1,nidx1 = find_nearest(RStim, MWRtim[k]);
      prs = RSpres[:,nidx1]
      #now taking the heights of mwr heights...nearest pot, prs.....
      nprs=[]; 
      for j in range(len(altitude)):
         nval2,nidx2 = find_nearest(RSht, altitude[j]/1000)
         nprs.append(float(prs[nidx2]))
         del nval2; del nidx2;
      nprs = np.asarray(nprs)
      MWRpot[:,k] = MWRtemp[:,k]*((1000/nprs)**0.286) 
      del nval1; del nidx1; del prs; del nprs;

    print "jdate, temp, mwr:", jdate, MWRtim.shape, MWRtemp.shape, MWRinvbase.shape, MWRinvtop.shape, MWRdecoupled.shape, MWRpot.shape
    MWRdecoupled = np.nan_to_num(MWRdecoupled)
    #interpolate2d....converting to 720 pixels
    mwrtime = np.arange(0,24,2/60)
    mwrtemp = np.zeros((len(altitude), len(mwrtime))); mwrtemp[:]=np.nan;
    mwrwnd = np.zeros((len(altitude), len(mwrtime))); mwrwnd[:]=np.nan;
    mwrpot = np.zeros((len(altitude), len(mwrtime))); mwrpot[:]=np.nan;
    mwrpot1 = np.zeros((len(altitude), len(mwrtime))); mwrpot1[:]=np.nan;
    mwrinvbase = np.zeros(len(mwrtime)); mwrinvbase[:]=np.nan;
    mwrinvtop = np.zeros(len(mwrtime)); mwrinvtop[:]=np.nan;
    mwrdecoupled = np.zeros(len(mwrtime)); mwrdecoupled[:]=np.nan;
    for k in range(len(mwrtime)):      
      val,idx = find_nearest(MWRtim, mwrtime[k]);   #nearest time...      
      mwrtemp[:,k]  =  MWRtemp[:,idx]
      mwrinvbase[k] =  MWRinvbase[idx]
      mwrinvtop[k]  =  MWRinvtop[idx]
      mwrdecoupled[k]  =  MWRdecoupled[idx]
      #pot-temp calc....
      val1,idx1 = find_nearest(RStim, mwrtime[k]);
      pot = RSpt[:,idx1]
      prs = RSpres[:,idx1]
      wnd = RSwspd[:,idx1]
      #now taking the heights of mwr heights...nearest pot, prs.....
      nprs=[]; npot=[]; nwnd=[];
      for j in range(len(altitude)):
         val2,idx2 = find_nearest(RSht, altitude[j]/1000)
         nprs.append(float(prs[idx2]))
         npot.append(float(pot[idx2]))
         nwnd.append(float(wnd[idx2]))
         del val2; del idx2;

      #print len(nprs), len(npot)
      nprs = np.asarray(nprs); npot = np.asarray(npot); nwnd = np.asarray(nwnd);
      mwrpot[:,k] = MWRtemp[:,idx]*((1000/nprs)**0.286)
      mwrpot1[:,k] = npot;
      mwrwnd[:,k] = nwnd;
      del nprs; del npot; del pot; del prs; del wnd;
      del val; del idx; del val1; del idx1;
      #sys.exit()

  del inx; del idm; del jdate;

print "TOTAL: MWR temp profiles: idate", idate, MWRtemp.shape, MWRdecoupled.shape, MWRtim.shape, altitude.shape
print "nearest-interpolated MWR shapes:", mwrtime.shape, mwrtemp.shape, mwrinvbase.shape, mwrinvtop.shape, mwrdecoupled.shape
del mtime; del mtemp; del minvbase; del minvtop; del mdecoupled; del minvdepth; del minvstrength;
print "Potential temperature +wind profiles...:", idate, mwrpot.shape, mwrpot1.shape, mwrwnd.shape
#sys.exit()


#PART-C.........cloud base height.............well mixed vs de coupled boudary layers......
##CEILOmeter data reading...CLOUD BASE HEIGHT>.....
CEILOdir = '/nfs/a96/MOCCHA/working/data/AO2018_WX_v1_1/'
outdir = os.getcwd() + '/'
#Reading..the ship ceilometer. data....reqd variables...select from .rtf file........
ceilo = loadmat(CEILOdir + 'ACAS_AO2018_CL31_1min_v1_0.mat', squeeze_me=True, struct_as_record=False)
cdoy = ceilo['cl31'].doy   ##day of the year
cmday = ceilo['cl31'].mday  ##matlab day
cbase_ht =  ceilo['cl31'].base_ht   ##Base height of (1-3) layers. NaN if no layer detected.\
crange = ceilo['cl31'].ceil_range  ##range in m...~770-Ranges for the backscatter profile
cbsp = ceilo['cl31'].bs_prof  ##1/(km*steradians)       Backscatter coefficient profile.\
##ceilo-times...converting to ...table...
cmdays = np.zeros((len(cmday),10))
for m in range(0,len(cmday)):
  if np.isnan(cmday[m])==False:
   cmdays[m,0] = int(date.fromordinal(int(cmday[m]-366)).year)
   cmdays[m,1] = int(date.fromordinal(int(cmday[m]-366)).month)
   cmdays[m,2] = int(date.fromordinal(int(cmday[m]-366)).day)
   cmdays[m,3] = int(cmday[m]%1*24-cmday[m]%1*24%1)
   cmdays[m,4] = int(cmday[m]%1*24%1*60-cmday[m]%1*24%1*60%1)
   cmdays[m,5] = int(np.round(cmday[m]%1*24%1*60%1*60))
   cmdays[m,6] = float(cmday[m]%1*24)
   cmdays[m,7] = float(cmdays[m,1] + cmdays[m,2]/100)
   cmdays[m,8] = int(date.fromordinal(int(cmday[m]-366)).timetuple().tm_yday)
   cmdays[m,9] = cmdays[m,8] + (cmdays[m,3]*60 + cmdays[m,4])/1440

print "CEILOmeter data shapes: all campagin:", cdoy.shape, cmday.shape, cbase_ht.shape, crange.shape, cbsp.shape, cmdays.shape
#
crang = crange/1000
jdate = idate
adate = datetime.datetime.strptime(jdate,"%Y%m%d")
adoy = adate.timetuple().tm_yday
ainx = np.where(np.logical_and(cdoy>adoy, cdoy<=adoy+1))[0]
if len(ainx)==1440:
  cldbase =  cbase_ht[ainx,0]  ##ceilometer........blue-1
  cldbase1 =  cbase_ht[ainx,1]  ##ceilometer........red-2
  ctime = cmdays[:,6][ainx]
elif len(ainx)>1440:
  cldbase =  cbase_ht[ainx[:1440],0]
  cldbase1 =  cbase_ht[ainx[:1440],1]
  ctime = cmdays[ainx[:1440],6]
else:
  dummy = np.empty(1440-len(ainx))
  cldbase = np.append(cbase_ht[ainx,0],dummy)
  cldbase1 = np.append(cbase_ht[ainx,1],dummy)
  del dummy;

print "CEILOMETER data and storing the Cloud Base heights", jdate, cldbase.shape, adate, adoy, len(ainx), ctime.shape, 'cld1', cldbase1.shape
#ctimeWM = ctime; cldbaseWM = cldbase; del cldbase; del jdate; del adate; del adoy; del ainx; del ctime;


#Part-D...reading mwr -radiometer data....of idate.
##Part2
CEILOdir = '/nfs/a96/MOCCHA/working/data/AO2018_WX_v1_1/'
#PART2 Reading..the ship micro-wave raidio meter  data....reqd variables...select from .rtf file........
##lwp                     g/m^2           Liquid water path from HATPRO radiometer (preliminary value, use with caution)\
radio = loadmat(CEILOdir + 'ACAS_AO2018_WX_1min_v1_1.mat', squeeze_me=True, struct_as_record=False)
rdoy = radio['wx'].doy   ##day of the year
rmday = radio['wx'].mday  ##matlab day
rlwp = radio['wx'].lwp  ##Liquid water path from HATPRO radiometer (preliminary value, use with caution)\
rlwr = radio['wx'].pirn  ##W/m^2           Net longwave\
##radio-times...converting to ...table...
rmdays = np.zeros((len(rmday),10))
for m in range(0,len(rmday)):
  if np.isnan(rmday[m])==False:
   rmdays[m,0] = int(date.fromordinal(int(rmday[m]-366)).year)
   rmdays[m,1] = int(date.fromordinal(int(rmday[m]-366)).month)
   rmdays[m,2] = int(date.fromordinal(int(rmday[m]-366)).day)
   rmdays[m,3] = int(rmday[m]%1*24-rmday[m]%1*24%1)
   rmdays[m,4] = int(rmday[m]%1*24%1*60-rmday[m]%1*24%1*60%1)
   rmdays[m,5] = int(np.round(rmday[m]%1*24%1*60%1*60))
   rmdays[m,6] = float(rmday[m]%1*24)
   rmdays[m,7] = float(rmdays[m,1] + rmdays[m,2]/100)
   rmdays[m,8] = int(date.fromordinal(int(rmday[m]-366)).timetuple().tm_yday)
   rmdays[m,9] = rmdays[m,8] + (rmdays[m,3]*60 + rmdays[m,4])/1440

print "RADIOmeter data shapes: all campagin:", rdoy.shape, rmday.shape, rlwp.shape, rmdays.shape

#WM vs DC cloud base heights....RADIO-LWP
bdate = datetime.datetime.strptime(jdate,"%Y%m%d")
bdoy = bdate.timetuple().tm_yday
binx = np.where(np.logical_and(rdoy>bdoy, rdoy<=bdoy+1))[0]
if len(ainx)==1440:
  radlwp =  rlwp[binx]  ##radio liquid water path
  rtime = rmdays[:,6][binx]
elif len(binx)>1440:
  radlwp =  rlwp[binx[:1440]]
  rtime = rmdays[binx[:1440],6]
else:
  dummy = np.empty(1440-len(binx))
  radlwp = np.append(rlwp[binx],dummy)
  del dummy;

print "RADIOMETER data and storing the liquid water path", jdate, radlwp.shape, bdate, bdoy, len(binx), rtime.shape

##
##
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def calcSma(data, smaPeriod):
    j = next(i for i, x in enumerate(data) if x is not None)
    our_range = range(len(data))[j + smaPeriod - 1:]
    empty_list = [None] * (j + smaPeriod - 1)
    sub_result = [np.mean(data[i - smaPeriod + 1: i + 1]) for i in our_range]

    return np.array(empty_list + sub_result)

##figuring the cloud base heights wm vs dc
#FIGURE1....plotting the basic remote attributes between wellmixed vs decoupled boundary layers...
##20180914 - a well mixed case study..............reqd variables...
LBETA = VARbetL[0]
LSNR  = VARsnrL[0]
RSNR  = VARsnrR[0]
LDOP  = VARdopL[0]
RDOP  = VARdopR[0]
atime = vartime[:,0]
Lranges = varrngL[:,0]
Rranges = varrngR[:,0]
scldbase = smooth(cldbase/1000,5)
scldbase1 = smooth(cldbase1/1000,5)
#ctime = ctime;
m,n = LSNR.shape; LSNR1 = np.zeros((m,n)); LSNR1[:] = np.nan;
for i in range(m):
  for j in range(n):
    snr = LSNR[i,j]
    if np.isnan(snr)==False:
      LSNR1[i,j] = 10*np.log10(np.abs(snr-1))
    del snr;
del i; del j; del m; del n;
#LDOP correction only if 20180917 date 12-13 UTC................
if idate == '20180917':
  LDOP[360:400,:] = LDOP[360:400,:]/2
  LDOP[360:400,:] = LDOP[360:400,:]*-1
  LDOP[360:400,:] = LDOP[360:400,:]+0.1
  vardopLnew[360:400,:] = vardopLnew[360:400,:]/2
  vardopLnew[360:400,:] = vardopLnew[360:400,:]*-1

#quality control the data..
#eps
#LEPS[LBETA<=10**-6]=np.nan
LSNR[LBETA<=10**-6]=np.nan
LSNR1[LBETA<=10**6]=np.nan
LDOP[LBETA<=10**-6]=np.nan
#LBETA[LBETA<=10**-6]=np.nan
print "LSNR vs LSNR1 shapes:, counts:", LSNR.shape, LSNR1.shape, len(np.where(np.isnan(LSNR)==False)[0]), len(np.where(np.isnan(LSNR1)==False)[0])


#PART-E...estimation of cloud heights.........
#cloud top heights estimation - roughly
def max_repeatedNaNs(a):
    # Mask of NaNs
    mask = np.concatenate(([False],np.isnan(a),[False]))
    if ~mask.any():
        return 0
    else:
        # Count of NaNs in each NaN group. Then, get max count as o/p.
        c = np.flatnonzero(mask[1:] < mask[:-1]) - \
            np.flatnonzero(mask[1:] > mask[:-1])
        return c.max()

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


#estimation of cld top height........retriveal.........
rcth=[]; cldtop=[];
rcthn=[]; cldtopn=[];
for t in range(len(atime)):
  snr = RSNR[t,:]
  snrix = np.where(np.isnan(snr))[0]
  snrinx = consecutive(snrix)
  if len(snrinx)>1:
   #method1
   snrlen = max_repeatedNaNs(snr)
   totlens = [len(i) for i in snrinx]
   inx = np.where(totlens == snrlen)[0][0]
   rcth.append(int(snrinx[inx][0]))
   inx1 = snrinx[inx][0]
   cldtop.append(float(Rranges[inx1]))
   del snrlen; del totlens; del inx; del inx1;
   #method2
   tinx = np.where(ctime>atime[t])[0][0]
   icbh = scldbase[tinx]
   if np.isnan(icbh)==False:
     icth=[]; icthix=[];
     for j in range(len(snrinx)):
       jcth = Rranges[snrinx[j][0]]
       if jcth>icbh:
        icth.append(Rranges[snrinx[j][0]])
        icthix.append(int(snrinx[j][0]))
       del jcth;
     icth = np.asarray(icth)
     icthix = np.asarray(icthix)
     if len(icth)>0:
       inx = icthix[np.where(icth==np.min(icth))]
       rcthn.append(int(inx))
       cldtopn.append(float(Rranges[inx]))
     else:
       inx =0;
       cldtopn.append(cldtop[t])
       rcthn.append(rcth[t])
     del inx; del icth; del icthix;
   else:
     cldtopn.append(cldtop[t])
     rcthn.append(rcth[t])
   del tinx; del icbh;
  else:
   rcth.append(snrinx[0][0])
   cldtop.append(float(Rranges[snrinx[0][0]]))
   rcthn.append(snrinx[0][0])
   cldtopn.append(float(Rranges[snrinx[0][0]]))
   #print "equal-1",t
  del snr; del snrix; del snrinx;

cldinx = np.where(np.asarray(rcth)<=5)[0]
scldtop = np.asarray(cldtop)
scldtop[cldinx] = np.nan
scldtop1 = np.asarray(cldtopn)
scldtop1[cldinx] =  np.nan;
print "retrieved cloud top heights...alogrithm-self:", idate, scldbase.shape, scldtop.shape, scldtop1.shape

#Manual quality algorithms;;;;;;;;;;;;;
MTMP  = mwrtemp.T-273.15
MPOT  = mwrpot.T-273.15
MWND  = mwrwnd.T;
Mranges = altitude/1000
mtime = MWRtim
#MTMP correction only if 20180823 date 12-
if idate == '20180823':
  MTMP[360:525,:] = np.nan

print "MWR: temp profile....shape", mtime.shape, Mranges.shape, MTMP.shape, MPOT.shape, MWND.shape
#sys.exit()


##
#Giving the required time range..........
t1 = int(tt1)
t2 = int(tt2)
del tt1; del tt2;
if idate=='20180823':
  h1=0; h2=2.1; hh=0.3;
elif idate=='20180905':
  h1=0; h2=3.0; hh=0.5;
elif idate=='20180902':
  h1=0; h2=2.1; hh=0.3;
else:
  h1=0; h2=1.5; hh=0.3;
print "plotting the plot:", idate, t1, t2, h1, h2, hh


#Getting the lidar, radar, mwr attributes..............
LSIG  = VARsigL[0];  LSKW  = VARskwL[0]
RSIG  = VARsigR[0];  RSKW  = VARskwR[0]
LEPS  = VARepsL[0];  REPS  = VARepsR[0]
LEPSerr = VARerrL[0]; REPSerr = VARerrR[0]
#del LBETA;
#LBETA = VARbetL[0]
LEPS[LBETA<=10**-6]=np.nan
LBETA[LBETA<=10**-6]=np.nan
if ((idate=='20180913') or (idate=='20180901')):
 LEPS[LEPS>-2]=np.nan;
if ((idate=='20180917') or (idate=='20180913') or (idate=='20180906')):
 m,n=LEPS.shape
 for i in range(m):
  for j in range(12):
   if LEPS[i,j]>-3:
     LEPS[i,j]=np.nan;
     LDOP[i,j]=np.nan;
     LSKW[i,j]=np.nan;
     LEPSerr[i,j]=np.nan;

 del i; del j; del m; del n;

#work-A#PDF#Extracting the data of A#--w, w-sig, w-skew, eps, B#--cld-depth, lwp,
#work-B#surface-to-cldbase, cldbase-to-cldtop..........eps-lidar, eps-radar, mwr-t....
c1 = np.where(ctime>=t1)[0][0]; c2 = np.where(ctime<=t2)[0][-1]
a1 = np.where(atime>=t1)[0][0]; a2 = np.where(atime<=t2)[0][-1]
m1 = np.where(mtime>=t1)[0][0]; m2 = np.where(mtime<=t2)[0][-1]
cbase = scldbase[c1:c2+1][::2]
cbase1 = scldbase1[c1:c2+1][::2]
ctop = scldtop[a1:a2+1]
ctop1 = scldtop1[a1:a2+1]
slwp = radlwp[c1:c2+1]
print "cbase and ctop shapes:", cbase.shape, cbase1.shape, ctop.shape, slwp.shape, '----', idate, t1, t2
mwrdec = MWRdecoupled[m1:m2+1]
mwrtim = MWRtim[m1:m2+1]
BL0 = np.where(mwrdec==0)[0]
BL1 = np.where(mwrdec==1)[0]
print "Lengths of coupled(well-mixed) vs decoupled BLs:", jdate, len(BL0), len(BL1), '--perc--', 100*(len(BL0)/len(mwrdec)), 100*(len(BL1)/len(mwrdec))

##Bifurcation of coupled vs decoupled layers........
#well mixed vs decoupled boundary layer...temperature and pot.temp...height wise....
MTMPBL0mean = np.zeros(len(altitude)); MTMPBL0mean[:] = np.nan;
MTMPBL0medn = np.zeros(len(altitude)); MTMPBL0medn[:] = np.nan;
MTMPBL1mean = np.zeros(len(altitude)); MTMPBL1mean[:] = np.nan;
MTMPBL1medn = np.zeros(len(altitude)); MTMPBL1medn[:] = np.nan;
MPOTBL0mean = np.zeros(len(altitude)); MPOTBL0mean[:] = np.nan;
MPOTBL0medn = np.zeros(len(altitude)); MPOTBL0medn[:] = np.nan;
MPOTBL1mean = np.zeros(len(altitude)); MPOTBL1mean[:] = np.nan;
MPOTBL1medn = np.zeros(len(altitude)); MPOTBL1medn[:] = np.nan;
for r in range(len(altitude)):
  mtmp0 = MWRtemp[:,m1:m2+1][r,BL0]-273.15;  mtmp1 = MWRtemp[:,m1:m2+1][r,BL1]-273.15;
  mpot0 = MWRpot[:,m1:m2+1][r,BL0]-273.15;  mpot1 = MWRpot[:,m1:m2+1][r,BL1]-273.15;
  #tmp
  MTMPBL0mean[r] = mtmp0[~np.isnan(mtmp0)].mean()
  MTMPBL0medn[r] = np.median(mtmp0[~np.isnan(mtmp0)])
  MTMPBL1mean[r] = mtmp1[~np.isnan(mtmp1)].mean()
  MTMPBL1medn[r] = np.median(mtmp1[~np.isnan(mtmp1)])
  #pot
  MPOTBL0mean[r] = mpot0[~np.isnan(mpot0)].mean()
  MPOTBL0medn[r] = np.median(mpot0[~np.isnan(mpot0)])
  MPOTBL1mean[r] = mpot1[~np.isnan(mpot1)].mean()
  MPOTBL1medn[r] = np.median(mpot1[~np.isnan(mpot1)])
  del mtmp0; del mtmp1; del mpot0; del mpot1;

print "MTMP and MPOT of coupled vs decoupled BL's:", len(MTMPBL0mean), len(MTMPBL1mean), len(MPOTBL0mean), len(MPOTBL1mean)
#

#Case-1;; 20180917
#q1
if idate=='20180917':
  cbase[:30] = cbase1[:30]
  cbase[0+30:2+30] = cbase[~np.isnan(cbase)].mean()
  cbase[289+30:291+30] = cbase[~np.isnan(cbase)].mean()
  cbase[318+30:326+30] = cbase[~np.isnan(cbase)].mean()
#q2
if idate=='20180917':
  scldbase[720:780] = scldbase1[720:780]
  scldbase[780:783] = cbase[~np.isnan(cbase)].mean()
  scldbase[1269:1278] = cbase[~np.isnan(cbase)].mean()
  scldbase[1356:1362] = cbase[~np.isnan(cbase)].mean()
  scldbase[1438:1440] = cbase[~np.isnan(cbase)].mean()
  scldbase[1417] = scldbase[1418]
  scldbase[1419] = scldbase[1420]
  scldbase[1421] = scldbase[1422]
  scldbase[1423] = scldbase[1424]
  scldbase[1425] = scldbase[1426]
  scldbase[1427] = scldbase[1428]
  scldbase[1429] = scldbase[1430]
  scldtop1[scldtop1<0.1]=np.nan
  scldtop1[a1:a2+1] = ctop

#q3
if idate=='20180917':
 m,n = REPSerr.shape
 for i in range(m):
  for j in range(n):
    if ((np.isnan(REPSerr[i,j])==False) and (REPSerr[i,j]>300)):
       REPSerr[i,j]=REPSerr[i,j]-75
    if ((np.isnan(REPSerr[i,j])==False) and (REPSerr[i,j]>250)):
       REPSerr[i,j]=REPSerr[i,j]-50
    if ((np.isnan(REPSerr[i,j])==False) and (REPSerr[i,j]>200)):
       REPSerr[i,j]=REPSerr[i,j]-50

 del i; del j; del m; del n;
 #rerr = REPSerr[520:650,:][:,:58]
 for k in np.arange(520,650):
  err = REPSerr[k,:58]
  val = err[~np.isnan(err)].mean()
  inx = np.where(err>val)[0]
  err[inx] = val;
  REPSerr[k,:58] = err;
  del err; del val; del inx;

#q4
#interpolate scldbase and scldtop1...............i.e. cbh and cth....20180917
if idate=='20180917':
  scldtop1[715] = (scldtop1[714]+scldtop1[716])/2
  scldtop1[717] = (scldtop1[716]+scldtop1[718])/2
  #cbase
  abc =  scldbase[720:780]  #12-13UTC
  val = abc[~np.isnan(abc)].mean()
  inx = np.where(np.isnan(abc))[0]
  abc[inx] = val;
  scldbase[720:780] = abc;
  del inx; del abc; del val;
  abc =  scldbase[780:840]  #13-14UTC
  val = abc[~np.isnan(abc)].mean()
  inx = np.where(np.isnan(abc))[0]
  abc[inx] = val;
  scldbase[780:840] = abc;
  del inx; del abc; del val;
  abc =  scldbase[1320:1380]  #22-23UTC
  val = abc[~np.isnan(abc)].mean()
  inx = np.where(np.isnan(abc))[0]
  abc[inx] = val;
  scldbase[1320:1380] = abc;
  del inx; del abc; del val;
  abc =  scldbase[1380:1440]  #23-24UTC
  val = abc[~np.isnan(abc)].mean()
  inx = np.where(np.isnan(abc))[0]
  abc[inx] = val;
  scldbase[1380:1440] = abc;
  del inx; del abc; del val;

#Case2...20180911...quality
#q1
if idate=='20180911':
  i1 = np.where(~np.isnan(cbase1))[0][0]
  i2 = np.where(~np.isnan(cbase1))[0][-1]
  cbase[i1:i2+1] = cbase1[i1:i2+1]
  cbase[cbase<0.1]=np.nan

#q2..........
if idate=='20180911':
  i1=np.where(ctime>6)[0][0]
  i2=np.where(ctime<10)[0][-1]
  scldbase[i1:i2+1]=np.nan
  scldbase[i1:i2+1] = scldbase1[i1:i2+1]
  del scldtop1; scldtop1 = scldtop;
  scldtop1[224] = (scldtop1[223] + scldtop1[225])/2
  scldbase[0:2] =  scldbase[3]
  del ctop; ctop = scldtop1[a1:a2+1]

#q3
#quality-casestudies.....cmap..contourf
if idate=='20180911':
  top1 = scldtop1[330:360]
  scldtop1[333:336] = top1[~np.isnan(top1)].mean()
  scldtop1[349] = top1[~np.isnan(top1)].mean()
  del top1;
  scldtop1[328:338] = scldtop1[328:338]+0.1
  scldtop1[328:333] = scldtop1[328:333]+0.05
  #6-7UTC
  xp = np.array([353,354,355,356,357,372,373,374,375]); fp = scldbase[xp]
  nxp = np.arange(358,372);   nfp = np.interp(nxp, xp, fp)
  scldbase[nxp] = nfp; del nxp; del nfp;
  nxp = np.arange(376,386); nfp = np.interp(nxp, xp, fp);
  scldbase[nxp] = nfp;
  del xp; del fp; del nxp; del nfp;
  #7-8UTC
  xp = np.array([435,436,437,438,439,449,450,451,452,453]); fp = scldbase[xp]
  nxp = np.arange(440,449);   nfp = np.interp(nxp, xp, fp)
  scldbase[nxp] = nfp;
  del xp; del fp; del nxp; del nfp;
  #8-9UTC
  xp = np.array([498,499,500,501,502,516,517,518,519,520]); fp = scldbase[xp]
  nxp = np.arange(503,516);   nfp = np.interp(nxp, xp, fp)
  scldbase[nxp] = nfp;
  del xp; del fp; del nxp; del nfp;
  xp = np.array([519,520,521,522,523,539,540,541,542,543]); fp = scldbase[xp]
  nxp = np.arange(524,539);   nfp = np.interp(nxp, xp, fp)
  scldbase[nxp] = nfp;
  del xp; del fp; del nxp; del nfp;
  #9-10UTC
  xp = np.array([540,541,542,543,544,545,546,547,548,611,612,613,614,615,616,617,618,619])
  fp = scldbase[xp]
  #nxp = np.arange(549,611); nfp = np.interp(nxp, xp, fp)
  nxp = np.arange(585,611); nfp = np.interp(nxp, xp, fp)
  scldbase[nxp] = nfp;
  del xp; del fp; del nxp; del nfp;

#Case3---20180913...0to14UTC
#q1
if idate=='20180913':
  cbase[cbase<0.1]=np.nan
  cix = np.where(cbase>1.2)[0]
  cbase[cix]=np.mean(cbase)
  ctop = ctop1;
  ctop[ctop>2]=np.nan
  ctop[246] = (ctop[245] + ctop[247])/2
  ctop[255] = (ctop[254] + ctop[256])/2

#q2
if idate=='20180913':
  scldbase[807] = scldbase[808]
  scldbase[809] = scldbase[808]
  scldbase[811] = scldbase[808]
  scldbase[813] = (scldbase[812] + scldbase[814])/2
  scldtop1[246] = (scldtop1[245] + scldtop1[247])/2
  scldtop1[255] = (scldtop1[254] + scldtop1[256])/2

#Case-4...20180831...0to4 UTC...
#q1
if idate=='20180831':
  cbase[0:10] = cbase[~np.isnan(cbase)].mean()
  cbase[30] = cbase[~np.isnan(cbase)].mean()
  cbase[37:39] = cbase[~np.isnan(cbase)].mean()
#q2
if idate=='20180831':
  scldbase[0:20] = cbase[~np.isnan(cbase)].mean()
  scldbase[60:62] = cbase[~np.isnan(cbase)].mean()
  scldbase[72:78] =  cbase[~np.isnan(cbase)].mean()

#Case-5...20180823...12to24 UTC...
#q1
if idate=='20180823':
  cbase = cbase1;
  base = cbase[:100]  #considering 12-15.5UTC time..
  cldtim = atime[a1:a2+1]
  nbase = cbase[270:];
  base1 = np.zeros(len(base)); base1[:] = np.nan
  for i in range(len(base)):
    if ((np.isnan(base[i])==True) or (base[i]<1)):
      base1[i] = nbase[~np.isnan(nbase)].mean()
    else:
      base1[i] = base[i]
  cbase[:100] = base1;
  del base; del base1; del i; del nbase;
  inx = np.where(ctop<0.5)[0]
  ctop[inx] = ctop[~np.isnan(ctop)].mean()
  del inx;

#q2
if idate=='20180823':
  n2 = np.where(mtime>=atime[525])[0][0]
  MWRtemp[:, m1:n2] = np.nan;
  mwrdec1 = MWRdecoupled[n2:]
  BLs0 = np.where(mwrdec1==0)[0]
  BLs1 = np.where(mwrdec1==1)[0]
  print "Lengths of coupled(well-mixed) vs decoupled BLs:", jdate, len(BLs0), len(BLs1), '--perc--', 100*(len(BLs0)/len(mwrdec1)), 100*(len(BLs1)/len(mwrdec1))

#q3
if idate=='20180823':
  del scldbase; scldbase = scldbase1; #scldbase[a1:a2+1] = scldbase1[a1:a2+1]
  base = scldbase[c1:c2+1][:200]
  base1 = np.zeros(len(base)); base1[:] = np.nan
  nbase = scldbase[1260:]
  for i in range(len(base)):
    if ((np.isnan(base[i])==True) or (base[i]<1)):
      base1[i] = nbase[~np.isnan(nbase)].mean()
    else:
      base1[i] = base[i]
  scldbase[c1:c2+1][:200] = smooth(base1,5);
  scldbase[919] = cbase[~np.isnan(cbase)].mean()
  del base; del base1; del i; del nbase;
  scldtop1[a1:a2+1] = ctop

#Case-6...20180903....19to22 UTC.........
#q1
if idate=='20180903':
  rdop = RDOP[a1:a2+1,:][:,:58]
  cbase2 = np.zeros(len(cbase)); cbase2[:] = np.nan;
  for i in range(len(cbase)):
    if ((np.isnan(cbase1[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<50)):
      cbase2[i] = cbase1[~np.isnan(cbase1)].mean()
    else:
      cbase2[i] = cbase1[i]
  cbase = cbase2; del cbase2; del i; del rdop;
  ctop[ctop>2]=np.nan;
#q2
if idate=='20180903':
  base = scldbase1[c1:c2+1]
  base1 = np.zeros(len(base)); base1[:] = np.nan
  for i in range(len(base)):
    if np.isnan(base[i])==True:
      base1[i] = cbase[~np.isnan(cbase)].mean()
    else:
      base1[i] = base[i]
  inx =  np.where(np.isnan(cbase))[0]
  base1[inx*2] = np.nan; base1[inx*2+1] = np.nan;
  scldbase[c1:c2+1] = base1; del base1;  del i; del inx;
  scldtop1 = scldtop;

#Case-7....20180905...5to15 UTC....
#q1
#quality of cld base and cld top heights....for analyzing norm heights
if idate=='20180905':
  ctop[ctop>2.6]=np.nan;
  rdop = RDOP[a1:a2+1,:][:,:100]
  ctop2 = np.zeros(len(cbase)); ctop2[:] = np.nan;
  cbase2 = np.zeros(len(cbase)); cbase2[:] = np.nan;
  for i in range(len(cbase)):
    if ((np.isnan(ctop[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<=80)):
      ix = np.where(~np.isnan(rdop[i]))[0][-1]
      ctop2[i] = Rranges[ix]; del ix;
    else:
      ctop2[i] = ctop[i]
  ctop = ctop2; del ctop2; del i;
  #1
  for i in range(len(cbase)):
    if ((np.isnan(cbase1[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<80)):
      ix = np.where(~np.isnan(rdop[i]))[0][-1]
      rdop1 = rdop[i,:][:ix][::-1]
      ix1 = np.where(~np.isnan(rdop1))[0][-1]
      cbase2[i] = Rranges[ix-ix1]
      del ix; del ix1; del rdop1;
    else:
      cbase2[i] = cbase1[i]
  cbase = cbase2; del cbase2; del i;
  #2
  cbase[cbase<1]=np.nan;
  cbase2 = np.zeros(len(cbase)); cbase2[:] = np.nan;
  for i in range(len(cbase)):
    if ((np.isnan(cbase[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<80)):
      cbase2[i] = cbase[~np.isnan(cbase)].mean()
    else:
      cbase2[i] = cbase[i]
  cbase = cbase2; del cbase2; del i;
  #3
  base = cbase[:30]; #5-6 UTC
  cbase[12:19] = base[~np.isnan(base)].mean() + 0.2
  cbase[24:30] = base[~np.isnan(base)].mean() + 0.2
  del base;
  base = cbase[30:60]; #6-7 UTC
  cbase[30:34] = base[~np.isnan(base)].mean()
  cbase[37] = base[~np.isnan(base)].mean()
  cbase[41] = base[~np.isnan(base)].mean()
  del base;
  base = cbase[60:90] #7-8UTC
  cbase[83:87] = base[~np.isnan(base)].mean()
  cbase[88] = base[~np.isnan(base)].mean()
  del base;
  base = cbase[120:150] #9-10UTC
  cbase[117:124] = base[~np.isnan(base)].mean()
  cbase[144:147] = base[~np.isnan(base)].mean()
  cbase[194] = base[~np.isnan(base)].mean()
  cbase[204] = base[~np.isnan(base)].mean()
  cbase[213:215] = base[~np.isnan(base)].mean()
  cbase[220:222] = base[~np.isnan(base)].mean()
  cbase[224:226] = base[~np.isnan(base)].mean()
  cbase[229:231] = base[~np.isnan(base)].mean()
  cbase[270] = base[~np.isnan(base)].mean()
  cbase[282] = base[~np.isnan(base)].mean()
  cbase[286] = base[~np.isnan(base)].mean()
  del base;

#q2
if idate=='20180905':
  scldtop1[a1:a2+1] = ctop;
  inx = np.where(np.logical_and(atime>=t1, atime<=t2))[0]
  scldbase[inx*2] = cbase;
  scldbase[inx*2+1] = cbase;
  del inx;

##Case-8...20180901...3to18.UTC...
if idate=='20180901':
  scldtop1[a1:a2+1]=ctop;
  scldtop1[scldtop1>1]=np.nan;  

#Case-9...20180902...0to15.UTC...
#q1
if idate=='20180902':
  cldtim = atime[a1:a2+1]
  ctop1[ctop1>2]=np.nan;
  top = ctop1[420:450]  #14-15UTC
  ctop1[441:445] = top[~np.isnan(top)].mean()
  ctop1[433:438] = top[~np.isnan(top)].mean()
  del top;
  top = ctop[390:420]  #13-14UTC
  ctop1[393] = top[~np.isnan(top)].mean()
  ctop1[411] = top[~np.isnan(top)].mean()
  ctop1[414:421] = top[~np.isnan(top)].mean()
  del top;
  top = ctop1[330:360]  #11-12UTC
  ctop1[331:333] = top[~np.isnan(top)].mean()
  ctop1[327:329] = top[~np.isnan(top)].mean()
  del top;
  ctop1[:181] = ctop[:181];
  del ctop; ctop = ctop1; ctop[ctop>2]=np.nan;
  top = ctop[:30]  ##00-01UTC
  ctop[11] = top[~np.isnan(top)].mean()
  ctop[15:17] = top[~np.isnan(top)].mean()
  ctop[22:26] = top[~np.isnan(top)].mean()
  ctop[27:29] = top[~np.isnan(top)].mean()
  del top;

#Case-10...20180906...14to22.UTC...
#q1
if idate=='20180906':
  base = cbase[30:60]  #15to16UTC.
  cbase[30:37]=base[~np.isnan(base)].mean()
  cbase[40:47]=base[~np.isnan(base)].mean()
  del base;
  base = cbase[90:120]  #17to18UTC.
  cbase[90:94]=base[~np.isnan(base)].mean()
  cbase[116:122]=base[~np.isnan(base)].mean()
  del base;
  
#q2
if idate=='20180906':
  scldtop1[a1:a2+1]=ctop;

#Case-11...20180908..5to10.UTC..
#q1
if idate == '20180908':
  cbase = cbase1;
  ctop[ctop>1.25]=np.nan;
  ctop2 = np.zeros(len(cbase)); ctop2[:] = np.nan;
  cbase2 = np.zeros(len(cbase)); cbase2[:] = np.nan;
  rdop = RDOP[a1:a2+1,:][:,:58]
  for i in range(len(cbase)):
    if ((np.isnan(cbase[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<=50)):
      cbase2[i] = cbase[~np.isnan(cbase)].mean()
    else:
      cbase2[i] = cbase[i]

    if ((np.isnan(ctop[i])==True) and (len(np.where(np.isnan(rdop[i,:]))[0])<=50)):
      ctop2[i] = ctop[~np.isnan(ctop)].mean()
    else:
      ctop2[i] = ctop[i]

  cbase = cbase2; ctop = ctop2; del cbase2; del ctop2; del i; del rdop;
#q2
if idate == '20180908':
  scldbase = scldbase1;
  scldtop1[a1:a2+1] = ctop;



##MEAN plot of lidar eps and radar eps with mdecoupled...overlay.............
outdir1 = os.getcwd() + '/cases/'
##Plotting the radar attributes.................
#Figure of case study only radar attribtues....
dt=2/60
#colorbar....
eps_levs = [10**-7, 10**-6/2, 10**-6, 10**-5/2, 10**-5, 10**-4/2, 10**-4, 10**-3/2, 10**-3, 10**-2/2, 10**-2, 10**-1/2, 10**-1]
#eps_cmap = matplotlib.colors.ListedColormap(['black','rosybrown','darkred','red','orange','gold','yellow','forestgreen','limegreen','lime','skyblue','cyan','gray'])
eps_cmap = matplotlib.cm.get_cmap('jet')
eps_norm = matplotlib.colors.BoundaryNorm(eps_levs,eps_cmap.N,clip=True)
#
scldbase = smooth(scldbase, 5)
if idate=="20180905":
  cldbase[cldbase>1500]=np.nan

#24-h picture of epsilon, doppler, temperature profile.....................
plt.clf(); plt.close()
fig = plt.figure(figsize=(11,6)); fig.subplots_adjust(wspace=0.15, hspace=0.3)
fig.suptitle('LIDARvsRADAR: Dissipation rates:' + idate ,x=0.2, y=0.965,horizontalalignment='left',verticalalignment='top', fontsize=12)
mdecoupled1 = MWRdecoupled; mdecoupled1[mdecoupled1==0]=np.nan;
ax3 = fig.add_subplot(211)
cf3 = ax3.contourf(atime, Rranges, 10**REPS.T, cmap=eps_cmap, levels=eps_levs, norm=eps_norm)
#ax3.plot(ctime, scldbase, linewidth=1.0, color='k', label='CBH-C')
ax3.plot(ctime, cldbase/1000, linewidth=1.0, color='k', label='CBH-C')
ax3.plot(atime, scldtop1, linewidth=1.0, color='r', label='CTH-R')
#ax3.plot(atime[a1:a2+1], cbase, linewidth=1.0, color='k', label='CBH-R')
ax3.plot(MWRtim, mdecoupled1+0.5, linewidth=3, color='m', marker='*', markeredgecolor='m')
ax3.plot(MWRtim, MWRdecbase/1000, linewidth=2.0, color='r')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(a) RADAR: log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=12)
plt.ylabel('Range (km)',size=12)
plt.colorbar(cf3,shrink=0.95, format='%1.0e')

ax3 = fig.add_subplot(212)
cf3 = ax3.contourf(atime, Lranges, 10**LEPS.T, cmap=eps_cmap, levels=eps_levs, norm=eps_norm)
#ax3.plot(ctime, scldbase, linewidth=1.0, color='k', label='CBH-C')
ax3.plot(ctime, cldbase/1000, linewidth=1.0, color='k', label='CBH-C')
ax3.plot(atime, scldtop1, linewidth=1.0, color='r', label='CTH-R')
#ax3.plot(atime[a1:a2+1], cbase, linewidth=1.0, color='k', label='CBH-R')
ax3.plot(MWRtim, mdecoupled1+0.5, linewidth=3, color='m', marker='*', markeredgecolor='m')
ax3.plot(MWRtim, MWRdecbase/1000, linewidth=2.0, color='r')
#ax3.plot(MWRtim, MWRinvbase/1000, linewidth=2, color='g', linestyle='--', label='invbase')
#ax3.plot(MWRtim, MWRinvtop/1000, linewidth=2, color='g', linestyle='--', label='invtop')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(b) LIDAR: log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=12)
plt.ylabel('Range (km)',size=12); plt.xlabel('Time [UTC]',size=14);
plt.colorbar(cf3,shrink=0.95, format='%1.0e')

plt.savefig(outdir1 + idate + 'LIDARvsRADAReps_casespl.png', format='png', dpi=100)
del fig; del ax3;
#sys.exit()


##
dop_levs = np.arange(-1.5,1.51,0.2)
dop_cmap = matplotlib.cm.get_cmap('jet')
dop_norm = matplotlib.colors.BoundaryNorm(dop_levs,dop_cmap.N,clip=False)
temp_levs = np.arange(-15,5.1,0.5);
temp_cmap = matplotlib.cm.get_cmap('jet');
temp_norm = matplotlib.colors.BoundaryNorm(temp_levs,temp_cmap.N,clip=True)
pot_levs = np.arange(-10,10.1,0.5);
pot_cmap = matplotlib.cm.get_cmap('jet');
pot_norm = matplotlib.colors.BoundaryNorm(pot_levs,pot_cmap.N,clip=True)

#24-h picture of epsilon, doppler, temperature profile.....................
plt.clf(); plt.close()
fig = plt.figure(figsize=(7,12)); fig.subplots_adjust(wspace=0.15, hspace=0.285)
fig.suptitle('LIDAR vs RADAR: Dissipation rate & Doppler velocity :' + idate ,x=0.1, y=0.955,horizontalalignment='left',verticalalignment='top', fontsize=12)
ax1 = fig.add_subplot(511)
cf1 = ax1.contourf(atime, Lranges, 10**LEPS.T, cmap=eps_cmap, levels=eps_levs, norm=eps_norm)
ax1.plot(ctime, cldbase/1000, linewidth=1.5, color='k', label='CBH-C')
#ax1.plot(ctime, scldbase1, linewidth=1.5, color='b', label='CBH-C')
#ax1.plot(atime, scldtop, linewidth=1.5, color='g', label='CTH-R')
ax1.plot(atime, scldtop1, linewidth=1.5, color='r', label='CTH-R')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(a) LIDAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=12)
plt.ylabel('Range (km)',size=14)
cb=plt.colorbar(cf1,shrink=0.95, format='%1.0e')
#plt.legend(loc='center left',borderpad=0.2, bbox_to_anchor=(0.2, 1.22), ncol=3, prop={'size': 12});

ax2 = fig.add_subplot(512)
cf2 = ax2.contourf(atime, Rranges, 10**REPS.T, cmap=eps_cmap, levels=eps_levs, norm=eps_norm)
ax2.plot(ctime, cldbase/1000, linewidth=1.5, color='k', label='CBH-C')
#ax2.plot(atime, scldtop, linewidth=1.5, color='g', label='CTH-R')
ax2.plot(atime, scldtop1, linewidth=1.5, color='r', label='CTH-R')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(b) RADAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=12)
plt.ylabel('Range (km)',size=14)
cb=plt.colorbar(cf2,shrink=0.95, format='%1.0e')

ax3 = fig.add_subplot(513)
cf3 = ax3.contourf(atime, Lranges, LDOP.T, cmap=dop_cmap, levels=dop_levs, norm=dop_norm, extend='both')
ax3.plot(ctime, cldbase/1000, linewidth=1.5, color='k', label='CBH-C')
#ax3.plot(ctime, scldbase1, linewidth=1.5, color='b', label='CBH-C')
#ax3.plot(atime, scldtop, linewidth=1.5, color='g', label='CTH-R')
ax3.plot(atime, scldtop1, linewidth=1.5, color='r', label='CTH-R')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(c) LIDAR Doppler velocity [m s$\mathregular{^{-1}}$]',size=12)
plt.ylabel('Range (km)',size=14)
plt.colorbar(cf3,shrink=0.95,extend='both')
cf3.cmap.set_under('gray');   cf3.cmap.set_over('black')

ax4 = fig.add_subplot(514)
cf4 = ax4.contourf(atime, Rranges, RDOP.T, cmap=dop_cmap, levels=dop_levs, norm=dop_norm, extend='max')
ax4.plot(ctime, cldbase/1000, linewidth=1.5, color='k', label='CBH-C')
#ax4.plot(ctime, scldbase1, linewidth=1.5, color='b', label='CBH-C')
#ax4.plot(atime, scldtop, linewidth=1.5, color='g', label='CTH-R')
ax4.plot(atime, scldtop1, linewidth=1.5, color='r', label='CTH-R')
plt.xlim(t1,t2); plt.ylim(h1,h2); plt.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(d) RADAR: Doppler velocity [m s$\mathregular{^{-1}}$]',size=12)
plt.ylabel('Range (km)',size=14)
plt.colorbar(cf4,shrink=0.95,extend='both')
cf4.cmap.set_under('gray');   cf4.cmap.set_over('black')

mdecoupled1 = MWRdecoupled; mdecoupled1[mdecoupled1==1]=np.nan;
ax6 = fig.add_subplot(515)
cf6 = ax6.contourf(MWRtim, Mranges, MWRtemp-273, cmap=temp_cmap, levels=temp_levs, norm=temp_norm, extend='both')
#ax6.plot(MWRtim, mdecoupled1+1.5, linewidth=3, color='m', marker='*', markeredgecolor='m')
ax6.plot(ctime, cldbase/1000, linewidth=1.5, color='k', label='CBH-C')
#ax6.plot(ctime, scldbase1, linewidth=1.5, color='b', label='CBH-C')
#ax6.plot(atime, scldtop, linewidth=1.5, color='g', label='CTH-R')
ax6.plot(atime, scldtop1, linewidth=1.5, color='r', label='CTH-R')
#ax6.plot(MWRtim, MWRinvbase/1000, linewidth=2, color='k', linestyle='--', label='invbase')
#ax6.plot(MWRtim, MWRinvtop/1000, linewidth=2, color='m', linestyle='--', label='invtop')
plt.xlim(t1,t2); plt.ylim(h1,h2); ax6.grid(True)
plt.yticks(np.arange(h1,h2+0.1,hh)); plt.xticks(np.arange(t1,t2+0.1,1))
plt.title('(e) MWR: Temperature Profile [degC]',size=12)
plt.ylabel('Altitude [km]',size=14);
plt.xlabel('Time [UTC]',size=14);
plt.colorbar(cf6, shrink=0.85)
cf6.cmap.set_under('gray');   cf6.cmap.set_over('black')
del ax6; del cf6;
del ax1; del ax2; del ax3; del ax4; del cf1; del cf2; del cf3; del cf4; del cb;
plt.savefig(outdir1 + idate + 'LIDARvsRADAR_Cloud_epsdopT.png', format='png', dpi=100)
del fig;

toc()
#sys.exit()


###Post-processing.algorithm starts here....
#'''
##Retriveing the profiles of sub-cloud and cloud....layers.........then combining.......both.
#cloud base to cloud top -lidar -norm
cldlen = np.arange(a1,a2+1,1)
sfcb_tmpM = np.zeros((len(cldlen),10)); sfcb_tmpM[:]=np.nan; sfcb_potM = np.zeros((len(cldlen),10)); sfcb_potM[:]=np.nan;
sfcb_wndM = np.zeros((len(cldlen),10)); sfcb_wndM[:]=np.nan;
sfcb_epsL = np.zeros((len(cldlen),10)); sfcb_epsL[:]=np.nan; sfcb_epsR = np.zeros((len(cldlen),10)); sfcb_epsR[:]=np.nan;
sfcb_dopL = np.zeros((len(cldlen),10)); sfcb_dopL[:]=np.nan; sfcb_dopR = np.zeros((len(cldlen),10)); sfcb_dopR[:]=np.nan;
sfcb_sigL = np.zeros((len(cldlen),10)); sfcb_sigL[:]=np.nan; sfcb_sigR = np.zeros((len(cldlen),10)); sfcb_sigR[:]=np.nan;
sfcb_skwL = np.zeros((len(cldlen),10)); sfcb_skwL[:]=np.nan; sfcb_skwR = np.zeros((len(cldlen),10)); sfcb_skwR[:]=np.nan;
#
for k in range(len(np.arange(a1,a2+1,1))):
  cld1 = cbase[k];   cld2 = ctop[k]
  if np.isnan(cld1)==True:
    cld1 = cbase[~np.isnan(cbase)].mean()
  if np.isnan(cld2)==True:
    cld2 = np.nanmax(ctop)
  if cld2-cld1>1:
    cld1 = np.mean(cbase)
  #lidar
  v1,ix1 = find_nearest(Lranges, cld1);
  v2,ix2 = find_nearest(Lranges, cld2);
  ldop = LDOP[k+a1,:]
  if ((np.isnan(ldop[ix1])==False) and (san==1)):
    ldop1 = ldop[:ix1][::-1]
    ix = np.where(np.isnan(ldop1))[0][0]
    ix1 = ix1-ix
    del ldop1;
  del ldop;
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     ldop = LDOP[k+a1,:][:ix1+1]
     lsig = LSIG[k+a1,:][:ix1+1]
     lskw = LSKW[k+a1,:][:ix1+1]
     leps = LEPS[k+a1,:][:ix1+1]
     lrng = Lranges[:ix1+1]
     #
     #surface to cloud base height.......normalized height.....lidar
     if ix1>=10:
        lht = Lranges[:ix1+1]
        lhtn = lht/max(lht)
        t=0;
        for j in range(10):
          linx = np.where(np.logical_and(lhtn>=t, lhtn<t+0.1))[0]
          eps = leps[linx];  dop = ldop[linx];  sig = lsig[linx];  skw = lskw[linx]
          sfcb_epsL[k,j] = eps[~np.isnan(eps)].mean()
          sfcb_dopL[k,j] = dop[~np.isnan(dop)].mean()
          sfcb_sigL[k,j] = sig[~np.isnan(sig)].mean()
          sfcb_skwL[k,j] = skw[~np.isnan(skw)].mean()
          del eps; del linx; del dop; del sig; del skw;
          t=t+0.1
        del j; del t; del lht; del lhtn;
     del ldop; del lsig; del lskw; del leps; del lrng;
  #radar
  del ix1; del ix2; del v1; del v2;
  v1,ix1 = find_nearest(Rranges, cld1);
  v2,ix2 = find_nearest(Rranges, cld2);
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     rdop = RDOP[k+a1,:][:ix1+1]
     rsig = RSIG[k+a1,:][:ix1+1]
     rskw = RSKW[k+a1,:][:ix1+1]
     reps = REPS[k+a1,:][:ix1+1]
     rrng = Rranges[:ix1+1]
     #
     #surface to cloud base....normalized height.....radar
     if ix1>=10:
        rht = Rranges[:ix1+1]
        rhtn = rht/max(rht)
        t=0;
        for j in range(10):
          rinx = np.where(np.logical_and(rhtn>=t, rhtn<t+0.1))[0]
          eps = reps[rinx];  dop = rdop[rinx];  sig = rsig[rinx];  skw = rskw[rinx]
          sfcb_epsR[k,j] = eps[~np.isnan(eps)].mean()
          sfcb_dopR[k,j] = dop[~np.isnan(dop)].mean()
          sfcb_sigR[k,j] = sig[~np.isnan(sig)].mean()
          sfcb_skwR[k,j] = skw[~np.isnan(skw)].mean()
          del eps; del rinx; del dop; del sig; del skw;
          t=t+0.1
        del j; del t; del rht; del rhtn;
     del rdop; del rsig; del rskw; del reps; del rrng;
     del ix1; del ix2; del v1; del v2;
  #mwr-temp--pot...
  v1,ix1 = find_nearest(Mranges, cld1);
  v2,ix2 = find_nearest(Mranges, cld2);
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     mtmp = MTMP[k+a1,:][:ix1+1]
     mpot = MPOT[k+a1,:][:ix1+1]
     mwnd = MWND[k+a1,:][:ix1+1]
     mrng = Mranges[:ix1+1]
     #
     #cloud base to cloud top....normalized height.....radar
     if ix1>=10:
        mht = Mranges[:ix1+1]
        mhtn = mht/max(mht)
        t=0;
        for j in range(10):
          minx = np.where(np.logical_and(mhtn>=t, mhtn<t+0.1))[0]
          tmp = mtmp[minx]; pot = mpot[minx]; wnd = mwnd[minx]
          sfcb_tmpM[k,j] = tmp[~np.isnan(tmp)].mean()
          sfcb_potM[k,j] = pot[~np.isnan(pot)].mean()
          sfcb_wndM[k,j] = wnd[~np.isnan(wnd)].mean()
          del tmp; del pot; del minx; del wnd;
          t=t+0.1
        del j; del t; del mht; del mhtn;
     del mtmp; del mrng; del mpot; del mwnd;
  del ix1; del ix2; del v1; del v2;
  del cld1; del cld2;
  #sys.exit()

print "----------------------------------SUB-CLOUDS------------------------------------------"
print "Surface to Cloud base height....profile......shapes", sfcb_epsL.shape, sfcb_epsR.shape, sfcb_tmpM.shape, sfcb_potM.shape, sfcb_wndM.shape
print "Surface to Cloud base height: LIDAR, DOP, SIG, SKW ", sfcb_dopL.shape, sfcb_sigL.shape, sfcb_skwL.shape
print "Surface to Cloud base height: RADAR, DOP, SIG, SKW ", sfcb_dopR.shape, sfcb_sigR.shape, sfcb_skwR.shape


#
mtmpdata=np.zeros(2); mtmpdata[:]=np.nan;
ldopdata=np.zeros(2); ldopdata[:]=np.nan;  ldopndata=np.zeros(2); ldopndata[:]=np.nan;  rdopdata=np.zeros(2); rdopdata[:]=np.nan;
lsigdata=np.zeros(2); lsigdata[:]=np.nan;  lsigndata=np.zeros(2); lsigndata[:]=np.nan;  rsigdata=np.zeros(2); rsigdata[:]=np.nan;
lskwdata=np.zeros(2); lskwdata[:]=np.nan;  lskwndata=np.zeros(2); lskwndata[:]=np.nan;  rskwdata=np.zeros(2); rskwdata[:]=np.nan;
lepsdata=np.zeros(2); lepsdata[:]=np.nan;  lepsndata=np.zeros(2); lepsndata[:]=np.nan;  repsdata=np.zeros(2); repsdata[:]=np.nan;
lrngdata=np.zeros(2); lrngdata[:]=np.nan;  rrngdata=np.zeros(2);  rrngdata[:]=np.nan;   mrngdata=np.zeros(2); mrngdata[:]=np.nan;
#cloud base to cloud top -lidar -norm
cldlen = np.arange(a1,a2+1,1)
cbct_tmpM = np.zeros((len(cldlen),10)); cbct_tmpM[:]=np.nan; cbct_potM = np.zeros((len(cldlen),10)); cbct_potM[:]=np.nan;
cbct_wndM = np.zeros((len(cldlen),10)); cbct_wndM[:]=np.nan;
cbct_epsL = np.zeros((len(cldlen),10)); cbct_epsL[:]=np.nan; cbct_epsR = np.zeros((len(cldlen),10)); cbct_epsR[:]=np.nan;
cbct_dopL = np.zeros((len(cldlen),10)); cbct_dopL[:]=np.nan; cbct_dopR = np.zeros((len(cldlen),10)); cbct_dopR[:]=np.nan;
cbct_sigL = np.zeros((len(cldlen),10)); cbct_sigL[:]=np.nan; cbct_sigR = np.zeros((len(cldlen),10)); cbct_sigR[:]=np.nan;
cbct_skwL = np.zeros((len(cldlen),10)); cbct_skwL[:]=np.nan; cbct_skwR = np.zeros((len(cldlen),10)); cbct_skwR[:]=np.nan;
#
for k in range(len(np.arange(a1,a2+1,1))):
  cld1 = cbase[k];   cld2 = ctop[k]
  if np.isnan(cld1)==True:
    cld1 = cbase[~np.isnan(cbase)].mean()
  if np.isnan(cld2)==True:
    cld2 = np.nanmax(ctop)
  if cld2-cld1>1:
    cld1 = np.mean(cbase)
  #lidar
  v1,ix1 = find_nearest(Lranges, cld1);
  v2,ix2 = find_nearest(Lranges, cld2);
  ldop = LDOP[k+a1,:]
  if ((np.isnan(ldop[ix1])==False) and (san==1)):
    ldop1 = ldop[:ix1][::-1]
    ix = np.where(np.isnan(ldop1))[0][0]
    ix1 = ix1-ix
    del ldop1;
  del ldop;
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     ldop = LDOP[k+a1,:][ix1:ix2+1]
     lsig = LSIG[k+a1,:][ix1:ix2+1]
     lskw = LSKW[k+a1,:][ix1:ix2+1]
     leps = LEPS[k+a1,:][ix1:ix2+1]
     lrng = Lranges[ix1:ix2+1]
     if k==0:
       ldopdata=ldop;  lsigdata=lsig;  lskwdata=lskw;  lepsdata=leps;  lrngdata=lrng;
     else:
       ldopdata = np.concatenate([ldopdata,ldop])
       lsigdata = np.concatenate([lsigdata,lsig])
       lskwdata = np.concatenate([lskwdata,lskw])
       lepsdata = np.concatenate([lepsdata,leps])
       lrngdata = np.concatenate([lrngdata,lrng])
     #
     #cloud base to cloud top....normalized height.....lidar
     if ix2-ix1>=10:
        lht = Lranges[ix1:ix2+1]
        lhtn = (lht-min(lht))/max(lht-min(lht))
        t=0;
        for j in range(10):
          linx = np.where(np.logical_and(lhtn>=t, lhtn<t+0.1))[0]
          eps = leps[linx];  dop = ldop[linx];  sig = lsig[linx];  skw = lskw[linx]
          cbct_epsL[k,j] = eps[~np.isnan(eps)].mean()
          cbct_dopL[k,j] = dop[~np.isnan(dop)].mean()
          cbct_sigL[k,j] = sig[~np.isnan(sig)].mean()
          cbct_skwL[k,j] = skw[~np.isnan(skw)].mean()
          del eps; del linx; del dop; del sig; del skw;
          t=t+0.1
        del j; del t; del lht; del lhtn;
     del ldop; del lsig; del lskw; del leps; del lrng;
  #radar
  del ix1; del ix2; del v1; del v2;
  v1,ix1 = find_nearest(Rranges, cld1);
  v2,ix2 = find_nearest(Rranges, cld2);
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     rdop = RDOP[k+a1,:][ix1:ix2+1]
     rsig = RSIG[k+a1,:][ix1:ix2+1]
     rskw = RSKW[k+a1,:][ix1:ix2+1]
     reps = REPS[k+a1,:][ix1:ix2+1]
     rrng = Rranges[ix1:ix2+1]
     if k==0:
       rdopdata=rdop;  rsigdata=rsig;  rskwdata=rskw;  repsdata=reps;  rrngdata=rrng;
     else:
       rdopdata = np.concatenate([rdopdata,rdop])
       rsigdata = np.concatenate([rsigdata,rsig])
       rskwdata = np.concatenate([rskwdata,rskw])
       repsdata = np.concatenate([repsdata,reps])
       rrngdata = np.concatenate([rrngdata,rrng])
     #
     #cloud base to cloud top....normalized height.....radar
     if ix2-ix1>=10:
        rht = Rranges[ix1:ix2+1]
        rhtn = (rht-min(rht))/max(rht-min(rht))
        t=0;
        for j in range(10):
          rinx = np.where(np.logical_and(rhtn>=t, rhtn<t+0.1))[0]
          eps = reps[rinx];  dop = rdop[rinx];  sig = rsig[rinx];  skw = rskw[rinx]
          cbct_epsR[k,j] = eps[~np.isnan(eps)].mean()
          cbct_dopR[k,j] = dop[~np.isnan(dop)].mean()
          cbct_sigR[k,j] = sig[~np.isnan(sig)].mean()
          cbct_skwR[k,j] = skw[~np.isnan(skw)].mean()
          del eps; del rinx; del dop; del sig; del skw;
          t=t+0.1
        del j; del t; del rht; del rhtn;
     del rdop; del rsig; del rskw; del reps; del rrng;
  del ix1; del ix2; del v1; del v2;
  #mwr-temp--pot
  v1,ix1 = find_nearest(Mranges, cld1);
  v2,ix2 = find_nearest(Mranges, cld2);
  if ((cld2>cld1) and (ix2-ix1>5) and (cld2-cld1<1)):
     mtmp = MTMP[k+a1,:][ix1:ix2+1]
     mpot = MPOT[k+a1,:][ix1:ix2+1]
     mwnd = MWND[k+a1,:][ix1:ix2+1]
     mrng = Mranges[ix1:ix2+1]
     #
     #cloud base to cloud top....normalized height.....radar
     if ix2-ix1>=10:
        mht = Mranges[ix1:ix2+1]
        mhtn = (mht-min(mht))/max(mht-min(mht))
        t=0;
        for j in range(10):
          minx = np.where(np.logical_and(mhtn>=t, mhtn<t+0.1))[0]
          tmp = mtmp[minx];  pot = mpot[minx]; wnd = mwnd[minx]
          cbct_tmpM[k,j] = tmp[~np.isnan(tmp)].mean()
          cbct_potM[k,j] = pot[~np.isnan(pot)].mean()
          cbct_wndM[k,j] = wnd[~np.isnan(wnd)].mean()
          del tmp; del minx; del pot; del wnd;
          t=t+0.1
        del j; del t; del mht; del mhtn;
     del mtmp; del mrng; del mpot; del mwnd;
  del ix1; del ix2; del v1; del v2;
  del cld1; del cld2;
  #sys.exit()

print "----------------------------------CLOUDS------------------------------------------"
print "Cloud base to Cloud top height....profile......shapes", cbct_epsL.shape, cbct_epsR.shape, cbct_tmpM.shape, cbct_potM.shape, cbct_wndM.shape
print "Cloud base to Cloud top height: LIDAR, DOP, SIG, SKW ", cbct_dopL.shape, cbct_sigL.shape, cbct_skwL.shape
print "Cloud base to Cloud top height: RADAR, DOP, SIG, SKW ", cbct_dopR.shape, cbct_sigR.shape, cbct_skwR.shape
print "Cloud: lidar cld base to cld top shapes:", ldopdata.shape, lsigdata.shape, lskwdata.shape, lepsdata.shape
print "Cloud: radar cld base to cld top shapes:", rdopdata.shape, rsigdata.shape, rskwdata.shape, repsdata.shape

#sys.exit()
##sfcb and cbct.....seggregate the data......12-15UTC, 15-20.5UTC, 20.5-24UTC....2categoreis....
SubCloudmean_epsL = np.zeros(10);  SubCloudmedn_epsL = np.zeros(10);
SubCloudmean_epsR = np.zeros(10);  SubCloudmedn_epsR = np.zeros(10);
Cloudmean_epsL = np.zeros(10);  Cloudmedn_epsL = np.zeros(10);
Cloudmean_epsR = np.zeros(10);  Cloudmedn_epsR = np.zeros(10);
SubCloudmean_epsL1 = np.zeros(10); SubCloudmedn_epsL1 = np.zeros(10);
SubCloudmean_epsL2 = np.zeros(10); SubCloudmedn_epsL2 = np.zeros(10);
SubCloudmean_epsR1 = np.zeros(10); SubCloudmedn_epsR1 = np.zeros(10);
SubCloudmean_epsR2 = np.zeros(10); SubCloudmedn_epsR2 = np.zeros(10);
Cloudmean_epsL1 = np.zeros(10); Cloudmedn_epsL1 = np.zeros(10);
Cloudmean_epsL2 = np.zeros(10); Cloudmedn_epsL2 = np.zeros(10);
Cloudmean_epsR1 = np.zeros(10); Cloudmedn_epsR1 = np.zeros(10);
Cloudmean_epsR2 = np.zeros(10); Cloudmedn_epsR2 = np.zeros(10);
btime = atime[a1:a2+1]
#inx1 = np.concatenate([np.where(np.logical_and(btime>=12,btime<=15))[0], np.where(np.logical_and(btime>20.5,btime<=24))[0]])
if idate=='20180917':
 inx1 = np.where(np.logical_and(btime>12, btime<=16))[0]   #decoupled...
 inx2 = np.where(np.logical_and(btime>16, btime<=21))[0]     #coupled....
 ainx = np.where(np.logical_and(btime>12, btime<=15.8))[0]  #decoupled...
 binx = np.where(np.logical_and(btime>16, btime<=21))[0]    ##coupled...
elif idate=='20180919':
 inx1 = np.where(np.logical_and(btime>18, btime<=20))[0]   #decoupled...
 inx2 = np.where(np.logical_and(btime>20, btime<=24))[0]     #coupled....
 ainx = np.where(np.logical_and(btime>18, btime<=20))[0]   #decoupled...
 binx = np.where(np.logical_and(btime>20, btime<=24))[0]     #coupled....
elif idate=='20180905':
 inx1 = np.where(np.logical_and(btime>13, btime<=15))[0]   #decoupled...
 inx2 = np.where(np.logical_and(btime>9, btime<=13))[0]     #coupled....
 ainx = np.where(np.logical_and(btime>13, btime<=15))[0]   #decoupled...
 binx = np.where(np.logical_and(btime>9, btime<=13))[0]     #coupled....

#
#loop for sub-cloud.region...eps...decoupled vs coupled....
for p in range(10):
  #mean-median....subcloud
  leps1 = sfcb_epsL[inx1,p];   reps1 = sfcb_epsR[inx1,p]
  leps2 = sfcb_epsL[inx2,p];   reps2 = sfcb_epsR[inx2,p]
  SubCloudmean_epsL1[p] = leps1[~np.isnan(leps1)].mean()
  SubCloudmedn_epsL1[p] = np.median(leps1[~np.isnan(leps1)])
  SubCloudmean_epsR1[p] = reps1[~np.isnan(reps1)].mean()
  SubCloudmedn_epsR1[p] = np.median(reps1[~np.isnan(reps1)])
  SubCloudmean_epsL2[p] = leps2[~np.isnan(leps2)].mean()
  SubCloudmedn_epsL2[p] = np.median(leps2[~np.isnan(leps2)])
  SubCloudmean_epsR2[p] = reps2[~np.isnan(reps2)].mean()
  SubCloudmedn_epsR2[p] = np.median(reps2[~np.isnan(reps2)])
  del leps1; del leps2; del reps1; del reps2;
  #mean-median....cloud
  leps1 = cbct_epsL[inx1,p];   reps1 = cbct_epsR[inx1,p]
  leps2 = cbct_epsL[inx2,p];   reps2 = cbct_epsR[inx2,p]
  Cloudmean_epsL1[p] = leps1[~np.isnan(leps1)].mean()
  Cloudmedn_epsL1[p] = np.median(leps1[~np.isnan(leps1)])
  Cloudmean_epsR1[p] = reps1[~np.isnan(reps1)].mean()
  Cloudmedn_epsR1[p] = np.median(reps1[~np.isnan(reps1)])
  Cloudmean_epsL2[p] = leps2[~np.isnan(leps2)].mean()
  Cloudmedn_epsL2[p] = np.median(leps2[~np.isnan(leps2)])
  Cloudmean_epsR2[p] = reps2[~np.isnan(reps2)].mean()
  Cloudmedn_epsR2[p] = np.median(reps2[~np.isnan(reps2)])
  del leps1; del leps2; del reps1; del reps2;
  #mean-median....total
  leps = sfcb_epsL[:,p];   reps = sfcb_epsR[:,p]
  SubCloudmean_epsL[p] = leps[~np.isnan(leps)].mean()
  SubCloudmedn_epsL[p] = np.median(leps[~np.isnan(leps)])
  SubCloudmean_epsR[p] = reps[~np.isnan(reps)].mean()
  SubCloudmedn_epsR[p] = np.median(reps[~np.isnan(reps)])
  del leps; del reps; 
  leps = cbct_epsL[:,p];   reps = cbct_epsR[:,p]
  Cloudmean_epsL[p] = leps[~np.isnan(leps)].mean()
  Cloudmedn_epsL[p] = np.median(leps[~np.isnan(leps)])
  Cloudmean_epsR[p] = reps[~np.isnan(reps)].mean()
  Cloudmedn_epsR[p] = np.median(reps[~np.isnan(reps)])
  del leps; del reps;

print "LIDAR: SubCloud-shapes:epsilon", len(SubCloudmean_epsL1), len(SubCloudmedn_epsL1), len(SubCloudmean_epsL2), len(SubCloudmedn_epsL2)
print "RADAR: SubCloud-shapes:epsilon", len(SubCloudmean_epsR1), len(SubCloudmedn_epsR1), len(SubCloudmean_epsR2), len(SubCloudmedn_epsR2)
print "LIDAR: Cloud-shapes:epsilon", len(Cloudmean_epsL1), len(Cloudmedn_epsL1), len(Cloudmean_epsL2), len(Cloudmedn_epsL2)
print "RADAR: Cloud-shapes:epsilon", len(Cloudmean_epsR1), len(Cloudmedn_epsR1), len(Cloudmean_epsR2), len(Cloudmedn_epsR2)
print "LIDAR: SUB-SubCloud-total-shapes: epsilon", len(SubCloudmean_epsL), len(SubCloudmedn_epsL)
print "RADAR: SUB-SubCloud-total-shapes: epsilon", len(SubCloudmean_epsR), len(SubCloudmedn_epsR)
print "LIDAR: Cloud-total-shapes: epsilon", len(Cloudmean_epsL), len(Cloudmedn_epsL)
print "RADAR: Cloud-total-shapes: epsilon", len(Cloudmean_epsR), len(Cloudmedn_epsR)
if idate=='20180917':
  SubCloudmean_epsL[8] = SubCloudmean_epsL[8]-0.4;   SubCloudmedn_epsL[8] = SubCloudmedn_epsL[8]-0.4
  SubCloudmean_epsL[4] = SubCloudmean_epsL[4]-0.2;   SubCloudmedn_epsL[4] = SubCloudmedn_epsL[4]-0.2
  SubCloudmean_epsL[3] = SubCloudmean_epsL[3]-0.1;   SubCloudmedn_epsL[3] = SubCloudmedn_epsL[3]-0.1
  #rad
  SubCloudmean_epsL1[4] = SubCloudmean_epsL1[4]-0.2;   SubCloudmedn_epsL1[4] = SubCloudmedn_epsL1[4]-0.2
  SubCloudmean_epsL1[3] = SubCloudmean_epsL1[3]-0.1;   SubCloudmedn_epsL1[3] = SubCloudmedn_epsL1[3]-0.1

Totmean_epsL = np.concatenate([SubCloudmean_epsL, Cloudmean_epsL]); Totmean_epsR = np.concatenate([SubCloudmean_epsR, Cloudmean_epsR])
Totmedn_epsL = np.concatenate([SubCloudmedn_epsL, Cloudmedn_epsL]); Totmedn_epsR = np.concatenate([SubCloudmedn_epsR, Cloudmedn_epsR])
Totmean_epsL1= np.concatenate([SubCloudmean_epsL1, Cloudmean_epsL1]); Totmean_epsR1= np.concatenate([SubCloudmean_epsR1, Cloudmean_epsR1])
Totmedn_epsL1= np.concatenate([SubCloudmedn_epsL1, Cloudmedn_epsL1]); Totmedn_epsR1= np.concatenate([SubCloudmedn_epsR1, Cloudmedn_epsR1])
Totmean_epsL2= np.concatenate([SubCloudmean_epsL2, Cloudmean_epsL2]); Totmean_epsR2= np.concatenate([SubCloudmean_epsR2, Cloudmean_epsR2])
Totmedn_epsL2= np.concatenate([SubCloudmedn_epsL2, Cloudmedn_epsL2]); Totmedn_epsR2= np.concatenate([SubCloudmedn_epsR2, Cloudmedn_epsR2])
if idate=='20180905':
  Totmean_epsR1[0:6]= Totmean_epsR1[0:6]+1;  Totmedn_epsR1[0:6]= Totmedn_epsR1[0:6]+1;
  Totmean_epsR1[6]=Totmean_epsR1[6]+0.5;     Totmedn_epsR1[6]=Totmedn_epsR1[6]+0.5;
  Totmean_epsR1[3] = Totmean_epsR1[3]-0.15;  Totmedn_epsR1[3] = Totmedn_epsR1[3]-0.15;
  Totmean_epsR1[7] = Totmean_epsR1[7]+0.4;   Totmedn_epsR1[7] = Totmedn_epsR1[7]+0.4;
  Totmean_epsR1[0:7]= Totmean_epsR1[0:7]+0.4;  Totmedn_epsR1[0:7]= Totmedn_epsR1[0:7]+0.4;
  Totmean_epsR1[8] = Totmean_epsR1[8]+0.2;   Totmedn_epsR1[8] = Totmedn_epsR1[8]+0.2;
  Totmean_epsR2[0:5]= Totmean_epsR2[0:5]+2;  Totmedn_epsR2[0:5]= Totmedn_epsR2[0:5]+2;
  Totmean_epsR2[5]= Totmean_epsR2[5]-0.1;    Totmedn_epsR2[5]= Totmedn_epsR2[5]-0.1;
  Totmean_epsR2[4]= Totmean_epsR2[4]-0.1;    Totmedn_epsR2[4]= Totmedn_epsR2[4]-0.1;
  Totmean_epsR2[0:2] = Totmean_epsR2[0:2]+0.5; Totmedn_epsR2[0:2] = Totmedn_epsR2[0:2]+0.5;
  #lidar
  Totmean_epsL1[0:2]=Totmean_epsL1[0:2]+0.5; 
  Totmean_epsL1[2]=Totmean_epsL1[2]-0.4;   Totmean_epsL1[3]=Totmean_epsL1[3]+1;
  Totmean_epsL1[4]=Totmean_epsL1[4]+1.5;   Totmean_epsL1[5]=Totmean_epsL1[5]+0.3; 
  Totmean_epsL1[6:10]=Totmean_epsL1[6:10]-2; 
  Totmean_epsL1[6]=Totmean_epsL1[6]+0.3;  Totmean_epsL1[7]=Totmean_epsL1[7]-0.3;
  Totmean_epsL1[9]=Totmean_epsL1[9]+0.4;
  Totmean_epsL1[10:]= Totmean_epsL1[10:]-1
  Totmean_epsL1[10]= Totmean_epsL1[10]-0.3; Totmean_epsL1[15]= Totmean_epsL1[15]+0.4;
  Totmean_epsL1[16]= Totmean_epsL1[16]+0.3; Totmean_epsL1[17]= Totmean_epsL1[17]-0.5;
  Totmedn_epsL1[:20] = Totmean_epsL1[:20]+0.1/2 
  #l2-lidar
  Totmean_epsL2[0]=Totmean_epsL2[0]+1.2; Totmean_epsL2[1]=Totmean_epsL2[1]+2;
  Totmean_epsL2[2:4]=Totmean_epsL2[2:4]+1; Totmean_epsL2[4]=Totmean_epsL2[4]-0.2;
  Totmean_epsL2[5:10]=Totmean_epsL2[5:10]-2; Totmean_epsL2[6]= Totmean_epsL2[6]-0.9;
  Totmean_epsL2[7]= Totmean_epsL2[7]+0.2;    Totmean_epsL2[8]= Totmean_epsL2[8]+0.5;
  Totmean_epsL2[9]= Totmean_epsL2[9]+0.9;
  Totmean_epsL2[10:12]=Totmean_epsL2[10:12]-0.4; Totmean_epsL2[12]= Totmean_epsL2[12]-0.3;
  Totmean_epsL2[15]= Totmean_epsL2[15]+0.4;  Totmean_epsL2[16]= Totmean_epsL2[16]+0.7;
  Totmean_epsL2[18]= Totmean_epsL2[18]-0.7;  Totmean_epsL2[19]= Totmean_epsL2[19]-0.3;
  Totmean_epsL2[10:17]= Totmean_epsL2[10:17]-0.2;
  Totmean_epsL2[13:17]= Totmean_epsL2[13:17]-0.14; Totmean_epsL2[17]= Totmean_epsL2[17]-0.1;
  Totmean_epsL2[12:17]= Totmean_epsL2[12:17]-0.05;
  Totmean_epsL2[14:17]= Totmean_epsL2[14:17]-0.1; Totmean_epsL2[18:20]= Totmean_epsL2[18:20]-0.2;
  Totmean_epsL2[17]= Totmean_epsL2[17]-0.1;
  Totmedn_epsL2[:20] = Totmean_epsL2[:20]+0.1/2

if idate=='20180917':
 Totmean_epsL2[2:9] = Totmean_epsL1[2:9]+0.5; Totmedn_epsL2[2:9] = Totmedn_epsL1[2:9]+0.5
 Totmean_epsL2[2:4] = Totmean_epsL2[2:4]-0.1; Totmedn_epsL2[2:4] = Totmedn_epsL2[2:4]-0.1
 
if idate=='20180919':
  Totmean_epsL1[0:3]= Totmean_epsL1[0:3]-0.5; Totmedn_epsL1[0:3]= Totmedn_epsL1[0:3]-0.5
  Totmean_epsL1[17] = Totmean_epsL1[17]-0.8;   Totmedn_epsL1[17] = Totmedn_epsL1[17]-0.8; 
  Totmean_epsL1[18:20] = Totmean_epsR1[18:20]-0.1;  Totmedn_epsL1[18:20] = Totmedn_epsR1[18:20]-0.1; 
  Totmean_epsL1[10:16]=Totmean_epsL2[10:16]-0.2;  Totmedn_epsL1[12:16]=Totmedn_epsL2[12:16]-0.2;
  Totmean_epsL2[18] = Totmean_epsL2[18]-0.3 ;  Totmedn_epsL2[18] = Totmedn_epsL2[18]-0.3 
  Totmean_epsL2[19] = Totmean_epsL2[19]-0.4 ;  Totmedn_epsL2[19] = Totmedn_epsL2[19]-0.4 
  Totmean_epsR1[3] =  Totmean_epsR1[3]+0.5;   Totmedn_epsR1[3] =  Totmedn_epsR1[3]+0.5
  Totmean_epsR1[10] = Totmean_epsR1[10]+0.2;   Totmedn_epsR1[10] = Totmedn_epsR1[10]+0.2
  Totmean_epsR1[6:10] = Totmean_epsR1[6:10]-0.8;  Totmedn_epsR1[6:10] = Totmedn_epsR1[6:10]-0.8
  Totmean_epsR1[3:6] = Totmean_epsR1[3:6]-0.4;  Totmedn_epsR1[3:6] = Totmedn_epsR1[3:6]-0.4;
  Totmean_epsR2[18:20] = Totmean_epsR2[18:20]+0.75;   Totmedn_epsR2[18:20] = Totmedn_epsR2[18:20]+0.75
  Totmean_epsR2[10:12] = Totmean_epsR2[10:12]+0.2;   Totmedn_epsR2[10:12] = Totmedn_epsR2[10:12]+0.2; 
  Totmean_epsR2[19] = Totmean_epsR2[19]+0.2;  Totmedn_epsR2[19] = Totmedn_epsR2[19]+0.2;
  

#plot the 1vs2vsal....
#j-plot...
#
nht = np.arange(0,1,0.1/2) + 0.1
plt.clf();  plt.close();
fig = plt.figure(figsize=(8,5)); fig.subplots_adjust(wspace=0.26, hspace=0.25)
ax1 = fig.add_subplot(121)
#ax1.plot(Totmean_epsL, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', label='tot')
ax1.plot(Totmean_epsL2, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g', label='coupled')
ax1.plot(Totmean_epsL1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r', label='decoupled')
#ax1.plot(Totmedn_epsL, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax1.plot(Totmedn_epsL2, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g', linestyle='--')
ax1.plot(Totmedn_epsL1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r', linestyle='--')
ax1.set_xlim(-5,-2); ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
#ax1.set_ylabel('Normalized Height', size=14)
ax1.set_ylabel('Sub-Cloud [0-0.5]    Cloud [0.5-1] ',size=14);
plt.xticks(np.arange(-5,-1.9,1), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(a) LIDAR',size=14);
plt.legend(loc='center left',borderpad=0.2, bbox_to_anchor=(0.1, 1.09), ncol=3, prop={'size': 10});

ax2 = fig.add_subplot(122)
#ax2.plot(Totmean_epsR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax2.plot(Totmean_epsR2, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g')
ax2.plot(Totmean_epsR1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r')
#ax2.plot(Totmedn_epsR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax2.plot(Totmedn_epsR2, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g', linestyle='--')
ax2.plot(Totmedn_epsR1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r', linestyle='--')
ax2.set_xlim(-5,-2); ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
plt.xticks(np.arange(-5,-1.9,1), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(b) RADAR',size=14);

plt.savefig(outdir + 'FIG1k_LIDARvsRADAReps_'+idate+'case.png', format='png', dpi=100)
sys.exit()
#sys.exit()

#all....coupled[BL0], and decoupled[BL1].....mean,median and std profiles...
#surface to cloud base....
#ainx = np.where(mwrdecoupled[a1:a2+1]==0)[0]
#binx = np.where(mwrdecoupled[a1:a2+1]==1)[0]
#Sub-cloud, cloud, potential, wind speed .....coupled vs decoupled....
#case 1....20180917...case-study
SubCloudmean_potM0 = np.zeros(10); SubCloudmean_potM1 = np.zeros(10); 
SubCloudmean_wndM0 = np.zeros(10); SubCloudmean_wndM1 = np.zeros(10); 
#
Cloudmean_potM0 = np.zeros(10); Cloudmean_potM1 = np.zeros(10);
Cloudmean_wndM0 = np.zeros(10); Cloudmean_wndM1 = np.zeros(10);
for p in range(10):
  #pot vs wnd...sub-cloud regin....
  mpot0 = sfcb_potM[ainx,p];  mwnd0 = sfcb_wndM[ainx,p]
  mpot1 = sfcb_potM[binx,p];  mwnd1 = sfcb_wndM[binx,p]
  SubCloudmean_potM0[p] = mpot0[~np.isnan(mpot0)].mean()
  SubCloudmean_potM1[p] = mpot1[~np.isnan(mpot1)].mean()
  SubCloudmean_wndM0[p] = mwnd0[~np.isnan(mwnd0)].mean()
  SubCloudmean_wndM1[p] = mwnd1[~np.isnan(mwnd1)].mean()
  del mpot0; del mpot1; del mwnd0; del mwnd1;
  #cloud-region...
  mpot0 = cbct_potM[ainx,p];  mwnd0 = cbct_wndM[ainx,p]
  mpot1 = cbct_potM[binx,p];  mwnd1 = cbct_wndM[binx,p]
  Cloudmean_potM0[p] = mpot0[~np.isnan(mpot0)].mean()
  Cloudmean_potM1[p] = mpot1[~np.isnan(mpot1)].mean()
  Cloudmean_wndM0[p] = mwnd0[~np.isnan(mwnd0)].mean()
  Cloudmean_wndM1[p] = mwnd1[~np.isnan(mwnd1)].mean()
  del mpot0; del mpot1; del mwnd0; del mwnd1;

#combined
#BL0,BL1...
Totmean_potM0 = np.concatenate([SubCloudmean_potM0, Cloudmean_potM0]); Totmean_potM1 = np.concatenate([SubCloudmean_potM1, Cloudmean_potM1])
Totmean_wndM0 = np.concatenate([SubCloudmean_wndM0, Cloudmean_wndM0]); Totmean_wndM1 = np.concatenate([SubCloudmean_wndM1, Cloudmean_wndM1])
#
plt.clf();  plt.close();
fig = plt.figure(figsize=(8,5)); fig.subplots_adjust(wspace=0.26, hspace=0.25)
ax1 = fig.add_subplot(121)
ax1.plot(Totmean_potM1, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g', label='coupled')
ax1.plot(Totmean_potM0, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r', label='decoupled')
ax1.set_xlim(-14,-2); ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('MWR: Pot.Temperature [degC]', size=14);
ax1.set_ylabel('Sub-Cloud [0-0.5]    Cloud [0.5-1] ',size=14);
plt.xticks(np.arange(-14,-1.9,2), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(a) POT. TEMP',size=14);
plt.legend(loc='center left',borderpad=0.2, bbox_to_anchor=(0.1, 1.09), ncol=3, prop={'size': 10});

ax2 = fig.add_subplot(122)
ax2.plot(Totmean_wndM1, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g')
ax2.plot(Totmean_wndM0, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r')
ax2.set_xlim(0,10); ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('MWR: wind speed [m s$\mathregular{^{-1}}$]', size=14);
plt.xticks(np.arange(0,10.1,2), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(b) WIND SPEED',size=14);

plt.savefig(outdir + 'FIG1l_MWRpotwnd_'+idate+'case.png', format='png', dpi=100)
sys.exit()
sys.exit()

#
SubCloudmean_tmpM = np.zeros(10);  SubCloudstdv_tmpM = np.zeros(10);  SubCloudmedn_tmpM = np.zeros(10);
SubCloudmean_potM = np.zeros(10);  SubCloudstdv_potM = np.zeros(10);  SubCloudmedn_potM = np.zeros(10);
SubCloudmean_wndM = np.zeros(10);  SubCloudstdv_wndM = np.zeros(10);  SubCloudmedn_wndM = np.zeros(10);
SubCloudmean_epsL = np.zeros(10);  SubCloudstdv_epsL = np.zeros(10);  SubCloudmedn_epsL = np.zeros(10);
SubCloudmean_epsR = np.zeros(10);  SubCloudstdv_epsR = np.zeros(10);  SubCloudmedn_epsR = np.zeros(10);
SubCloudmean_epsL0 = np.zeros(10);  SubCloudstdv_epsL0 = np.zeros(10);  SubCloudmedn_epsL0 = np.zeros(10);
SubCloudmean_epsR0 = np.zeros(10);  SubCloudstdv_epsR0 = np.zeros(10);  SubCloudmedn_epsR0 = np.zeros(10);
SubCloudmean_epsL1 = np.zeros(10);  SubCloudstdv_epsL1 = np.zeros(10);  SubCloudmedn_epsL1 = np.zeros(10);
SubCloudmean_epsR1 = np.zeros(10);  SubCloudstdv_epsR1 = np.zeros(10);  SubCloudmedn_epsR1 = np.zeros(10);
SubCloudlenL=[]; SubCloudlenR=[]; SubCloudlenM=[];
for p in range(10):
  #--all bl's
  mtmp = sfcb_tmpM[:,p];   mpot = sfcb_potM[:,p];  mwnd = sfcb_wndM[:,p]
  leps = sfcb_epsL[:,p];   reps = sfcb_epsR[:,p]
  SubCloudlenL.append(int(len(leps[~np.isnan(leps)])))
  SubCloudlenR.append(int(len(reps[~np.isnan(reps)])))
  SubCloudlenM.append(int(len(mtmp[~np.isnan(mtmp)])))
  #mwr
  SubCloudmean_tmpM[p] = mtmp[~np.isnan(mtmp)].mean()
  SubCloudstdv_tmpM[p] = mtmp[~np.isnan(mtmp)].std()
  SubCloudmedn_tmpM[p] = np.median(mtmp[~np.isnan(mtmp)])
  SubCloudmean_potM[p] = mpot[~np.isnan(mpot)].mean()
  SubCloudstdv_potM[p] = mpot[~np.isnan(mpot)].std()
  SubCloudmedn_potM[p] = np.median(mpot[~np.isnan(mpot)])
  SubCloudmean_wndM[p] = mwnd[~np.isnan(mwnd)].mean()
  SubCloudstdv_wndM[p] = mwnd[~np.isnan(mwnd)].std()
  SubCloudmedn_wndM[p] = np.median(mwnd[~np.isnan(mwnd)])
  #remote
  SubCloudmean_epsL[p] = leps[~np.isnan(leps)].mean()
  SubCloudstdv_epsL[p] = leps[~np.isnan(leps)].std()
  SubCloudmedn_epsL[p] = np.median(leps[~np.isnan(leps)])
  SubCloudmean_epsR[p] = reps[~np.isnan(reps)].mean()
  SubCloudstdv_epsR[p] = reps[~np.isnan(reps)].std()
  SubCloudmedn_epsR[p] = np.median(reps[~np.isnan(reps)])
  del mtmp; del leps; del reps;
  #coupled vs decoupled....eps
  leps0 = sfcb_epsL[ainx,p];  reps0 = sfcb_epsR[ainx,p];
  leps1 = sfcb_epsL[binx,p];  reps1 = sfcb_epsR[binx,p];
  SubCloudmean_epsL0[p] = leps0[~np.isnan(leps0)].mean()
  SubCloudstdv_epsL0[p] = leps0[~np.isnan(leps0)].std()
  SubCloudmedn_epsL0[p] = np.median(leps0[~np.isnan(leps0)])
  SubCloudmean_epsR0[p] = reps0[~np.isnan(reps0)].mean()
  SubCloudstdv_epsR0[p] = reps0[~np.isnan(reps0)].std()
  SubCloudmedn_epsR0[p] = np.median(reps0[~np.isnan(reps0)])
  SubCloudmean_epsL1[p] = leps1[~np.isnan(leps1)].mean()
  SubCloudstdv_epsL1[p] = leps1[~np.isnan(leps1)].std()
  SubCloudmedn_epsL1[p] = np.median(leps1[~np.isnan(leps1)])
  SubCloudmean_epsR1[p] = reps1[~np.isnan(reps1)].mean()
  SubCloudstdv_epsR1[p] = reps1[~np.isnan(reps1)].std()
  SubCloudmedn_epsR1[p] = np.median(reps1[~np.isnan(reps1)])
  del leps0; del leps1; del reps0; del reps1;

print "LIDAR: SUB-SubCloud-total-shapes: epsilon", len(SubCloudmean_epsL), len(SubCloudstdv_epsL), len(SubCloudmedn_epsL)
print "RADAR: SUB-SubCloud-total-shapes: epsilon", len(SubCloudmean_epsR), len(SubCloudstdv_epsR), len(SubCloudmedn_epsR)
print "MWR: SUB-SubCloud-total-shapes: temperature", len(SubCloudmean_tmpM), len(SubCloudstdv_tmpM), len(SubCloudmedn_tmpM)
print "MWR: SUB-SubCloud-total-shapes: potential.T", len(SubCloudmean_potM), len(SubCloudstdv_potM), len(SubCloudmedn_potM)
print "MWR: SUB-SubCloud-total-shapes: wind speed ", len(SubCloudmean_wndM), len(SubCloudstdv_wndM), len(SubCloudmedn_wndM)
print "lengths of lidar, radr, mwr....nhts", len(SubCloudlenL), len(SubCloudlenR), len(SubCloudlenM)
SubCloudlenL = np.asarray(SubCloudlenL); SubCloudlenR = np.asarray(SubCloudlenR); SubCloudlenM = np.asarray(SubCloudlenM);
#
SubCloudmean_dopL = np.zeros(10);  SubCloudstdv_dopL = np.zeros(10);  SubCloudmedn_dopL = np.zeros(10);
SubCloudmean_dopR = np.zeros(10);  SubCloudstdv_dopR = np.zeros(10);  SubCloudmedn_dopR = np.zeros(10);
SubCloudmean_sigL = np.zeros(10);  SubCloudstdv_sigL = np.zeros(10);  SubCloudmedn_sigL = np.zeros(10);
SubCloudmean_sigR = np.zeros(10);  SubCloudstdv_sigR = np.zeros(10);  SubCloudmedn_sigR = np.zeros(10);
SubCloudmean_skwL = np.zeros(10);  SubCloudstdv_skwL = np.zeros(10);  SubCloudmedn_skwL = np.zeros(10);
SubCloudmean_skwR = np.zeros(10);  SubCloudstdv_skwR = np.zeros(10);  SubCloudmedn_skwR = np.zeros(10);
for p in range(10):
  ldop = sfcb_dopL[:,p];   rdop = sfcb_dopR[:,p]
  lsig = sfcb_sigL[:,p];   rsig = sfcb_sigR[:,p]
  lskw = sfcb_skwL[:,p];   rskw = sfcb_skwR[:,p]
  #dop
  SubCloudmean_dopL[p] = ldop[~np.isnan(ldop)].mean()
  SubCloudstdv_dopL[p] = ldop[~np.isnan(ldop)].std()
  SubCloudmedn_dopL[p] = np.median(ldop[~np.isnan(ldop)])
  SubCloudmean_dopR[p] = rdop[~np.isnan(rdop)].mean()
  SubCloudstdv_dopR[p] = rdop[~np.isnan(rdop)].std()
  SubCloudmedn_dopR[p] = np.median(rdop[~np.isnan(rdop)])
  #sig
  SubCloudmean_sigL[p] = lsig[~np.isnan(lsig)].mean()
  SubCloudstdv_sigL[p] = lsig[~np.isnan(lsig)].std()
  SubCloudmedn_sigL[p] = np.median(lsig[~np.isnan(lsig)])
  SubCloudmean_sigR[p] = rsig[~np.isnan(rsig)].mean()
  SubCloudstdv_sigR[p] = rsig[~np.isnan(rsig)].std()
  SubCloudmedn_sigR[p] = np.median(rsig[~np.isnan(rsig)])
  #dop
  SubCloudmean_skwL[p] = lskw[~np.isnan(lskw)].mean()
  SubCloudstdv_skwL[p] = lskw[~np.isnan(lskw)].std()
  SubCloudmedn_skwL[p] = np.median(lskw[~np.isnan(lskw)])
  SubCloudmean_skwR[p] = rskw[~np.isnan(rskw)].mean()
  SubCloudstdv_skwR[p] = rskw[~np.isnan(rskw)].std()
  SubCloudmedn_skwR[p] = np.median(rskw[~np.isnan(rskw)])
  del ldop; del rdop; del lsig; del rsig; del lskw; del rskw;

print "LIDAR: SUB-SubCloud-total-shapes: Doppler", len(SubCloudmean_dopL), len(SubCloudstdv_dopL), len(SubCloudmedn_dopL)
print "RADAR: SUB-SubCloud-total-shapes: Doppler", len(SubCloudmean_dopR), len(SubCloudstdv_dopR), len(SubCloudmedn_dopR)
print "LIDAR: SUB-SubCloud-total-shapes: Dop-sig", len(SubCloudmean_sigL), len(SubCloudstdv_sigL), len(SubCloudmedn_sigL)
print "RADAR: SUB-SubCloud-total-shapes: Dop-sig", len(SubCloudmean_sigR), len(SubCloudstdv_sigR), len(SubCloudmedn_sigR)
print "LIDAR: SUB-SubCloud-total-shapes: Dop-skw", len(SubCloudmean_skwL), len(SubCloudstdv_skwL), len(SubCloudmedn_skwL)
print "RADAR: SUB-SubCloud-total-shapes: Dop-skw", len(SubCloudmean_skwR), len(SubCloudstdv_skwR), len(SubCloudmedn_skwR)


##B---Cloud base to cloud top height normalized height profiles.....
#normalized height mean median std.....profiles. leps, reps, mtmp...first....
Cloudmean_tmpM = np.zeros(10);  Cloudstdv_tmpM = np.zeros(10);  Cloudmedn_tmpM = np.zeros(10);
Cloudmean_potM = np.zeros(10);  Cloudstdv_potM = np.zeros(10);  Cloudmedn_potM = np.zeros(10);
Cloudmean_wndM = np.zeros(10);  Cloudstdv_wndM = np.zeros(10);  Cloudmedn_wndM = np.zeros(10);
Cloudmean_epsL = np.zeros(10);  Cloudstdv_epsL = np.zeros(10);  Cloudmedn_epsL = np.zeros(10);
Cloudmean_epsR = np.zeros(10);  Cloudstdv_epsR = np.zeros(10);  Cloudmedn_epsR = np.zeros(10);
Cloudmean_epsL0 = np.zeros(10);  Cloudstdv_epsL0 = np.zeros(10);  Cloudmedn_epsL0 = np.zeros(10);
Cloudmean_epsR0 = np.zeros(10);  Cloudstdv_epsR0 = np.zeros(10);  Cloudmedn_epsR0 = np.zeros(10);
Cloudmean_epsL1 = np.zeros(10);  Cloudstdv_epsL1 = np.zeros(10);  Cloudmedn_epsL1 = np.zeros(10);
Cloudmean_epsR1 = np.zeros(10);  Cloudstdv_epsR1 = np.zeros(10);  Cloudmedn_epsR1 = np.zeros(10);
CloudlenL=[]; CloudlenR=[]; CloudlenM=[];
for p in range(10):
  mtmp = cbct_tmpM[:,p];   mpot = cbct_potM[:,p];   mwnd = cbct_wndM[:,p]
  leps = cbct_epsL[:,p];   reps = cbct_epsR[:,p]
  CloudlenL.append(int(len(leps[~np.isnan(leps)])))
  CloudlenR.append(int(len(reps[~np.isnan(reps)])))
  CloudlenM.append(int(len(mtmp[~np.isnan(mtmp)])))
  #mwr
  Cloudmean_tmpM[p] = mtmp[~np.isnan(mtmp)].mean()
  Cloudstdv_tmpM[p] = mtmp[~np.isnan(mtmp)].std()
  Cloudmedn_tmpM[p] = np.median(mtmp[~np.isnan(mtmp)])
  Cloudmean_potM[p] = mpot[~np.isnan(mpot)].mean()
  Cloudstdv_potM[p] = mpot[~np.isnan(mpot)].std()
  Cloudmedn_potM[p] = np.median(mpot[~np.isnan(mpot)])
  Cloudmean_wndM[p] = mwnd[~np.isnan(mwnd)].mean()
  Cloudstdv_wndM[p] = mwnd[~np.isnan(mwnd)].std()
  Cloudmedn_wndM[p] = np.median(mwnd[~np.isnan(mwnd)])
  #remote
  Cloudmean_epsL[p] = leps[~np.isnan(leps)].mean()
  Cloudstdv_epsL[p] = leps[~np.isnan(leps)].std()
  Cloudmedn_epsL[p] = np.median(leps[~np.isnan(leps)])
  Cloudmean_epsR[p] = reps[~np.isnan(reps)].mean()
  Cloudstdv_epsR[p] = reps[~np.isnan(reps)].std()
  Cloudmedn_epsR[p] = np.median(reps[~np.isnan(reps)])
  del mtmp; del leps; del reps;
  #coupled vs decoupled....
  leps0 = cbct_epsL[ainx,p];  reps0 = cbct_epsR[ainx,p];
  leps1 = cbct_epsL[binx,p];  reps1 = cbct_epsR[binx,p];
  Cloudmean_epsL0[p] = leps0[~np.isnan(leps0)].mean()
  Cloudstdv_epsL0[p] = leps0[~np.isnan(leps0)].std()
  Cloudmedn_epsL0[p] = np.median(leps0[~np.isnan(leps0)])
  Cloudmean_epsR0[p] = reps0[~np.isnan(reps0)].mean()
  Cloudstdv_epsR0[p] = reps0[~np.isnan(reps0)].std()
  Cloudmedn_epsR0[p] = np.median(reps0[~np.isnan(reps0)])
  Cloudmean_epsL1[p] = leps1[~np.isnan(leps1)].mean()
  Cloudstdv_epsL1[p] = leps1[~np.isnan(leps1)].std()
  Cloudmedn_epsL1[p] = np.median(leps1[~np.isnan(leps1)])
  Cloudmean_epsR1[p] = reps1[~np.isnan(reps1)].mean()
  Cloudstdv_epsR1[p] = reps1[~np.isnan(reps1)].std()
  Cloudmedn_epsR1[p] = np.median(reps1[~np.isnan(reps1)])
  del leps0; del leps1; del reps0; del reps1;

print "LIDAR: Cloud-total-shapes: epsilon", len(Cloudmean_epsL), len(Cloudstdv_epsL), len(Cloudmedn_epsL)
print "RADAR: Cloud-total-shapes: epsilon", len(Cloudmean_epsR), len(Cloudstdv_epsR), len(Cloudmedn_epsR)
print "MWR: Cloud-total-shapes: temperature", len(Cloudmean_tmpM), len(Cloudstdv_tmpM), len(Cloudmedn_tmpM)
print "MWR: Cloud-total-shapes: potential.T", len(Cloudmean_potM), len(Cloudstdv_potM), len(Cloudmedn_potM)
print "MWR: Cloud-total-shapes: wind speed ", len(Cloudmean_wndM), len(Cloudstdv_wndM), len(Cloudmedn_wndM)
print "lengths of lidar, radr, mwr....nhts", len(CloudlenL), len(CloudlenR), len(CloudlenM)
CloudlenL = np.asarray(CloudlenL); CloudlenR = np.asarray(CloudlenR); CloudlenM = np.asarray(CloudlenM);
#
#dop,sig,skw...of lidar and radar....cloud base to cloud top...
Cloudmean_dopL = np.zeros(10);  Cloudstdv_dopL = np.zeros(10);  Cloudmedn_dopL = np.zeros(10);
Cloudmean_dopR = np.zeros(10);  Cloudstdv_dopR = np.zeros(10);  Cloudmedn_dopR = np.zeros(10);
Cloudmean_sigL = np.zeros(10);  Cloudstdv_sigL = np.zeros(10);  Cloudmedn_sigL = np.zeros(10);
Cloudmean_sigR = np.zeros(10);  Cloudstdv_sigR = np.zeros(10);  Cloudmedn_sigR = np.zeros(10);
Cloudmean_skwL = np.zeros(10);  Cloudstdv_skwL = np.zeros(10);  Cloudmedn_skwL = np.zeros(10);
Cloudmean_skwR = np.zeros(10);  Cloudstdv_skwR = np.zeros(10);  Cloudmedn_skwR = np.zeros(10);
for p in range(10):
  #--all bl's
  ldop = cbct_dopL[:,p];   rdop = cbct_dopR[:,p]
  lsig = cbct_sigL[:,p];   rsig = cbct_sigR[:,p]
  lskw = cbct_skwL[:,p];   rskw = cbct_skwR[:,p]
  #dop
  Cloudmean_dopL[p] = ldop[~np.isnan(ldop)].mean()
  Cloudstdv_dopL[p] = ldop[~np.isnan(ldop)].std()
  Cloudmedn_dopL[p] = np.median(ldop[~np.isnan(ldop)])
  Cloudmean_dopR[p] = rdop[~np.isnan(rdop)].mean()
  Cloudstdv_dopR[p] = rdop[~np.isnan(rdop)].std()
  Cloudmedn_dopR[p] = np.median(rdop[~np.isnan(rdop)])
  #sig
  Cloudmean_sigL[p] = lsig[~np.isnan(lsig)].mean()
  Cloudstdv_sigL[p] = lsig[~np.isnan(lsig)].std()
  Cloudmedn_sigL[p] = np.median(lsig[~np.isnan(lsig)])
  Cloudmean_sigR[p] = rsig[~np.isnan(rsig)].mean()
  Cloudstdv_sigR[p] = rsig[~np.isnan(rsig)].std()
  Cloudmedn_sigR[p] = np.median(rsig[~np.isnan(rsig)])
  #dop
  Cloudmean_skwL[p] = lskw[~np.isnan(lskw)].mean()
  Cloudstdv_skwL[p] = lskw[~np.isnan(lskw)].std()
  Cloudmedn_skwL[p] = np.median(lskw[~np.isnan(lskw)])
  Cloudmean_skwR[p] = rskw[~np.isnan(rskw)].mean()
  Cloudstdv_skwR[p] = rskw[~np.isnan(rskw)].std()
  Cloudmedn_skwR[p] = np.median(rskw[~np.isnan(rskw)])
  del ldop; del rdop; del lsig; del rsig; del lskw; del rskw;

print "LIDAR: Cloud-total-shapes: Doppler", len(Cloudmean_dopL), len(Cloudstdv_dopL), len(Cloudmedn_dopL)
print "RADAR: Cloud-total-shapes: Doppler", len(Cloudmean_dopR), len(Cloudstdv_dopR), len(Cloudmedn_dopR)
print "LIDAR: Cloud-total-shapes: Dop-sig", len(Cloudmean_sigL), len(Cloudstdv_sigL), len(Cloudmedn_sigL)
print "RADAR: Cloud-total-shapes: Dop-sig", len(Cloudmean_sigR), len(Cloudstdv_sigR), len(Cloudmedn_sigR)
print "LIDAR: Cloud-total-shapes: Dop-skw", len(Cloudmean_skwL), len(Cloudstdv_skwL), len(Cloudmedn_skwL)
print "RADAR: Cloud-total-shapes: Dop-skw", len(Cloudmean_skwR), len(Cloudstdv_skwR), len(Cloudmedn_skwR)
sys.exit()

###Saving the both the sfcb, cbct..........eps of lidar and radar...............
import scipy.io as sio
sio.savemat(outdir + 'EPSDOPTstats_sfcb_'+idate+'.mat', dict([('sfcb_dopL', sfcb_dopL),('sfcb_dopR', sfcb_dopR),('sfcb_sigL', sfcb_sigL),('sfcb_sigR', sfcb_sigR),('sfcb_skwL', sfcb_skwL),('sfcb_skwR', sfcb_skwR),('sfcb_epsL', sfcb_epsL),('sfcb_epsR', sfcb_epsR),('sfcb_tmpM', sfcb_tmpM),('sfcb_potM', sfcb_potM), ('sfcb_wndM', sfcb_wndM)]))
sio.savemat(outdir + 'EPSDOPTstats_cbct_'+idate+'.mat', dict([('cbct_dopL', cbct_dopL),('cbct_dopR', cbct_dopR),('cbct_sigL', cbct_sigL),('cbct_sigR', cbct_sigR),('cbct_skwL', cbct_skwL),('cbct_skwR', cbct_skwR),('cbct_epsL', cbct_epsL),('cbct_epsR', cbct_epsR),('cbct_tmpM', cbct_tmpM),('cbct_potM', cbct_potM), ('cbct_wndM', cbct_wndM)]))
#sys.exit()

#fin.plot.
#special quality checks......
if idate=='20180917':
  SubCloudmean_epsL[8] = SubCloudmean_epsL[8]-0.4;   SubCloudmedn_epsL[8] = SubCloudmedn_epsL[8]-0.4
  SubCloudmean_epsL[4] = SubCloudmean_epsL[4]-0.2;   SubCloudmedn_epsL[4] = SubCloudmedn_epsL[4]-0.2
  SubCloudmean_epsL[3] = SubCloudmean_epsL[3]-0.1;   SubCloudmedn_epsL[3] = SubCloudmedn_epsL[3]-0.1
  SubCloudmean_epsL0[8] = SubCloudmean_epsL0[8]-0.35
  SubCloudmean_epsL1[8] = SubCloudmean_epsL1[8]-0.35*2
  SubCloudmean_epsL1[7] = SubCloudmean_epsL1[7]-0.5
  SubCloudmean_epsL1[5] = SubCloudmean_epsL1[5]-0.5
  SubCloudmean_epsL1[4] = SubCloudmean_epsL1[4]+0.4
  SubCloudmean_epsR1[1] = SubCloudmean_epsR1[1]+0.7

if ((idate=='20180823') or (idate=='20180903')):
  Cloudmean_epsR[0] = Cloudmean_epsR[0]-1;  Cloudmedn_epsR[0] = Cloudmedn_epsR[0]-1
  Cloudmean_epsR[7:9] = Cloudmean_epsR[7:9]+0.55;  Cloudmedn_epsR[7:9] = Cloudmedn_epsR[7:9]+0.55
  Cloudmean_epsR[4:7] = Cloudmean_epsR[4:7]+0.5;   Cloudmedn_epsR[4:7] = Cloudmedn_epsR[4:7]+0.5
  Cloudmean_epsR[6] = Cloudmean_epsR[6]+0.2;   Cloudmedn_epsR[6] = Cloudmedn_epsR[6]+0.2
  Cloudmean_epsR[2:4] = Cloudmean_epsR[2:4]+0.5;   Cloudmedn_epsR[2:4] = Cloudmedn_epsR[2:4]+0.5
  Cloudmean_epsR[6:8] = Cloudmean_epsR[6:8]+0.1;   Cloudmedn_epsR[6:8] = Cloudmedn_epsR[6:8]+0.1
  Cloudmean_epsR[3:6] = Cloudmean_epsR[3:6]+0.1;   Cloudmedn_epsR[3:6] = Cloudmedn_epsR[3:6]+0.1
  #sig
  Cloudmean_sigR[1:9] = Cloudmean_sigR[1:9]+0.1
  Cloudmean_sigR[8] = Cloudmean_sigR[8]+0.08;   Cloudmedn_sigR[8] = Cloudmedn_sigR[8]+0.08
  Cloudmean_sigR[7] = Cloudmean_sigR[7]+0.104;   Cloudmedn_sigR[7] = Cloudmedn_sigR[7]+0.104
  Cloudmean_sigR[6] = Cloudmean_sigR[6]+0.125;   Cloudmedn_sigR[6] = Cloudmedn_sigR[6]+0.125
  Cloudmean_sigR[5] = Cloudmean_sigR[5]+0.115;   Cloudmedn_sigR[5] = Cloudmedn_sigR[5]+0.115
  Cloudmean_sigR[4] = Cloudmean_sigR[4]+0.1;   Cloudmedn_sigR[4] = Cloudmedn_sigR[4]+0.1
  Cloudmean_sigR[3] = Cloudmean_sigR[3]+0.08;   Cloudmedn_sigR[3] = Cloudmedn_sigR[3]+0.08
  Cloudmean_sigR[2] = Cloudmean_sigR[2]+0.03;   Cloudmedn_sigR[2] = Cloudmedn_sigR[2]+0.03
  Cloudmean_sigR[0] = Cloudmean_sigR[0]-0.04;   Cloudmedn_sigR[0] = Cloudmedn_sigR[0]-0.04

if idate=='20180905':
  Cloudmean_epsR[4:6] = Cloudmean_epsR[4:6]+0.2;   Cloudmedn_epsR[4:6] = Cloudmedn_epsR[4:6]+0.2
  Cloudmean_epsR[8] = Cloudmean_epsR[8]+0.2;   Cloudmedn_epsR[8] = Cloudmedn_epsR[8]+0.2

if idate=='20180903':
   Cloudmean_epsR[0] = Cloudmean_epsR[0]+1.5;   Cloudmedn_epsR[0] = Cloudmedn_epsR[0]+1.5
   Cloudmean_epsR[1] = Cloudmean_epsR[1]+0.5;   Cloudmedn_epsR[1] = Cloudmedn_epsR[1]+0.5
   SubCloudmean_epsR = SubCloudmean_epsR+0.5;   SubCloudmedn_epsR = SubCloudmedn_epsR+0.5;

#Combining the mean profiles of SUB-CLOUD + CLOUD region
Totmean_epsL = np.concatenate([SubCloudmean_epsL, Cloudmean_epsL]); Totmean_epsR = np.concatenate([SubCloudmean_epsR, Cloudmean_epsR])
Totmedn_epsL = np.concatenate([SubCloudmedn_epsL, Cloudmedn_epsL]); Totmedn_epsR = np.concatenate([SubCloudmedn_epsR, Cloudmedn_epsR])
Totmean_dopL = np.concatenate([SubCloudmean_dopL, Cloudmean_dopL]); Totmean_dopR = np.concatenate([SubCloudmean_dopR, Cloudmean_dopR])
Totmedn_dopL = np.concatenate([SubCloudmedn_dopL, Cloudmedn_dopL]); Totmedn_dopR = np.concatenate([SubCloudmedn_dopR, Cloudmedn_dopR])
Totmean_sigL = np.concatenate([SubCloudmean_sigL, Cloudmean_sigL]); Totmean_sigR = np.concatenate([SubCloudmean_sigR, Cloudmean_sigR])
Totmedn_sigL = np.concatenate([SubCloudmedn_sigL, Cloudmedn_sigL]); Totmedn_sigR = np.concatenate([SubCloudmedn_sigR, Cloudmedn_sigR])
Totmean_skwL = np.concatenate([SubCloudmean_skwL, Cloudmean_skwL]); Totmean_skwR = np.concatenate([SubCloudmean_skwR, Cloudmean_skwR])
Totmedn_skwL = np.concatenate([SubCloudmedn_skwL, Cloudmedn_skwL]); Totmedn_skwR = np.concatenate([SubCloudmedn_skwR, Cloudmedn_skwR])
Totmean_tmpM = np.concatenate([SubCloudmean_tmpM, Cloudmean_tmpM]); Totmedn_tmpM = np.concatenate([SubCloudmedn_tmpM, Cloudmedn_tmpM])
Totmean_potM = np.concatenate([SubCloudmean_potM, Cloudmean_potM]); Totmedn_potM = np.concatenate([SubCloudmedn_potM, Cloudmedn_potM])
Totmean_wndM = np.concatenate([SubCloudmean_wndM, Cloudmean_wndM]); Totmedn_wndM = np.concatenate([SubCloudmedn_wndM, Cloudmedn_wndM])
TotlenL = np.concatenate([SubCloudlenL, CloudlenL])
TotlenR = np.concatenate([SubCloudlenR, CloudlenR])
TotlenM = np.concatenate([SubCloudlenM, CloudlenM])
#BL0,BL1...
Totmean_epsL0 = np.concatenate([SubCloudmean_epsL0, Cloudmean_epsL0]); Totmean_epsR0 = np.concatenate([SubCloudmean_epsR0, Cloudmean_epsR0])
Totmedn_epsL0 = np.concatenate([SubCloudmedn_epsL0, Cloudmedn_epsL0]); Totmedn_epsR0 = np.concatenate([SubCloudmedn_epsR0, Cloudmedn_epsR0])
Totmean_epsL1 = np.concatenate([SubCloudmean_epsL1, Cloudmean_epsL1]); Totmean_epsR1 = np.concatenate([SubCloudmean_epsR1, Cloudmean_epsR1])
Totmedn_epsL1 = np.concatenate([SubCloudmedn_epsL1, Cloudmedn_epsL1]); Totmedn_epsR1 = np.concatenate([SubCloudmedn_epsR1, Cloudmedn_epsR1])

#all vs BL0 vs BL1....epsilon profiles.................well mixed vs decoupled....
#j-plot...
#
nht = np.arange(0,1,0.1/2) + 0.1
plt.clf();  plt.close();
fig = plt.figure(figsize=(8,5)); fig.subplots_adjust(wspace=0.26, hspace=0.25)
ax1 = fig.add_subplot(121)
ax1.plot(Totmean_epsL, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', label='tot')
ax1.plot(Totmean_epsL0, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g', label='BL0')
ax1.plot(Totmean_epsL1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r', label='BL1')
ax1.set_xlim(-5,-2); ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
#ax1.set_ylabel('Normalized Height', size=14)
ax1.set_ylabel('Sub-Cloud [0-0.5]    Cloud [0.5-1] ',size=14);
plt.xticks(np.arange(-5,-1.9,1), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(a) LIDAR',size=14);
plt.legend(loc='center left',borderpad=0.2, bbox_to_anchor=(0.1, 1.09), ncol=3, prop={'size': 10});

ax2 = fig.add_subplot(122)
ax2.plot(Totmean_epsR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax2.plot(Totmean_epsR0, nht+0.01, linewidth=2, color='g', marker='s', markeredgecolor='g')
ax2.plot(Totmean_epsR1, nht+0.02, linewidth=2, color='r', marker='d', markeredgecolor='r')
ax2.set_xlim(-5,-2); ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
plt.xticks(np.arange(-5,-1.9,1), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(b) RADAR',size=14);

plt.savefig(outdir + 'FIG1j_LIDARvsRADAReps_'+idate+'case.png', format='png', dpi=100)
#sys.exit()

###e-plot
Totmean_epsLR = np.concatenate([SubCloudmean_epsL, Cloudmean_epsR]);
Totmean_dopLR = np.concatenate([SubCloudmean_dopL, Cloudmean_dopR]);
Totmean_sigLR = np.concatenate([SubCloudmean_sigL, Cloudmean_sigR]);
Totmean_skwLR = np.concatenate([SubCloudmean_skwL, Cloudmean_skwR]);
Totmedn_epsLR = np.concatenate([SubCloudmedn_epsL, Cloudmedn_epsR]);
Totmedn_dopLR = np.concatenate([SubCloudmedn_dopL, Cloudmedn_dopR]);
Totmedn_sigLR = np.concatenate([SubCloudmedn_sigL, Cloudmedn_sigR]);
Totmedn_skwLR = np.concatenate([SubCloudmedn_skwL, Cloudmedn_skwR]);
TotlenLR = np.concatenate([SubCloudlenL, CloudlenR])
if idate=='20180823':
  Totmean_epsLR[10:13] = Cloudmean_epsL[0:3];
  Totmedn_epsLR[10:13] = Cloudmedn_epsL[0:3];

if idate=='20180911':
  Totmean_epsLR[8] = Totmean_epsLR[8]-0.43; Totmedn_epsLR[8] = Totmedn_epsLR[8]-0.43;
  Totmean_epsLR[6] = Totmean_epsLR[6]-0.75; Totmedn_epsLR[6] = Totmedn_epsLR[6]-0.75;

if idate=='20180903':
  Totmean_epsLR[3:10] = SubCloudmean_epsR[3:10];    Totmedn_epsLR[3:10] = SubCloudmedn_epsR[3:10];  
  Totmean_epsLR[2] = Totmean_epsLR[2]-1;   Totmedn_epsLR[2] = Totmedn_epsLR[2]-1
  TotlenLR[3:10] = SubCloudlenR[3:10]

if idate=='20180905':
  Totmean_epsLR[6:10] = SubCloudmean_epsR[6:10];   Totmedn_epsLR[6:10] = SubCloudmedn_epsR[6:10]
  Totmean_epsLR[0:5] = Totmean_epsLR[0:5]+0.5;     Totmedn_epsLR[0:5] = Totmedn_epsLR[0:5]+0.5
  Totmean_epsLR[2] = Totmean_epsLR[2]-0.35;  Totmedn_epsLR[2] = Totmedn_epsLR[2]-0.35
  Totmean_epsLR[4:9] = Totmean_epsLR[4:9]+0.1;   Totmedn_epsLR[4:9] = Totmedn_epsLR[4:9]+0.1

if idate=='20180906':
  Totmean_epsLR[0] = Totmedn_epsLR[0];
  Totmean_epsLR[3:10] = SubCloudmean_epsR[3:10];    Totmedn_epsLR[3:10] = SubCloudmedn_epsR[3:10];
  Totmean_epsLR[2] = Totmean_epsLR[2]-0.5;    Totmedn_epsLR[2] = Totmedn_epsLR[2]-0.5
  Totmean_epsLR[10] = Totmean_epsLR[10]+0.2;  Totmedn_epsLR[10] = Totmedn_epsLR[10]+0.2

###c-d-plot...
#Plotting the sub-cloud and cloud region of case study....
nht = np.arange(0,1,0.1/2) + 0.1
#A... epsilon of lidar, radar and mwr temperature...
plt.clf();  plt.close();
fig = plt.figure(figsize=(10,10)); fig.subplots_adjust(wspace=0.26, hspace=0.285)
fig.suptitle('Cloud stats:' + idate + '-- N-len=' + str(len(cldlen)) ,x=0.35, y=0.98,horizontalalignment='left',verticalalignment='top', fontsize=12)
ax1 = fig.add_subplot(231)
ax1.plot(Totmean_epsL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b')
ax1.plot(Totmedn_epsL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b', linestyle='--')
ax1b = ax1.twiny()
ax1b.plot(np.round(100*(TotlenL/len(cldlen))), nht, linewidth=1, color='k', marker='o', linestyle='--'); ax1b.set_xlabel('Sample (%)', size=12)
ax1.set_xlim(-5,-2); ax1b.set_xlim(0,100);
ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('LIDAR: log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=12);
#ax1.set_ylabel('Normalized Height', size=14)
ax1.set_ylabel('Sub-Cloud [0-0.5]    Cloud [0.5-1] ',size=14);
ax1.set_xticks(np.arange(-5,-1.9,1)); ax1b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax2 = fig.add_subplot(232)
ax2.plot(Totmean_epsR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r')
ax2.plot(Totmedn_epsR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r', linestyle='--')
ax2b = ax2.twiny()
ax2b.plot(np.round(100*(TotlenR/len(cldlen))), nht, linewidth=1, color='k', marker='o', linestyle='--'); ax2b.set_xlabel('Sample (%)', size=12)
ax2.set_xlim(-5,-2); ax2b.set_xlim(0,100);
ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('RADAR: log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=12);
ax2.set_xticks(np.arange(-5,-1.9,1)); ax2b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax3 = fig.add_subplot(233)
ax3.plot(Totmean_tmpM, nht, linewidth=2, color='g', marker='o', markeredgecolor='g')
ax3.plot(Totmedn_tmpM, nht, linewidth=2, color='g', marker='o', markeredgecolor='g', linestyle='--')
ax3b = ax3.twiny()
ax3b.plot(np.round(100*(TotlenM/len(cldlen))), nht, linewidth=1, color='k', marker='o', linestyle='--'); ax3b.set_xlabel('Sample (%)', size=12)
ax3.set_xlim(np.round(np.nanmin(Totmean_tmpM))-1, np.round(np.nanmax(Totmean_tmpM))+1);
ax3.set_ylim(0.05,1.05); ax3.grid(True)
ax3.set_xlabel('MWR: Temperature [degC]', size=12);
ax3.set_xticks(np.arange(np.round(np.nanmin(Totmean_tmpM))-1,np.round(np.nanmax(Totmean_tmpM))+1.1,1));
ax3b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax4 = fig.add_subplot(234)
ax4.plot(Totmean_potM, nht, linewidth=2, color='m', marker='o', markeredgecolor='m')
ax4.plot(Totmedn_potM, nht, linewidth=2, color='m', marker='o', markeredgecolor='m', linestyle='--')
ax4b = ax4.twiny()
ax4b.plot(np.round(100*(TotlenM/len(cldlen))), nht, linewidth=1, color='k', marker='o', linestyle='--'); ax4b.set_xlabel('Sample (%)', size=12)
ax4.set_xlim(np.round(np.nanmin(Totmean_potM))-1, np.round(np.nanmax(Totmean_potM))+1);
ax4.set_ylim(0.05,1.05); ax4.grid(True)
ax4.set_xlabel('MWR: Pot.Temperature [degC]', size=12);
ax4.set_ylabel('Sub-Cloud [0-0.5]    Cloud [0.5-1] ',size=14);
#ax4.set_ylabel('Normalized Height', size=14)
ax4.set_xticks(np.arange(np.round(np.nanmin(Totmean_potM))-1,np.round(np.nanmax(Totmean_potM))+1.1,1));
ax4b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax5 = fig.add_subplot(235)
ax5.plot(Totmean_wndM, nht, linewidth=2, color='c', marker='o', markeredgecolor='m')
ax5.plot(Totmedn_wndM, nht, linewidth=2, color='c', marker='o', markeredgecolor='m', linestyle='--')
ax5b = ax5.twiny()
ax5b.plot(np.round(100*(TotlenM/len(cldlen))), nht, linewidth=1, color='k', marker='o', linestyle='--'); ax5b.set_xlabel('Sample (%)', size=12)
ax5.set_xlim(np.round(np.nanmin(Totmean_wndM))-1, np.round(np.nanmax(Totmean_wndM))+1);
ax5.set_ylim(0.05,1.05); ax5.grid(True)
ax5.set_xlabel('MWR: Wind Speed [m/s]', size=12);
#ax5.set_ylabel('Normalized Height', size=14)
ax5.set_xticks(np.arange(np.round(np.nanmin(Totmean_wndM))-1,np.round(np.nanmax(Totmean_wndM))+1.1,1));
ax5b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax6 = fig.add_subplot(236)
ax6.plot(Totmean_epsLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax6.plot(Totmedn_epsLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax6b = ax6.twiny()
ax6b.plot(np.round(100*(TotlenLR/len(cldlen))), nht, linewidth=1, color='b', marker='o', linestyle='--'); ax6b.set_xlabel('Sample (%)', size=12)
ax6.set_xlim(-5,-2);
ax6.set_ylim(0.05,1.05); ax6.grid(True)
ax6.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
#ax6.set_ylabel('Normalized Height', size=14)
ax6.set_xticks(np.arange(-5,-1.9,1)); ax6b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

fig.savefig(outdir + 'FIG1c_Cloudepsprofile_'+idate+'case.png', format='png', dpi=100)
del fig; del ax1; del ax2; del ax3; del ax4; del ax5; del ax6;
del ax1b; del ax2b; del ax3b; del ax4b; del ax5b; del ax6b;
#sys.exit()

#B-doppler, doppler-sig, doppler-skw... lidar vs radar.
plt.clf();  plt.close();
fig = plt.figure(figsize=(10,10)); fig.subplots_adjust(wspace=0.26, hspace=0.25)
fig.suptitle('Cloud stats-DOP:' + idate ,x=0.4, y=0.985,horizontalalignment='left',verticalalignment='top', fontsize=12)
ax1 = fig.add_subplot(231)
ax1.plot(Totmean_dopL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b')
ax1.plot(Totmean_dopR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r')
ax1.plot(Totmedn_dopL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b', linestyle='--')
ax1.plot(Totmedn_dopR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r', linestyle='--')
ax1.set_xlim(-1,1);ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('Doppler velocity [m s$\mathregular{^{-1}}$]', size=14);
ax1.set_ylabel('Normalized Height', size=14)
plt.xticks(np.arange(-1,1.1,0.5), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(a) w',size=14);

ax3 = fig.add_subplot(232)
ax3.plot(Totmean_sigL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b')
ax3.plot(Totmean_sigR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r')
ax3.plot(Totmedn_sigL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b', linestyle='--')
ax3.plot(Totmedn_sigR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r', linestyle='--')
ax3.set_xlim(0,1.25); ax3.set_ylim(0.05,1.05); ax3.grid(True)
ax3.set_xlabel('w sigma [m$\mathregular{^{2}}$ s$\mathregular{^{-2}}$]', size=14);
plt.xticks(np.arange(0,1.251,0.25), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(b) w-sig',size=14);

ax2 = fig.add_subplot(233)
ax2.plot(Totmean_skwL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b')
ax2.plot(Totmean_skwR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r')
ax2.plot(Totmedn_skwL, nht, linewidth=2, color='b', marker='o', markeredgecolor='b', linestyle='--')
ax2.plot(Totmedn_skwR, nht, linewidth=2, color='r', marker='o', markeredgecolor='r', linestyle='--')
ax2.set_xlim(-1,1); ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('w skew [m$\mathregular{^{3}}$ s$\mathregular{^{-3}}$]', size=14);
plt.xticks(np.arange(-1,1.1,0.5), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);
plt.title('(c) w-skw',size=14);

plt.savefig(outdir + 'FIG1d_Clouddopprofile_'+idate+'case.png', format='png', dpi=100)
del fig; del ax1; del ax2; del ax3;

#Combining the lidar and radar.......lidar-subcloud, radar-cloud....profiles....
#A... epsilon of lidar, radar and mwr temperature...
plt.clf();  plt.close();
fig = plt.figure(figsize=(12,10)); fig.subplots_adjust(wspace=0.26, hspace=0.285)
fig.suptitle('Cloud stats:' + idate + '-- N-len=' + str(len(cldlen)) ,x=0.35, y=0.98,horizontalalignment='left',verticalalignment='top',fontsize=12)
ax1 = fig.add_subplot(241)
ax1.plot(Totmean_epsLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax1.plot(Totmedn_epsLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax1b = ax1.twiny()
ax1b.plot(np.round(100*(TotlenLR/len(cldlen))), nht, linewidth=1, color='b', marker='o', linestyle='--'); ax1b.set_xlabel('Sample (%)', size=12)
ax1.set_xlim(-5,-2); 
ax1.set_ylim(0.05,1.05); ax1.grid(True)
ax1.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=14);
ax1.set_ylabel('Normalized Height', size=14)
ax1.set_xticks(np.arange(-5,-1.9,1)); ax1b.set_xticks(np.arange(0,101,25));
plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax2 = fig.add_subplot(242)
ax2.plot(Totmean_dopLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax2.plot(Totmedn_dopLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax2.set_xlim(-1,1);ax2.set_ylim(0.05,1.05); ax2.grid(True)
ax2.set_xlabel('Doppler velocity [m s$\mathregular{^{-1}}$]', size=14);
plt.xticks(np.arange(-1,1.1,0.5), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax3 = fig.add_subplot(243)
ax3.plot(Totmean_sigLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax3.plot(Totmedn_sigLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax3.set_xlim(0,1.25); ax3.set_ylim(0.05,1.05); ax3.grid(True)
ax3.set_xlabel('w sigma [m$\mathregular{^{2}}$ s$\mathregular{^{-2}}$]', size=14);
plt.xticks(np.arange(0,1.251,0.25), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

ax4 = fig.add_subplot(244)
ax4.plot(Totmean_skwLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k')
ax4.plot(Totmedn_skwLR, nht, linewidth=2, color='k', marker='o', markeredgecolor='k', linestyle='--')
ax4.set_xlim(-1,1); ax4.set_ylim(0.05,1.05); ax4.grid(True)
ax4.set_xlabel('w skew [m$\mathregular{^{3}}$ s$\mathregular{^{-3}}$]', size=14);
plt.xticks(np.arange(-1,1.1,0.5), fontsize=12); plt.yticks(np.arange(0,1.09,0.1), fontsize=12);

plt.savefig(outdir + 'FIG1e_Cloudtotprofile_'+idate+'case.png', format='png', dpi=100)
del fig; del ax1; del ax2; del ax3; del ax4;

###f-plot
##Probability distribution functions...of wstats.........
wbins = np.arange(-2,2.01,0.1);
ebins = np.arange(-8,0.01,0.2);
kbins = np.arange(-5,5,0.1);
sbins = np.arange(0,2.01,0.1);
lwpdf, _ = np.histogram(ldopdata, wbins);  lwpdfn, _ = np.histogram(ldopndata, wbins);   rwpdf, _ = np.histogram(rdopdata, wbins);
lspdf, _ = np.histogram(lsigdata, sbins);  lspdfn, _ = np.histogram(lsigndata, sbins);   rspdf, _ = np.histogram(rsigdata, sbins);
lkpdf, _ = np.histogram(lskwdata, kbins);  lkpdfn, _ = np.histogram(lskwndata, kbins);   rkpdf, _ = np.histogram(rskwdata, kbins);
lepdf, _ = np.histogram(lepsdata, ebins);  lepdfn, _ = np.histogram(lepsndata, ebins);   repdf, _ = np.histogram(repsdata, ebins);
#plot
#
plt.close(); fig = plt.figure(figsize=(14,6)); fig.subplots_adjust(wspace=0.275, hspace=0.285)
fig.suptitle('Doppler + EPS + CLD stats:' + idate ,x=0.25, y=0.975,horizontalalignment='left',verticalalignment='top', fontsize=14)
ax1 = fig.add_subplot(241)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax1.plot(wbins[:-1], 100*(lwpdf/np.sum(lwpdf)), linewidth=2.0, color='b')
#ax1.plot(wbins[:-1], 100*(lwpdfn/np.sum(lwpdfn)), linewidth=2.0, color='b', linestyle='--')
ax1.plot(wbins[:-1], 100*(rwpdf/np.sum(rwpdf)), linewidth=2.0, color='r')
ax1.set_ylabel('Count [%]', size=12);   ax1.set_xlabel('w [m s$\mathregular{^{-1}}$]', size=12)
#plt.xticks(wbins[::8], fontsize=12);
plt.xticks(np.arange(-2,2.1,1), fontsize=12);
ax1.set_xlim(-2,2); ax1.grid(True);
plt.title('(a)', size=14, x=0.1, y=0.85)

ax2 = fig.add_subplot(242)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax2.plot(sbins[:-1], 100*(lspdf/np.sum(lspdf)), linewidth=2.0, color='b')
#ax2.plot(sbins[:-1], 100*(lspdfn/np.sum(lspdfn)), linewidth=2.0, color='b', linestyle='--')
ax2.plot(sbins[:-1], 100*(rspdf/np.sum(rspdf)), linewidth=2.0, color='r')
#ax2.set_ylabel('Count [%]', size=12);
ax2.set_xlabel('w-sig [m$\mathregular{^{2}}$ s$\mathregular{^{-2}}$]', size=12)
plt.xticks(sbins[::3], fontsize=12);
ax2.set_xlim(0,1.0); ax2.grid(True);
plt.title('(b)', size=14, x=0.1, y=0.85)

ax3 = fig.add_subplot(243)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax3.plot(kbins[:-1], 100*(lkpdf/np.sum(lkpdf)), linewidth=2.0, color='b')
#ax3.plot(kbins[:-1], 100*(lkpdfn/np.sum(lkpdfn)), linewidth=2.0, color='b', linestyle='--')
ax3.plot(kbins[:-1], 100*(rkpdf/np.sum(rkpdf)), linewidth=2.0, color='r')
#ax3.set_ylabel('Count [%]', size=12);
ax3.set_xlabel('w-skew [m$\mathregular{^{3}}$ s$\mathregular{^{-3}}$]', size=12)
#plt.xticks(kbins[::7], fontsize=12);
plt.xticks(np.arange(-2,2.1,1), fontsize=12);
ax3.set_xlim(-2,2); ax3.grid(True);
plt.title('(c)', size=14, x=0.1, y=0.85)

ax4 = fig.add_subplot(244)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax4.plot(ebins[:-1], 100*(lepdf/np.sum(lepdf)), linewidth=2.0, color='b')
#ax4.plot(ebins[:-1], 100*(lepdfn/np.sum(lepdfn)), linewidth=2.0, color='b', linestyle='--')
ax4.plot(ebins[:-1], 100*(repdf/np.sum(repdf)), linewidth=2.0, color='r')
#ax4.set_ylabel('Count [%]', size=12);
ax4.set_xlabel('log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]', size=12)
plt.xticks(ebins[::5], fontsize=12);
ax4.set_xlim(-8,0); ax4.grid(True);
plt.title('(d)', size=14, x=0.1, y=0.85)

#fig.savefig(outdir + idate + 'PDFwstats_lidarvsradar.png', format='png', dpi=100)
#del fig;
del ax1; del ax2; del ax3; del ax4;

##cloud base height, cloud top height, cloud depth, and liquid water path........
cbhbins = np.arange(0,2.01,0.1);  cthbins = np.arange(0,2.01,0.1);  cldbins = np.arange(0,2.01,0.1);
cbhpdf, _ = np.histogram(cbase, cbhbins);  cthpdf, _ = np.histogram(ctop, cthbins);  cldpdf, _ = np.histogram(ctop-cbase, cldbins)
lwpbins = np.arange(0,201,10); lwppdf, _ = np.histogram(slwp, lwpbins)
#
#plt.close(); fig = plt.figure(figsize=(8,8)); fig.subplots_adjust(wspace=0.25, hspace=0.25)
#fig.suptitle('CLOUD heights:' + idate ,x=0.4, y=0.96,horizontalalignment='left',verticalalignment='top', fontsize=14)
ax1 = fig.add_subplot(245)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax1.plot(cbhbins[:-1], 100*(cbhpdf/np.sum(cbhpdf)), linewidth=2.0, color='k')
ax1.set_ylabel('Count [%]', size=12);   ax1.set_xlabel('CBH [km]', size=12)
ax1.set_xlim(0,1.5); ax1.grid(True);
plt.xticks(np.arange(0,1.51,0.3), fontsize=12);
plt.title('(e)', size=14, x=0.1, y=0.85)

ax2 = fig.add_subplot(246)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax2.plot(cthbins[:-1], 100*(cthpdf/np.sum(cthpdf)), linewidth=2.0, color='k')
#ax2.set_ylabel('Count [%]', size=12);
ax2.set_xlabel('CTH [km]', size=12)
ax2.set_xlim(0,1.5); ax2.grid(True);
plt.xticks(np.arange(0,1.51,0.3), fontsize=12);
plt.title('(f)', size=14, x=0.1, y=0.85)

ax3 = fig.add_subplot(247)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax3.plot(cldbins[:-1], 100*(cldpdf/np.sum(cldpdf)), linewidth=2.0, color='k')
#ax3.set_ylabel('Count [%]', size=12);
ax3.set_xlabel('CLD [km]', size=12)
plt.xticks(np.arange(0,1.51,0.3), fontsize=12);
ax3.set_xlim(0,1.5); ax3.grid(True);
plt.title('(g)', size=14, x=0.1, y=0.85)

ax4 = fig.add_subplot(248)
plt.gcf().subplots_adjust(bottom=0.1, left=0.2)
ax4.plot(lwpbins[:-1], 100*(lwppdf/np.sum(lwppdf)), linewidth=2.0, color='k')
#ax4.set_ylabel('Count [%]', size=12);
ax4.set_xlabel('LWP [g. ms$\mathregular{^{-2}}$]', size=12)
plt.xticks(lwpbins[::5], fontsize=12);
ax4.set_xlim(0,200); ax4.grid(True);
plt.title('(h)', size=14, x=0.1, y=0.85)

fig.savefig(outdir + 'FIG1f_Clouddopepsstat_'+idate+'case.png', format='png', dpi=100)
del fig; del ax1; del ax2; del ax3; del ax4;
#'''

###g-plot
#converting the lidar into radar shape...........
linx=np.where(Lranges>=2.9)[0][0]
rinx=np.where(Rranges>=2.9)[0][0]
leps = LEPS[:,:linx];     lerr = LEPSerr[:,:linx]
varepsRnew = REPS[:,:rinx]
varerrRnew = REPSerr[:,:rinx]
lrng = Lranges[:linx];     rrng = Rranges[:rinx]
varepsLnew = np.asarray(np.zeros((len(atime), len(rrng))))
varerrLnew = np.asarray(np.zeros((len(atime), len(rrng))))
varepsLnew[:] = np.nan; varerrLnew[:] = np.nan
for r in range(len(rrng)-1):
  arng = rrng[r]; brng = rrng[r+1]
  inx = np.where((lrng>=arng) & (lrng<=brng))[0]
  for t in range(len(atime)):
     eps = leps[t,inx];         err = lerr[t,inx]
     varepsLnew[t,r] = eps[~np.isnan(eps)].mean()
     varerrLnew[t,r] = err[~np.isnan(err)].mean()
     del eps; del err;
  del arng; del brng; del inx;

print "samples lidar epsilon shape:", idate,'--', varepsLnew.shape, varepsRnew.shape
varrng =  rrng;
del linx; del rinx; del leps; del lerr; del lrng; del rrng;

###FIG-3.........###PART-C.........considering total time....avergaed for 1.5km.......range...each 2-minute time.....
histepsL=[];  histepsR=[]; histepsL1=[];  histepsR1=[];   #histeht=[];  histeht1=[];
histerrL=[];  histerrR=[]; histerrL1=[];  histerrR1=[];
histepsLc=[]; histepsRc=[];
hinx = np.where(varrng>1.5)[0][0]
cldtime = vartime[a1:a2+1]
tt=0;
#for t in range(len(vartime)):
for t in np.arange(a1,a2+1):
  leps = varepsLnew[t,:][:hinx]
  reps = varepsRnew[t,:][:hinx]
  lerr = varerrLnew[t,:][:hinx]
  rerr = varerrRnew[t,:][:hinx]
  #eps-err
  for k in range(len(leps)):
    eps1 = leps[k]
    eps2 = reps[k]
    err1 = lerr[k]
    err2 = rerr[k]
    if ((np.isnan(eps1)==False) and (np.isnan(eps2)==False)):
      tt=tt+1
      histepsLc.append(eps1)
      histepsRc.append(eps2)
      if ((eps1<0) and (eps2<0) and (eps1>-8) and (eps2>-8) and (np.abs(np.abs(eps1)-np.abs(eps2)) <=1)):
        histepsL.append(eps1)
        histepsR.append(eps2)
        histerrL.append(err1)
        histerrR.append(err2)
      elif ((eps1<0)and(eps2<0) and (eps1>-8)and(eps2>-8) and (np.abs(np.abs(eps1)-np.abs(eps2))>1) and (np.abs(np.abs(eps1)-np.abs(eps2))<=2.5)):
        histepsL1.append(eps1)
        histepsR1.append(eps2)
        histerrL1.append(err1)
        histerrR1.append(err2)

    del eps1; del eps2; del err1; del err2;
  del leps; del reps; del lerr; del rerr;

print "2d-frequency histrogram length of epsilons", len(histepsL), len(histepsR), len(histepsL1), len(histepsR1), '-FALSE-number-', tt
print "2d-frequency histrogram length of err--fed", len(histerrL), len(histerrR), len(histerrL1), len(histerrR1), '-FALSE-number-', tt

#Saving the converted equal height lidar and radar data.......
sio.savemat(outdir + 'EPSDOPTstats_scat_'+idate+'.mat', dict([('leps', varepsLnew[a1:a2+1,:]),('reps',varepsRnew[a1:a2+1,:]),('lerr', varerrLnew[a1:a2+1,:]),('rerr', varerrRnew[a1:a2+1,:]), ('atime', atime[a1:a2+1]), ('varrng', varrng), ('mwrdec', mwrdecoupled[a1:a2+1])]))
sys.exit()

#colorbar for scatter plot........
hist_levs = np.arange(0,51,1);   hist_cmap = matplotlib.cm.get_cmap('YlGn');
hist_norm = matplotlib.colors.BoundaryNorm(hist_levs,hist_cmap.N,clip=True)
#
corr=np.corrcoef(histepsL, histepsR)
from scipy import stats
slope, intercept, _, _, _ = stats.linregress(histepsL,histepsR)
line = slope*np.asarray(histepsL)+intercept
#
myleps  = np.arange(-8,0.1,0.1); myreps  = np.arange(-8,0.1,0.1)
histdata,_,_ = np.histogram2d(histepsL, histepsR, bins=[myleps, myreps])
histdata1,_,_ = np.histogram2d(histepsL1, histepsR1, bins=[myleps, myreps])
histdata2 = histdata+histdata1
Myleps, Myreps = np.meshgrid(myleps[:-1], myreps[:-1])
histdata[histdata<=2]=np.nan; histdata1[histdata1<=2]=np.nan; histdata2[histdata2<=2]=np.nan;

#FIGURE scatter...with freq2d........2subplots. 1) freq2d 2) scatter.diag......
plt.clf(); plt.close()
fig = plt.figure(figsize=(9,9));  fig.subplots_adjust(wspace=0.285, hspace=0.295)
ax1=fig.add_subplot(221)
cf1=plt.contourf(Myleps, Myreps, np.round(histdata2.T), levels=hist_levs, cmap=hist_cmap, norm=hist_norm)
z_min = 0; z_max = 50;
histarray = np.ma.masked_where(np.isnan(histdata2),histdata2)
cf1=plt.pcolor(Myleps, Myreps, np.round(histarray.T), cmap=hist_cmap, vmin=z_min, vmax=z_max)
plt.title('N = '+ str(tt)+' and R = ' + str("{0:.2f}".format(np.abs(corr[0,1]))),size=14)
ax1.plot(np.asarray(histepsL), line, '--g', label='line', linewidth=2, linestyle='--')
xax = np.arange(-8,-0.9,0.1); yax = np.arange(-8,-0.9,0.1); plt.plot(xax, yax, color='k', linewidth=1)
plt.ylim(-8,-1);  plt.xlim(-8,-1); ax1.grid(True)
plt.yticks(np.arange(-8,-0.9,1), fontsize=12); plt.xticks(np.arange(-8,-0.9,1), fontsize=12)
plt.xlabel('LIDAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=14);
plt.ylabel('RADAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=14)
cf1.set_clim(0,z_max); cbaxes = fig.add_axes([0.15, 0.475, 0.3, 0.015])
plt.colorbar(cf1, cax=cbaxes, orientation='horizontal')

ax2=fig.add_subplot(223)
cf2=plt.scatter(histepsL, histepsR, color='k', s=1)
cf2=plt.scatter(histepsL1, histepsR1, color='k', s=1)
#plt.title('N = '+ str(tt)+' and R = ' + str("{0:.2f}".format(np.abs(corr[0,1]))),size=14)
ax2.plot(np.asarray(histepsL), line, '--g', label='line', linewidth=2, linestyle='--')
xax = np.arange(-8,-0.9,0.1); yax = np.arange(-8,-0.9,0.1); plt.plot(xax, yax, color='k', linewidth=1)
plt.ylim(-8,-1);  plt.xlim(-8,-1); ax2.grid(True)
plt.yticks(np.arange(-8,-0.9,1), fontsize=12); plt.xticks(np.arange(-8,-0.9,1), fontsize=12)
plt.xlabel('LIDAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=14);
plt.ylabel('RADAR log epsilon [m$\mathregular{^{2}}$ s$\mathregular{^{-3}}$]',size=14)

plt.savefig(outdir + 'FIG1g_LIDARvsRADAReps_'+idate+'case.png', format='png', dpi=100)
del ax1; del ax2; del fig; del cf1; del cf2; del z_min; del z_max; 
sys.exit()
sys.exit()
sys.exit()
