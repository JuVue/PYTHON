"""
Steps to read in Matlab struct files (saved as .mat) and associated functions to make them easier to use
==============================

"""
from __future__ import print_function
from scipy.io import loadmat
import numpy as np
import scipy.io as sio
from IPython import embed

def readMatlabStruct(filename):

    #### EXAMPLE OF USE:
    #### data = readMatlabStruct('../../jutta/UserReadyData/radiosondes/SondeData_h10int_V02.mat')

    #### for reference:
        #### find struct name with:
            #### sio.whosmat(filename)

    ### ----------------------------------
    ### Find struct name from .mat file using sio
    ### ----------------------------------
    dat = sio.whosmat(filename)
    ### ----------------------------------
    ### Extract out struct name
    ### ----------------------------------
    structname=[]
    for i in range(0,len(dat)):
        structname.append(dat[i][0])

    #### --------------------------------------------------------------------
    #### LOAD MATLAB FILE USING SCIPY
    #### --------------------------------------------------------------------
    print ('Reading in .mat file including struct...')
    dat = loadmat(filename)
    print ('')
    dout={}
    #### --------------------------------------------------------------------
    #### USE STRUCT_NAME TO DEFINE INTERMEDIATE STRUCT ARRAY
    #### --------------------------------------------------------------------
    for i in range(0,len(structname)):
        print ('Dealing with intermediate data assignments...')
        struct = dat[structname[i]]

        print ('')

        #### --------------------------------------------------------------------
        #### IDENTIFY DATA AS FIRST ENTRY IN INTERMEDIATE STRUCT
        #### --------------------------------------------------------------------
        a = struct[0,0]
        print (a.dtype.names)
            #### data.dtype:
                #### returns keys of dictionary (normal python dictionary access
                #### commands don't quite work...). MATLAB structs come back as
                #### numpy structured arrays.

        #### --------------------------------------------------------------------
        #### CHANGE NUMPY STRUCTURED ARRAY TO DICTIONARY FOR EASE OF USE
        ####
        #### --------------------------------------------------------------------
        # aa = a.astype(float)
        # b = {name:a[name].astype(float) for name in a.dtype.names}
        b = {}
        for name in a.dtype.names:
            # print (name)
            # print (a[name].dtype)
            if a[name].dtype == 'object':
                # print ('yes this is an object')
                continue
            # elif a[name].dtype == 'dict':
            #     b[name] = a[name]
            else:
                b[name] = a[name].astype(float)

        dout
        print ('Finished! :)')
        print ('Reading out ' + structname[i] + ' struct within .mat file')
        print ('')
        if  len(structname)>1:
            dout[structname[i]]={}
            dout[structname[i]]=b
        else:
            dout=b
    return dout     #### returns structured numpy array containing matlab struct
