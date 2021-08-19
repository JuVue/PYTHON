###script to convert npy to matlab files

from scipy.io import savemat
import numpy as np
import glob
import os

m_out_dir = '22_control_20180913T0000Z_qinit2-800m_rand-800m_thForcing-0000-0600_12hTim/'
monc_exp_dir= '/gws/nopw/j04/ncas_radar_vol1/jutta/MONC/output/'


npzFiles = glob.glob(monc_exp_dir + m_out_dir + '*.npy')
for f in npzFiles:
    fm = os.path.splitext(f)[0]+'.mat'
    d = np.load(f, allow_pickle=True).items()
    savemat(fm, d)
    print('generated ', fm, 'from', f)
