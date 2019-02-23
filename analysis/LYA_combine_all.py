from __future__ import division, print_function, absolute_import
import sys
sys.path.append("/home/sourabh/ULIRG_package/NULIRG/") # go to parent dir
config_file = "/home/sourabh/ULIRG_package/config/ULIRG_params.cfg"
from astropy.io import fits
from matplotlib import pyplot as plt
import utils
from utils import UV_centers_drz
from utils import masks_circular
from utils import basic_params


from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from astropy.visualization import (MinMaxInterval, AsinhStretch,
								   ImageNormalize)
from scipy.ndimage.filters import gaussian_filter

import glob
import numpy as np
from matplotlib import pylab
plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
from utils import circle
import matplotlib
from collections import OrderedDict
matplotlib.__version__
import pandas as pd
from numpy import linalg as LA

def plot_rectangle(x,y, szx, szy, color, ax):
    utils.rectangle(x, y, szx, szy, color, ax)
def sum_ext(data, data_dir, i, ext ):
    files_all = glob.glob('%s/ULIRG%s_pos_m/*%s*'%(data_dir , i+1, ext ))
    s = [fits.getdata(f) for f in files_all]
    sum_all = np.zeros((data.shape[0], data.shape[1]))
    for s1 in s:
        sum_all = sum_all + s1
    return sum_all


data_dir = '/home/sourabh/ULIRG_package/analysis/lya_fit_data/'

def main():
    for i in range(5):

        file125 = data_dir + 'gal%s_UV_F125_scale_04_psfmatch.fits'%(i+1)
        hd = fits.open(file125)
        photflam = float(hd[1].header['PHOTFLAM'])
        data_125 = hd[1].data
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)
        cent_x, cent_y = UV_centers_drz(file125)


        print ('=============================================\n')
        print ('<<<<<<<<<<<<<<<<<<<<<ULIRG %s>>>>>>>>>>>>>>>>>>\n'%(i+1))
        print ('=============================================\n')

        print ('Original shape', data_125.shape)
        print ('UV 125 centers --> (%s, %s)'% (cent_x, cent_y))

    
        for extension in ext_list:

            data_ext = sum_ext(data_125, data_dir, i, extension) #'A_lmfit_x'
            fits.writeto('%s/ULIRG%s_pos_m/%s%s_no_mask.fits'%(extension, data_dir, i+1, i+1), data = data_lya, overwrite = True)
        
if __name__ == '__main__':
    main()