import drizzlepac
from drizzlepac import astrodrizzle
from drizzlepac import tweakreg
from drizzlepac import tweakback
import pylab as p
import numpy
from numpy import ma, finfo, float32
from stsci.tools import teal
teal.unlearn('astrodrizzle')
teal.unlearn('tweakreg')
adriz = astrodrizzle.AstroDrizzle
twreg = tweakreg.TweakReg
twback = tweakback.tweakback
from utils import basic_params
import numpy as np
from astropy.io import fits

import shutil
import glob
import os
filt = np.array(['f125', 'f140', 'f150', 'f165'])
filt_num = np.array(['125', '140', '150', '165'])


def drizzle_SBC(params, params_gal, i):
    """Entry point for combining different extensions of FLT images 

    Args:
        params: (dict) Dictionary of 'basic' section of config file (common for all ULIRGs) 
        params_gal: (dict) Dictionary of galaxy parameters from config file
        i : (int) ULIRG number (It goes from 0 to 4)

    """

    gal_name = params_gal['name']
    dark = params_gal['dark_frames']
    dark = (dark.split(','))

    primary_dir = params['data_dir'] + gal_name + '/'
    for k in range(4):
        frames = params_gal[filt[k]]
        c = frames.split(',')
        images = []

        for j in c:
            '''
            if (j == '11e6q' or j == '21req' or j == '41eeq'):
                images.append(primary_dir + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'drk_v2' + '_allext_flt.fits')
            '''
            if j in dark:
                images.append(primary_dir + 'UV_ALIGN/' + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'drk' + '_allext_flt.fits')
            else:
                images.append(primary_dir + 'UV_ALIGN/' + 'jcmc' + j + '_F' + filt_num[k] + 'LP_' + 'sky' + '_allext_flt.fits')
        output = primary_dir + 'UV_DRZ/' + 'gal%s_UV_F%s' % (i + 1, filt_num[k])
        logfile = primary_dir + 'UV_DRZ/TXT/astrodrizzle_gal%s_%s.log' % (i + 1, k)
        adriz(input=images, output=output,
              runfile=logfile, configobj=params['sbc_config'], clean=True)


def main(config_file):
 
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)
        drizzle_SBC(params, params_gal, i)


if __name__ == '__main__':
    main(config_file)


'''main function for drizzing images
    Args:
        config_file: (str) name of comfig file

    
'''