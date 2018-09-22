import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import scipy.ndimage.filters as g
from scipy.interpolate import interp1d
import sys
from scipy.signal import convolve as scipy_convolve
import os
from astropy.convolution import convolve

from utils import basic_params
a = ['HA_cont', 'HA']
filt = ['775', '782']


def psf_match(filename, fileout, ker):
    hdulist = fits.open(filename)
    data = hdulist[0].data

    hdu_ker = fits.open(ker)
    ker_data = hdu_ker[0].data
    ker_shift = np.pad(ker_data, ((0, 1), (0, 1)), mode='constant')

    # dat1 = convolve(dat, ker1, boundary='extend')

    data_out = scipy_convolve(data, ker_shift, mode='same')
    hdulist[0].data = data_out

    hdulist.writeto(fileout, overwrite=True, output_verify="ignore")


def main(config_file):
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)
        gal_name = params_gal['name']
        primary_dir = params['data_dir'] + gal_name + '/'

        for j in range(2):
            filename = primary_dir + "gal%s_%s.fits" % (i + 1, a[j])
            fileout = filename.replace('.fits', '_psfmatch.fits')
            ker = "ker%s_gal%s_ref165_rotate.fits" % (filt[j], i + 1)
            print (" i m running slowly. Takes 10 minutes for each galaxy filter = %s \n" % (filt[j]))
            print (filename, ker)
            psf_match(filename, fileout, ker)

        print ("convolution data done galaxy %s \n" % (i + 1))


if __name__ == '__main__':
    main(config_file)
