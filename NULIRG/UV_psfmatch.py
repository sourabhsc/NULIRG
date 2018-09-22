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
a = ['F125', 'F140', 'F150', 'F165', 'FR782N_NB', 'F775W_BB', 'CLEAR1L', 'FR716N']  # 7425 central
b = ['F125', 'F140', 'F150', 'F165', 'FR782N', 'F775W', 'CLEAR1L', 'FR716N']


def psf_match(params, params_gal, i, j):

    """Entry point for combining different extensions of FLT images 

    Args:
        params: (dict) Dictionary of 'basic' section of config file (common for all ULIRGs) 
        params_gal: (dict) Dictionary of galaxy parameters from config file
        i : (int) ULIRG number (It goes from 0 to 4)
        j : (int) Filter number (It goes from 0 to 3) and correspond to (F125, F140, F150, F165)

    """
    gal_name = params_gal['name']
    primary_dir = params['data_dir'] + gal_name + '/'
    UVPSF_DIR = params['uvpsf_dir']
    if j < 4:
        file_name = "%sgal%s_UV_%s_drz.fits" % (primary_dir + 'UV_DRZ/', i + 1, a[j])
        file_out = primary_dir + os.path.basename(file_name).replace("drz", "psfmatch")
        hdulist = fits.open(file_name)
        prihdr = hdulist[1].header  # the primary HDU header
        dat = hdulist[1].data
        ########## replacing nan with zeros #############
        where_are_NaNs1 = np.isnan(dat)
        dat[where_are_NaNs1] = 0.0

        file_ker = fits.open('%sker%s_ref165.fits' % (UVPSF_DIR, b[j].replace("F", "")))
        ker = file_ker[0].data
        ker1 = np.pad(ker, ((0, 1), (0, 1)), mode='constant')
        # plt.figure()
        # img = plt.imshow(ker1)
        # plt.show()
        print (" i m running slowly. Takes 10 minutes for each galaxy filter = %s" % (a[j]))

        # dat1 = convolve(dat, ker1, boundary='extend')
        dat1 = scipy_convolve(dat, ker1, mode='same')

        # !!!!Remember no convolution for error !!!

        hdulist[1].data = dat1
    else:  # gal1_HA_FR782N_NB_UV_align.fits
        if j == 4:
            file_name = "%sgal%s_HA_%s_UV_iraf_v2.fits" % (primary_dir, i + 1, a[j])
            file_out = file_name.replace("UV_iraf_v2", "psfmatch")
        else:
            file_name = "%sgal%s_HA_%s_UV_align_v2.fits" % (primary_dir, i + 1, a[j])
            file_out = file_name.replace("UV_align_v2", "psfmatch")

        hdulist = fits.open(file_name)

        dat = hdulist[0].data

        file_ker = fits.open('%sker%s_rotate_ref165_gal%s.fits' % (PSF_DIR, b[j], i + 1))
        ker = file_ker[0].data
        ker1 = np.pad(ker, ((0, 1), (0, 1)), mode='constant')

        print (" i m running slowly. Takes 10 minutes for each galaxy filter = %s" % (a[j]))

        # dat1 = convolve(dat, ker1, boundary='extend')
        dat1 = scipy_convolve(dat, ker1, mode='same')

        # !!!!Remember no convolution for error !!!

        hdulist[0].data = dat1
    '''
    hdulist[0].header["CREATOR"] = "SC March 2nd 2018"
    hdulist[0].header.insert("HISTORY", ("KERNEL", file_ker), after = True)# "kernel file use for convolution")
    hdulist[0].header.insert("HISTORY", ("REFPSF", "F165"), after = True)# "reference PSF for psfmatching"

    hdulist[0].header["COMMENTS"] = "PSF matched image for filter %s using scipy convolve with mode = 'same' "%(a[j], )
    hdulist[0].header.insert("HISTORY", ("WARNING", "Dont't forget to pad kernel\
     with zeros (np.pad(ker,((0,1),(0,1)),mode='constant')"), after = True)
    '''

    hdulist.writeto(file_out, overwrite=True, output_verify="ignore")


def main(config_file):

    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)

        for j in range(4):

            psf_match(params, params_gal, i, j)


if __name__ == '__main__':
    main(config_file)


