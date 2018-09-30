
import argparse
from scipy import interpolate
import numpy as np
from astropy.io import fits as fio
from astropy.io import fits
from utils import basic_params
from utils import FLC_centers
from utils import UV_centers
import pyraf
from pyraf import iraf
#-----------------------------------------------------------------------------
# Main
#-----------------------------------------------------------------------------
from collections import OrderedDict
import pandas as pd
import os

from astropy.table import Table, Column, MaskedColumn
from utils import file_remove

import photutils
from photutils import create_matching_kernel
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import scipy
from photutils import TopHatWindow
from photutils import CosineBellWindow


window = CosineBellWindow(alpha=0.35)
from scipy.signal import convolve as scipy_convolve
window = TopHatWindow(0.35)
filt = ['775', '782']


def psf_rotate(psf, psf_rot, angle):
    file_remove(psf_rot)
    iraf.rotate(psf + '[0]', psf_rot, -angle, interpolant='nearest')
    iraf.rotate.unlearn()


def psf_cutting(psf, angle, psf_rot, psf_cut, size):
    hdu = fits.open(psf_rot)
    data = hdu[0].data
    hdu.close()
    c = (np.unravel_index(data.argmax(), data.shape))
    print (c)

    data_stamp = data[int(c[0] - size / 2):int(c[0] + size / 2), int(c[1] - size / 2):int(c[1] + size / 2)]

    psf_cut = psf_rot.replace('rotate', 'rotate_cut')

    hdu = fits.PrimaryHDU(data=data_stamp)
    header = hdu.header
    header["ANGLE"] = (angle, 'rotation angle')
    header["ORIGIN"] = (psf, 'main input file')
    header["ROTATE"] = (psf_rot, 'roatated file')
    hdu.writeto(psf_cut, overwrite=True)


def psf_matching(i, fil, input, psf165):
    ker = os.path.dirname(input) + '/' + "ker%s_gal%s_ref165.fits" % (fil, i + 1)
    data_ref = fits.getdata(psf165)
    data_psf = fits.getdata(input)
    kernel = create_matching_kernel(data_psf, data_ref)  # , window = window )
    fits.writeto(ker, data=kernel, overwrite=True)
    return ker


def main(config_file):

    for i in range(5):

        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)
        tab = Table.read(params['flt_files'], format='csv')
        gal_name = params_gal['name']
        primary_dir = params['data_dir'] + gal_name + '/'

        for j in range(2):
            for k in range(len(tab)):
                if tab['galaxy'][k] == gal_name and tab['filtername'][k] == 'F775W':
                    angle = float(tab['orientation'][k])
            print (angle)
            psf = params['hapsf_dir'] + 'PSF_%s_gal%s.fits' % (filt[j], i + 1)
            psf_rot = psf.replace('gal%s' % (i + 1), 'gal%s_rotate' % (i + 1))
            psf_cut = psf_rot.replace("rotate", "rotate_cut")
            psf_rotate(psf, psf_rot, angle)
            size = 130
            psf_cutting(psf, angle, psf_rot, psf_cut, size)
            psf165 = params['hapsf_dir'] + 'f165psf.fits'
            data = fits.getdata(psf165)
            print (data.shape, data.sum())

            hdu = fits.open(psf_cut)
            data = hdu[0].data
            hdu[0].data = data / data.sum()
            hdu.writeto(psf_cut, overwrite=True)
            hdu.close()
            print (data.shape, data.sum())
    for j in range(2):
        psf_cut = params['hapsf_dir'] + 'PSF_%s_gal%s_rotate_cut.fits' % (filt[j], i + 1)
        psfref = params['hapsf_dir'] + 'PSF_775_gal%s_rotate_cut.fits' % (i + 1)  # 'f165psf.fits'  #
        ker = psf_matching(i, filt[j], psf_cut, psfref)
        ker_rot = ker.replace('165', '165_rotate')
        # rotate kernels
        psf_rotate(ker, ker_rot, -angle)


if __name__ == '__main__':
    main(config_file)
