"""This code contains many small utility functions that are constantly used in ULIRG reduction steps.
The main funciton creates the location list for FLC images which are used for optical psfmatching
"""

import configparser
from configparser import ExtendedInterpolation
import numpy as np
from astropy.io import fits
from numpy import unravel_index
# import ULIRG_params as param
from astropy.wcs import WCS
from astropy.io import ascii
import pandas as pd
# basic call for reading config


# import os
# os.environ["iraf"] = "/iraf27/iraf"
# os.environ["IRAFARCH"] = "linux64"


import os
os.environ["iraf"] = "/home/sourabh/bin/miniconda3/envs/iraf27/iraf"
os.environ["IRAFARCH"] = "linux"
from matplotlib import patches

from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke
from astropy.table import Table, Column, MaskedColumn


def rectangle(left_corner_x, left_corner_y, x_size, y_size, color_val, ax2):
    """Function to create rectange on an image given corners and size of the image
    Args:
        left_corner_x: (float) left x corner of rectangle
        left_corner_y: (float) left y corner of rectangle
        x_size: (float) size of the box along x axis
        y_size: (float) size of the box along y axis
        color_val: (str) color of the box lines
        ax2: (str) plot axes where the rectangle will be inserted
    """
    ax2.add_patch(patches.Rectangle((left_corner_x, left_corner_y), x_size, y_size,
                                    fill=False,
                                    #linestyle='dotted',
                                    color=color_val,
                                    linewidth=2.0))


def basic_params(configfile, section, section_gal):
    """Function to parse out the config file values. There are 6 sections in the config file:-
    basic, NULIRG1, NULIRG2, NULIRG3, NULIRG4, NULIRG5. Basic section contains all the parameters that are
    common for all ULIRGs, while all other paramters are ULIRG specific written in their respective sections

    Args:
        configfile:(str) the config file for the ULIRG analysis (Default file is DIR/config/params_config.cfg)
        section: (str) 'basic' section which contains paramters commont for all ULIRGs
        section_gal: (str) ULIRG specific paramters
    Returns:
        dict : params- Dictionary of 'basic' section of config file (common for all ULIRGs) \n
        params_gal: (dict) Dictionary of galaxy parameters from config file\n

    """
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    options_gal = config.options(section_gal)

    params = {}
    params_gal = {}

    for option in options:
        params[option] = config.get(section, option)
    for option in options_gal:
        params_gal[option] = config.get(section_gal, option)

    return params, params_gal


def pad_with(vector, pad_width, iaxis, kwargs):

    pad_value = kwargs.get('padder', 0)
    if pad_width[0] != 0 and pad_width[1] != 0:
        vector[:pad_width[0]] = pad_value
        vector[-pad_width[1]:] = pad_value

    elif pad_width[0] == 0 and pad_width[1] != 0:
        vector[-pad_width[1]:] = pad_value
    elif pad_width[1] == 0 and pad_width[0] != 0:
        vector[:pad_width[0]] = pad_value

    return vector


def UV_centers(params, params_gal, i):
    """Function for finding center of drizzled F125 image used as refrence for optical image alignment.
    centers are found out by finding the pixel with largest flux value.
    Args:
        params: (dict) Dictionary of 'basic' section of config file (common for all ULIRGs)
        params_gal: (dict) Dictionary of galaxy parameters from config file
        i : index of ULIRG
    Returns:
        float: cent_x- center x for F125 image\n
        cent_y: center y\n
        cent_r: RA of the center\n
        cent_dec: Dec of the center\n

    """
    pos = []
    primary_dir = params['data_dir'] + params_gal['name'] + '/'
    hdulist = fits.open("%sgal%s_UV_F125_drz.fits" % (primary_dir + 'UV_DRZ/', i + 1))
    data = hdulist[1].data
    where_are_NaNs = np.isnan(data)

    data[where_are_NaNs] = 0
    if i == 2:
        data[:, 700:-1] = 0.0
        data[data > 0.020] = 0.0

    c = (unravel_index(data.argmax(), data.shape))
    pos.append((c[1], c[0]))

    header = hdulist[1].header
    w = WCS(header)
    r1 = w.wcs_pix2world(c[1], c[0], 0)
    cent_ra = (r1[0])  # changed integer thing
    cent_dec = (r1[1])

    cent_x = c[1]
    cent_y = c[0]

    return cent_x, cent_y, cent_ra, cent_dec


def UV_centers_drz(filename):
    data = fits.getdata(filename)
    where_are_NaNs = np.isnan(data)
    data[where_are_NaNs] = 0
    data[0:200, :] = 0.0  # removing hot pixels
    data[:, 800:] = 0.0  # removing hot pixels
    c = (unravel_index(data.argmax(), data.shape))
    xc = c[1]
    yc = c[0]
    # xcenter is c[1] and ycenter is c[0]
    return xc, yc


def DRC_centers(filename):
    hdu = fits.open(filename)
    data = fits.getdata(filename)
    where_are_NaNs = np.isnan(data)
    data[where_are_NaNs] = 0
    header = hdu[1].header
    data[0:200, :] = 0.0  # removing hot pixels
    data[:, 800:] = 0.0  # removing hot pixels

    c = (unravel_index(data.argmax(), data.shape))
    w = WCS(header)
    r1 = w.wcs_pix2world(c[1], c[0], 0)

    print(c[1], c[0])
    drc_ra = (r1[0])  # changed integer thing
    drc_dec = (r1[1])
    return drc_ra, drc_dec


def FLC_centers(i, cent_ra, cent_dec, file):
    """Function to highest flux position around the positon of galaxy in FLC files.

    These are used for psfmatching in optical filters.

    Args:
        i : (int) index of galaxy
        cent_ra:(float) RA of central pixel in F125 image
        cent_dec:(float) DEC of central pixel in F125 image
        file : (str) FLC filename
    Note:
        Remember for ULIRG 5 we have to choose a different extension (4) to find the FLC center. 
        Because ULIRG lies on chip-2 while for all others its on chip 1.
        Which is why extension 1 works for first 4 ULIRGs

    Returns:
        float: pixxx- central x\n
        pixy- central y \n

    """
    hdu = fits.open(file)
    if i == 4:
        header = hdu[4].header
    else:
        header = hdu[1].header

    # w = WCS.dropaxis(header)
    w = WCS(header)
    r1 = w.wcs_world2pix(float(cent_ra), float(cent_dec), 0)
    pixx = (r1[0])  # changed integer thing
    pixy = (r1[1])
    return pixx, pixy


def masks_circular(cent_x, cent_y, width, aper_lim, nx, ny):
    """Function for creating circular aperture given center of circle.

    Args:
        cent_x:(float) x center pixel
        cent_y:(float) y center pixel
        width: (float) width of each aperture
        aper_lim:(float) maximum radius of aperture
        nx:(int) x width of total ULIRG image
        ny:(int) y width of total ULIRG image

    Returns:
        dict: masks- circular masks\n
        masks_annulus - masks for annuli\n

    """
    rad1 = np.arange(1., aper_lim, width)
    y, x = np.mgrid[0:ny, 0:nx]
    masks_annulus = [np.where(((x - cent_x)**2 + (y - cent_y)**2 >= rad1[k]**2)
                              & ((x - cent_x)**2 + (y - cent_y)**2 <= rad1[k + 1]**2)) for k in range(len(rad1) - 1)]

    masks = [np.where((x - cent_x)**2 + (y - cent_y)**2 < rad1[k]**2) for k in range(len(rad1))]
    rad_annulus = ([(a + b) / 2 for a, b in zip(rad1, rad1[1:])])

    return rad1, rad_annulus, masks, masks_annulus


def mkdirp(directory):
    """makes adirectory if it does not exist

    """
    if not os.path.isdir(directory):
        os.makedirs(directory)


def file_remove(filename):
    """ removes a file if it exists

    """
    if os.path.exists(filename):
        os.remove(filename)


def circle(x, y, rad, col_circ1, ax4):
    """Plots a circle given center, radius, color and axes to plot. 

    """
    circle = Circle((x, y), rad, clip_on=False, linewidth=5.5, edgecolor=col_circ1, facecolor=(0, 0, 0, .0125), label="%s" % (rad))  # ,
    # path_effects=[withStroke(linewidth=5, foreground='w')])
    ax4.add_artist(circle)


gal_id_all = ["IRASF10594+3818", "IRASF12447+3721", "IRASF13469+5833", "IRASF14202+2615", "IRASF22206-2715"]
name_wfc = ['jcmc13prq_flc.fits', 'jcmc23grq_flc.fits', 'jcmc33dgq_flc.fits', 'jcmc43nqq_flc.fits', 'jcmc53tnq_flc.fits']


if __name__ == '__main__':
    filt = ['775', '782']
    xx = []
    yy = []
    simpos_ra = []
    simpos_dec = []
    for i in range(5):

        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params('/home/sourabh/ULIRG_package/config/params_default.cfg', 'basic', section_gal)
        primary_dir = params['data_dir'] + params_gal['name'] + '/'
        cent_x, cent_y, UV_ra, UV_dec = UV_centers(params, params_gal, i)
        file_drc = '/home/sourabh/ULIRG_package/data/%s/HA_INTER/FITS/gal%s_HA_F775W_scale_04_cut_drc.fits' % (gal_id_all[i], i + 1)
        drc_ra, drc_dec = DRC_centers(file_drc)
        name = '/home/sourabh/ULIRG_package/data/%s/HA_FLC/FITS/%s' % (gal_id_all[i], name_wfc[i])
        # print (cent_x, cent_y, UV_ra, UV_dec)
        print (drc_ra, drc_dec, UV_ra, UV_dec)
        pixx, pixy = FLC_centers(i, drc_ra, drc_dec, name)
        xx.append(pixx)
        yy.append(pixy)
        simpos_ra.append(drc_ra)
        simpos_dec.append(drc_dec)
        config_dir = params['config_dir']
        dict_cen_ha = dict.fromkeys(['xcen', 'ycen'])
    dict_cen_ha['xcen'] = xx
    dict_cen_ha['ycen'] = yy

    df = pd.DataFrame.from_dict(dict_cen_ha)
    df.to_csv(config_dir + 'ha_cent.txt', header=True, index=False, sep=' ', mode='w')

    print (xx, yy)
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params('/home/sourabh/ULIRG_package/config/params_default.cfg', 'basic', section_gal)
        primary_dir = params['data_dir'] + params_gal['name'] + '/'
        config_dir = params['config_dir']
        # txt_dir = primary_dir + 'HA_PSF/TXT/'

        for j in range(2):
            loclist = config_dir + 'loclist%s_gal%s.txt' % (filt[j], i + 1)
            simpos = config_dir + '/' + 'simpos%s_gal%s.txt' % (filt[j], i + 1)

            with open(loclist, 'w') as f:
                if i == 4:

                    [f.write("4 \t %s \t %s\n" % (xx[i], yy[i]))]
                    [f.write("4 \t %s \t %s\n" % (xx[i], yy[i]))]
                else:
                    [f.write("1 \t %s \t %s\n" % (xx[i], yy[i]))]
                    [f.write("1 \t %s \t %s\n" % (xx[i], yy[i]))]
            with open(simpos, 'w') as f:
                [f.write("%s \t %s \n" % (simpos_ra[i], simpos_dec[i]))]
                [f.write("%s \t %s \n" % (simpos_ra[i], simpos_dec[i]))]


# if __name__ == '__main__':

#   rad1, rad_annulus, masks, masks_annulus = masks_circular( 512, 512, 1, 510, 1024, 1024)
    # fits.writeto()


# os.path.basename
# b = np.pad(a, 2, 'constant', constant_values = 0)
# c = np.roll(b, 1, axis =1)  ### shift along y axis output to right if positive
# d = np.roll(b, 1, axis =0)  ### shift along y axis output to right if positive
