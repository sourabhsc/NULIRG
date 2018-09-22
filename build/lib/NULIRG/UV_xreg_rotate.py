import pyraf
from pyraf import iraf
import os
import scp
from astropy.io import fits
from astropy.table import Table, Column, MaskedColumn
import subprocess
# os.environ["iraf"]="/iraf/iraf"
# os.environ["IRAFARCH"]="linux64"

# my function
from utils import basic_params
from utils import file_remove

from stsci.tools import teal
#import iraf
#from iraf import pyraf


import glob


def iraf_rotate_command(key, rw, x_filter, primary_dir):

    """This function takes the name of file as input and returns the names from IRAF rotation output 

    Args:
        key: (str) key value for galaxy. e.g. jcmc11ctq_flt.fits has key == ctq
        rw: (str) name of FLT file
        x_filter: filter name . e.g. F165LP
        primary_dir: Direcoty of galaxy e.g. dir/IRASF10594+3818/

    Returns:
        str: (sky_dark_iraf) file name to be used for xregister. (just adding extenstion[1])\n
        rotate_out: output of roating the image\n
        sky_dark_err: ERR extension input image\n
        rotate_out_err: ERR extension of rotated image\n
        sky_dark_DQ: DQ extension of input image\n
        rotate_out_DQ: DQ extension of rotated image\n
    Examples:
        Examples should be written in doctest format, and should illustrate how
        to use the function.

        >>> print([i for i in example_generator(4)])
        [0, 1, 2, 3]

    """
    
    sky_dark_flt = rw.replace('flt.fits', '%s_flt.fits' % (key))
    sky_dark_iraf = '%s%s[1]' % (primary_dir + 'UV_DARK_SUB/', sky_dark_flt)
    sky_dark_err = '%s%s[2]' % (primary_dir + 'UV_DARK_SUB/', sky_dark_flt)
    sky_dark_DQ = '%s%s[3]' % (primary_dir + 'UV_DARK_SUB/', sky_dark_flt)

    rotate_out = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt.fits' % (x_filter, key))
    rotate_out = '%s%s' % (primary_dir + 'UV_ALIGN/FITS/', rotate_out)

    rotate_out_err = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt_err.fits' % (x_filter, key))
    rotate_out_err = '%s%s' % (primary_dir + 'UV_ALIGN/FITS/', rotate_out_err)

    rotate_out_DQ = sky_dark_flt.replace('%s_flt.fits' % (key), '%s_%s_rotate_flt_DQ.fits' % (x_filter, key))
    rotate_out_DQ = '%s%s' % (primary_dir + 'UV_ALIGN/FITS/', rotate_out_DQ)

    return sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ


def iraf_xreg_run(key, rw, x_filter, primary_dir, params_gal):

    """This function cross correlates each image with a reference image to align them 

    Args:
        key: (str) key value for galaxy. e.g. jcmc11ctq_flt.fits has key ==
        rw: (str) name of FLT file
        x_filter: (str)filter name . e.g. F165LP
        primary_dir: (str)Direcoty of galaxy e.g. dir/IRASF10594+3818/
        params_gal: (dict) Dictionary of galaxy parameters from config file


    Returns:
        str: (xreg_out) shifted image with a reference image

    """

    xreg_ref = "%sjcmc%s_sky_rotate_flt.fits[0]" % (primary_dir + 'UV_ALIGN/FITS/', params_gal["xreg_ref"])
    xreg_ref_err = xreg_ref.replace('rotate_flt', 'rotate_flt_err')
    xreg_lim = params_gal["xreg_lim"]
    xreg_xwindow = float(params_gal["xreg_xwindow"])
    xreg_ywindow = float(params_gal["xreg_ywindow"])

    sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ =\
        iraf_rotate_command(key, rw, x_filter, primary_dir)

    rotate_out_data = rotate_out.replace('.fits', '.fits[0]')
    rotate_out_err = rotate_out_err.replace('.fits', '.fits[0]')

    xreg_shift = rw.replace("flt.fits", '%s_%s_shift.txt' % (x_filter, key))
    xreg_shift = "%s%s" % (primary_dir + 'UV_ALIGN/TXT/', xreg_shift)
    xreg_out = rotate_out_data.replace("rotate_flt.fits[0]", "xregflt.fits")
    xreg_out_err = xreg_out.replace("xregflt.fits", "xregflt_err.fits")
    xreg_shift_err = xreg_shift.replace("shift.txt", "shift_err.txt")

    iraf.xregister.unlearn()
    #print (rotate_out_data, xreg_ref,  xreg_lim, xreg_shift, xreg_out, 'sawtooth', xreg_xwindow, xreg_ywindow, 'nearest')

    iraf.xregister(rotate_out_data, xreg_ref, xreg_lim, '%s' % (xreg_shift),
                   out='%s' % (xreg_out), function='sawtooth', xwindow=xreg_xwindow, ywindow=xreg_ywindow, interp_type='nearest')
    # iraf.wcscopy
   # iraf.wcscopy(file2 + "[2]", file_err + "[0]"  )  ### ref input
   # iraf.wcscopy(file2 + "[3]", file + "[0]"  )

    iraf.xregister.unlearn()
    iraf.xregister(rotate_out_err, xreg_ref, xreg_lim, '%s' % (xreg_shift_err),
                   out='%s' % (xreg_out_err), function='sawtooth', xwindow=xreg_xwindow, ywindow=xreg_ywindow, interp_type='nearest')
    iraf.xregister.unlearn()


def rotate(params, params_gal, dark_perform):
    """This function rotates the image 

    Args:
        key: (str) key value for galaxy. e.g. jcmc11ctq_flt.fits has key ==
        rw: (str) name of FLT file
        x_filter: (str)filter name . e.g. F165LP
        params: (dict) Dictionary of 'basic' section of config file (common for all ULIRGs) 
        params_gal: (dict) Dictionary of galaxy parameters from config file


    Returns:
        str : (xreg_out) shifted image with a reference image

    """

    flt_files = params['flt_files']
    tab = Table.read(flt_files, format='csv')
    rw = list(tab["filename"])

    gal_name = params_gal['name']
    dark = params_gal['dark_frames']
    dark = (dark.split(','))

    bad = params_gal['bad_frames']
    bad = (bad.split(','))

    primary_dir = params['data_dir'] + gal_name + '/'
    t = 0

    for i in range(len(tab)):

        c = []
        for letter in rw[i]:
            c.append(letter)
        gal_key = c[4] + c[5] + c[6] + c[7] + c[8]
        x_filter = tab["filtername"][i]

        if t == 0:
            ref_angle = float(tab["orientation"][i])
        rotate_ang = float(tab["orientation"][i]) - ref_angle
        if gal_key not in bad  \
                and gal_key not in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["detector"][i] == "SBC":
            t = t + 1

            sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ = \
                iraf_rotate_command("sky", rw[i], x_filter, primary_dir)

            iraf.rotate(sky_dark_iraf, rotate_out, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()
            iraf.rotate(sky_dark_err, rotate_out_err, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()
            iraf.rotate(sky_dark_DQ, rotate_out_DQ, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

        if gal_key in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["detector"][i] == "SBC" and dark_perform:

            sky_dark_iraf, rotate_out, sky_dark_err, rotate_out_err, sky_dark_DQ, rotate_out_DQ = \
                iraf_rotate_command("drk", rw[i], x_filter, primary_dir)

            iraf.rotate(sky_dark_iraf, rotate_out, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

            iraf.rotate(sky_dark_err, rotate_out_err, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

            iraf.rotate(sky_dark_DQ, rotate_out_DQ, rotate_ang, interpolant='nearest')
            iraf.rotate.unlearn()

    return tab, bad, dark, primary_dir


def xreg(params, params_gal, tab, bad, dark, primary_dir, dark_perform):
    """MAin function for applying xreguster 

    Args:
        params: (dict) Dictionary of 'basic' section of config file (common for all ULIRGs) 
        params_gal: (dict) Dictionary of galaxy parameters from config file
        tab : (dict) table with FLT file names 
        bad : (list) bad exposures list
        dark : (list) exposures that need dark subtraction
        primary_dir : (str) data_dir _ galaxy_name
        dark_perform : (bool) default - True

    Returns:
        str : shifted image with a reference image

    
    """

    rw = list(tab["filename"])
    gal_name = params_gal["name"]

    for i in range(len(tab)):
        c = []
        for letter in rw[i]:
            c.append(letter)
        gal_key = c[4] + c[5] + c[6] + c[7] + c[8]
        x_filter = tab["filtername"][i]

        if gal_key not in bad  \
                and gal_key not in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["detector"][i] == "SBC":

            iraf_xreg_run("sky", rw[i], x_filter, primary_dir, params_gal)

        if gal_key in dark  \
                and tab["galaxy"][i] == gal_name \
                and tab["detector"][i] == "SBC" and dark_perform:

            iraf_xreg_run("drk", rw[i], x_filter, primary_dir, params_gal)


def remove_files(key):
    """Function to remove all intermediate files from an xregister and rotate Run

    Args:
        key : (str) keyword for an exposure 
    Note:
        Make sure to remove these old files. Otherwise outputs from xregister \n
        just keep on adding additional extension to output file

    """
    f1 = glob.glob('%s*%s_rotate_flt.fits' % (primary_dir + 'UV_ALIGN/FITS/', key))
    f2 = glob.glob('%s*%s_rotate_flt_err.fits' % (primary_dir + 'UV_ALIGN/FITS/', key))
    f3 = glob.glob('%s*%s_rotate_flt_DQ.fits' % (primary_dir + 'UV_ALIGN/FITS/', key))

    f4 = glob.glob('%s*%s_xregflt.fits' % (primary_dir + 'UV_ALIGN/FITS/', key))
    f5 = glob.glob('%s*%s_xregflt_err.fits' % (primary_dir + 'UV_ALIGN/FITS/', key))

    f6 = glob.glob('%s*%s_shift.txt' % (primary_dir + 'UV_ALIGN/TXT/', key))
    f7 = glob.glob('%s*%s_shift_err.txt' % (primary_dir + 'UV_ALIGN/TXT/', key))

    [file_remove(a) for a in f1]
    [file_remove(a) for a in f2]
    [file_remove(a) for a in f3]
    [file_remove(a) for a in f4]
    [file_remove(a) for a in f5]
    [file_remove(a) for a in f6]
    [file_remove(a) for a in f7]


def main(config_file):

    first_run = True
    dark_perform = True
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))

        params, params_gal = basic_params(config_file, 'basic', section_gal)
        gal_name = params_gal['name']
        global primary_dir
        primary_dir = params['data_dir'] + gal_name + '/'
        if first_run:
            remove_files('sky')

        if dark_perform:
            remove_files('drk')
        tab, bad, dark, primary_dir = rotate(params, params_gal, dark_perform)
        xreg(params, params_gal, tab, bad, dark, primary_dir, dark_perform)


if __name__ == '__main__':

    main(config_file)
