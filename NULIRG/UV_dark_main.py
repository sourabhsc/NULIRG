import os
import scp
from astropy.table import Table
import numpy as np
import matplotlib.pylab as pylab
from matplotlib import pyplot as plt
import configparser
from configparser import  ExtendedInterpolation
from astropy.io import ascii
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as patches
from subprocess import call
from astropy.table import Table, Column, MaskedColumn
from matplotlib.ticker import FormatStrFormatter
import glob
from scipy.ndimage.filters import gaussian_filter
import subprocess
import shutil
from termcolor import cprint 
from pyfiglet import figlet_format
from acstools import calacs
from subprocess import DEVNULL
from stsci.tools import teal
#from utilities_function import *
from utils import *

plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
import sys

def dark_from_ISR (configfile, section):
    """
    This module collects all the 'RAW' and 'FLT' dark files.

    jhhjhjhj.

    Parameters
    ----------
    configfile : str
        cofiguration file with information about all the important parameters used in ULIRG pipeline (e.g. ULIRG_params.cfg)
    section : str
        section of config file with basic info of file locations. (e.g. 'basic').

    
    Returns
    -------
    dark_RAW : list of str
        raw files.
    dark_FLT : list of str
        flt files.
    temp : list of str
        temp of all dark flt files.
    """
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    params = {}
    for option in options:
        params[option] = config.get(section, option)
    dark_dir = params['dark_dir']
    print (dark_dir)

    dark_RAW = glob.glob("%s*raw.fits"%(dark_dir))   
    ### DONOT FORGET GLOB CAN HAVE DIFFERENCT sequence of files in raw and flt
    dark_FLT = glob.glob("%s*flt.fits"%(dark_dir))
    temp_RAW = np.zeros(len(dark_RAW))
    temp_FLT = np.zeros(len(dark_RAW))

    for i in range(len(dark_RAW)):
        fits_dark_FLT =  fits.open(dark_FLT[i])
        fits_dark_RAW =  fits.open(dark_RAW[i])
        temp_FLT[i] = (float(fits_dark_FLT[1].header["MDECODT1"]) + float(fits_dark_FLT[1].header["MDECODT2"]))/2.
        temp_RAW[i] = (float(fits_dark_RAW[1].header["MDECODT1"]) + float(fits_dark_RAW[1].header["MDECODT2"]))/2.

    arg = np.argsort(temp_FLT)
    arg1 = np.argsort(temp_RAW)

    dark_RAW = list(dark_RAW[u] for u in arg1)
    dark_FLT = list(dark_FLT[u] for u in arg)
    temp_FLT = list(temp_FLT[u] for u in arg)
    #temp_RAW = list(temp_RAW[u] for u in arg)

    temp = temp_FLT


    return dark_RAW, dark_FLT, temp#, temp_FLT

def sky_dark_sub(gal_num, configfile, section, section_gal, show_output, dark_RAW, dark_FLT, temp , dark_perform):
    """
    This the main code for sky subtraction and dark subtraction.

    Extended description of function.

    Parameters
    ----------
    gal_num : int
        galaxy index (0,1,2,3,4)
    config_file : str
        the configuration file that contains most of the parameters required for ULIRG codes (e.g. ULIRG_params.conf)
    section : str
        the section of config file with basic parameters for the runs
        (e.g. 'basic')
    section_gal : str
        the section of the config file for a given galaxy 
        (e.g. 'NULIRG1')
    show_output : bool
        True if we want to see output in terminal
        False if not (default 'True') 
    dark_RAW : str
        list of all raw ISR dark files 
    dark_FLT : str
        list of all FLT ISR dark files
    temp : str
        list of temperatures corresponding to each dark
    dark_perform : bool
        True if dark subtraction is performed, False if not

    Returns 
    -------
    table_sky : table
    table_dark : table
    seg_map : arr
    smooth_gal : arr
    masks : arr
    mask_flt : arr
    output : sky subtracted image, dark subtracted image


    """
    config = configparser.ConfigParser(interpolation=ExtendedInterpolation())
    config.read(configfile)

    options = config.options(section)
    options_gal = config.options(section_gal)

    params = {}
    params_gal = {}
    # this sets all the parameter values from your config file
    
    for option in options:
        params[option] = config.get(section, option)
    for option in options_gal:
        params_gal[option] = config.get(section_gal, option)

    
    flt_files = params['flt_files']
    tab = Table.read(flt_files, format = 'csv')
    rw = list(tab["filename"])

    gal_name = params_gal['name']
    primary_dir = params['data_dir']+gal_name+ '/'

    dark = params_gal['dark_frames']
    dark = (dark.split(','))
    bad = params_gal['bad_frames']
    bad = (bad.split(','))
    global hot_pix_x, hot_pix_y

    hot_pix_x = np.array(params['hot_pix_x'].split(','))
    hot_pix_y = np.array(params['hot_pix_y'].split(','))
    print ( hot_pix_x[0],  hot_pix_y[0])
    print ( hot_pix_x[1],  hot_pix_y[1])

    print ("<<<<<Performing sky subtraction using NULIRG %s>>>>" %(gal_num+1))

    if show_output == "True":
        print ("Listing out bad frames")
        print (bad)
        print("Listing out hot frames T>25 that need separate dark subtraction")
        print (dark)



    sky_value = np.zeros(len(tab))
    table_sky = Table( names=('file_name', 'sky_value', 'exp_time', 'filter'), dtype = ('S100', 'f4', 'f4', 'S4'))
    table_sky.meta['comments'] = ['NULIRG %s  with name %s Output sky values calculated on FLT images for future preference \
    \n !Remember the corrresponding plots have sky value per exposure time '%(gal_num+1, gal_name)]

    table_dark = Table( names=('file_name', 'm1', 'dark_radii', 'ind', 'A_min', 'K_min', 'diff_ar', 'exp_time'), dtype = ('S100', 'f4', 'f4', 'f4', 'f4', 'f4','f4', 'f4'))
    table_dark.meta['comments'] = ['NULIRG %s  with name %s Output sky values calculated on FLT images for future preference \
    \n !Remember the corrresponding plots have sky value per exposure time '%(gal_num+1, gal_name)]

    rad1, rad_annulus, masks , masks_annulus = masks_circular (int(params_gal["cent_x"]),\
    int(params_gal["cent_x"]),\
    1.0,\
    int(params_gal["aper_lim"]),\
    int(params["nx"]),\
    int(params["ny"]))


    for i in range(len(tab)):

        c = []
        for letter in rw[i]:
            c.append(letter)

        gal_key = c[4]+c[5]+c[6]+c[7]+c[8]
        
        if gal_key not in bad  \
        and tab["galaxy"][i] == gal_name \
        and tab["detector"][i] =="SBC":
            file_name = primary_dir + rw[i]
            hdulist = fits.open(os.path.dirname(file_name)+'/UV_FLT/FITS/'+ os.path.basename(file_name))
            DQ = hdulist[3].data

            DQ[int(hot_pix_x[0]), int(hot_pix_y[0])] = float(params['dq'])
            DQ[int(hot_pix_x[1]), int(hot_pix_y[1])] = float(params['dq'])

            data  = hdulist[1].data
            data[DQ!=0] = float(params['dq_bad'])


            x_filter = tab["filtername"][i].replace("LP","")
            x_filter = x_filter.replace("F", "")
            x_filter = float(x_filter)
            aper_annulus = [np.mean(data[masks_annulus[k]]) for k in range(len(rad1)-1) ] 
            exp_time = hdulist[0].header["EXPTIME"]

            sky_value[i] = np.mean (aper_annulus[int(params_gal['sky_min']):int(params_gal['sky_max'])]) ##380:450
           
            print ("output plot for index = %s with keyword = %s filter = %s sky_value = %s"%(i, gal_key, x_filter, sky_value[i]))
            if show_output == "True":
                print ("output plot for index = %s with keyword = %s filter = %s sky_value = %s"%(i, gal_key, x_filter, sky_value[i]))

                #plt.show()
            ###### creating intermediate directories ####
            directory1 = os.path.dirname(primary_dir+"UV_DARK_SUB/PNG")
            directory2 = os.path.dirname(primary_dir+"UV_DARK_SUB/TXT")
            directory3 = os.path.dirname(primary_dir+"UV_DARK_SUB/FITS")

            if not os.path.exists(directory1):
                os.makedirs(directory1) 
            if not os.path.exists(directory2):
                os.makedirs(directory2) 
            if not os.path.exists(directory3):
                os.makedirs(directory3) 
            #####
            hdulist.close()
            hdulist = fits.open(os.path.dirname(file_name)+'/UV_FLT/FITS/'+os.path.basename(file_name))
            data = hdulist[1].data
            data[DQ!=0] = float(params['dq_bad'])
            hdulist[1].data = data - sky_value[i]
            hdulist[1].data[DQ!=0] = float(params['dq_bad'])

   

            table_sky.add_row((file_name.replace(primary_dir, ""), sky_value[i], exp_time, tab["filtername"][i]))
            table_name = '%sUV_DARK_SUB/TXT/FLT_sky_gal_%s.txt'%(primary_dir, gal_num+1)

            table_sky.write(table_name, format='ascii', overwrite = True)
            
            sky_sub_name = file_name.replace("flt.fits", "sky_flt.fits")
            sky_sub_name = primary_dir + 'UV_DARK_SUB/' + os.path.basename(sky_sub_name)
            hdulist.writeto(sky_sub_name, overwrite = True, output_verify="ignore")
            hdulist.close()
       
        if gal_key in dark  \
        and tab["galaxy"][i] == gal_name \
        and tab["detector"][i] == "SBC" and dark_perform:


            file_name = primary_dir + rw[i]
            print ("BHAAAAAAAAAAGOOOOOOOOOOOO>>>>DRK ENCOUNTERED>>>>>>")
            cprint(figlet_format('DARK', font='graffiti'),
       'white', attrs=['bold'])
            '''
            ____           ___
            |   \    /\    |  \  | / 
            |    |  /__\   | _|  |/
            |    | /    \  | \   |\
            |___/ /      \ |  \  | \
            '''

            print ("performing dark subtraction for ", os.path.basename(file_name).replace("flt", "raw"))

            print ("<<<<< smoothing FLT for sextractor segmentation maps >>>>")
            smooth_scale = float(params_gal['smooth_scale'])
            print (smooth_scale)
            sex_config = (params['sex_config'])

           
            file_smooth = smooth_gal(file_name, gal_key, smooth_scale, params )
            print ("<<<<< Getting SEXTRACTOR for finding segmentaion maps >>>>")

            seg_map = sextractor_seg(file_smooth, sex_config, params_gal)
            print ("<<<<<  Getting masked image for galaxy >>>>")

            file_mask = galaxy_mask(file_smooth, gal_key, smooth_scale, gal_num, seg_map, sex_config)


            file_name_flt = file_name#.replace("flt", "raw")
            #for m1 in range (-150, 50, 10):
            m1 = 0
            dark_sub(file_name_flt, file_mask, dark_FLT, dark_RAW, temp, params, params_gal, gal_num, m1, rad1, rad_annulus, masks, masks_annulus, table_dark )
            if show_output == "True":
                print ("dark subtraction done for ", file_name_flt, m1)
                print ("HHHHAAAAAPPPPYYY!!!")
       


def smooth_gal (file_name, gal_key, smooth_scale, params):
    """
    code for smoothing galaxy image for making ssextractor segemenatations

    Parameters
    ----------
    file_name : str
        input file
    gal_key : str
        input file key (e.g. for jcmc11ctq_flt.fits, gal_key = 'ctq')
    smooth_scale : str
        smoothing scale as provided in 'ULIRG_params.cfg'
    params: dict
        dictionary of all parameters form ULIRG_params.cfg

    Returns
    -------
    file_smooth : str
        name of the output smooth file 

    """
    hdulist =fits.open(os.path.dirname(file_name)+'/UV_FLT/FITS/'+os.path.basename(file_name))
    data = gaussian_filter(hdulist[1].data, float(smooth_scale))
    DQ = hdulist[3].data

    DQ[int(hot_pix_x[0]), int(hot_pix_y[0])] = float(params['dq'])
    DQ[int(hot_pix_x[1]), int(hot_pix_y[1])] = float(params['dq'])
    file_smooth = file_name.replace("flt.fits", "flt_smooth.fits")
    file_smooth = os.path.dirname(file_smooth) + '/UV_DARK_SUB/FITS/'+os.path.basename(file_smooth)

    data[DQ!=0] = float(params['dq_bad'])
    hdu = fits.PrimaryHDU(data=data)

    header = hdu.header
    header.add_history("smooth galaxy image to be used to create sextrator segmentation maps for %s"%(file_name))
    header['ROOTRAW'] = (gal_key, 'keyword for raw images')
    header['SMTHSCL'] = (smooth_scale, 'smoothing scale for gausssian filter')
    header['CREATOR'] = ('SSC', 'FEB 24 2018')
    hdu.writeto(file_smooth, overwrite = True)

    return file_smooth

def sextractor_seg(file_smooth, sex_config, params_gal):
    """
    code for creating segmentaion maps using sextractor

    Parameters
    ----------
    file_smooth : str
        smoothed  galaxy image
    sex_config : str
        default sextractor config file (e.g. sex_default.conf) 
    params_gal: dict
        dictionary of all parameters for a given ULIRG

    Returns
    -------
    seg_map : str
        name of the segmentation map FITS file 

    """
    catalog = file_smooth.replace("smooth.fits", "catalog.cat")
    catalog = catalog.replace('FITS', 'TXT')
    seg_map = file_smooth.replace("smooth.fits", "seg_map.fits")

    DETECT_THRESH = float(params_gal['detect_thresh'])
    DETECT_MINAREA = float(params_gal['detect_minarea'])
    ANALYSIS_THRESH = float(params_gal['analysis_thresh'])

    cmd = "sextractor %s \
    -c %s \
    -CATALOG_NAME %s             \
    -CHECKIMAGE_TYPE SEGMENTATION  \
    -CHECKIMAGE_NAME %s\
    -DETECT_MINAREA   %s\
    -DETECT_THRESH    %s\
    -ANALYSIS_THRESH  %s"\
    %(file_smooth, sex_config, catalog, seg_map, DETECT_MINAREA, DETECT_THRESH,  ANALYSIS_THRESH )
    print (cmd)
    call(cmd, shell = True)
    return seg_map

def galaxy_mask(file_smooth, gal_key, smooth_scale, gal_num, seg_map, sex_config ):
    """
    code for making a mask for galaxy image using ssextractor segemenatation maps

    Parameters
    ----------
    file_smooth : str
        smooth galaxy image
    gal_key : str
        input file key (e.g. for jcmc11ctq_flt.fits, gal_key = 'ctq')
    smooth_scale : str
        smoothing scale as provided in 'ULIRG_params.cfg'
    gal_num : int
        galaxy index (0,1,2,3,4)
    seg_map : str
        name of segmentation map FITS file
    sex_config: str
        default sextractor config file
    Returns
    -------
    file_mask : str
        name of the mask file 

    """

    hdu = fits.open(file_smooth)
    data = hdu[0].data
    seg = fits.open(seg_map)
    masks = np.where (seg[0].data==0)
    data[masks] = 0.0
    ### removing additional artifacts from corners etc.
    #file_mask = file_smooth.replace("smooth", "mask1")

    for j in range(1024): 
        for k in range(1024):
            if k<100 or j<100:
                data[j][k] =0.0
            if gal_num!=3 and gal_num!= 4:
                if j>650:
                    data[j][k] =0.0
    file_mask = file_smooth.replace("smooth", "mask")
    #hdu = fits.PrimaryHDU(data=data)
    header = hdu[0].header
    header['COMMENT'] = (" Image creation steps :- 1) smooth FLT image using gaussian filter of a smoothing scale %s,\
     2) use segmentation maps = %s, created by sextrator using config file = %s \
     3) Replace flux values in smooth FLT images at all the pixels with zeros in segmentaion map with zero\
     4) Also replace all pixels with j, k such that k<100 or j<100 set, data[j][k]==0.0 \
     5) ### j, k are opposite to as it appears on ds9 window\
     6) output = %s"%(smooth_scale, seg_map, sex_config, file_mask))
    
    header['SEXCFG'] = (sex_config, 'sextractor config file')
    hdu.writeto(file_mask, overwrite = True)

    return file_mask

def function_dark_sub(a,k, aper_diff_gal, aper_diff_dark, dark_radii, rad1, galaxy, dark, masks_annulus):
    """
    code for calculating minimizer value for a given set of parameters

    Parameters
    ----------
    a : float
        scaling parameter for dark
    k : float
        parameter for sky in dark fitting method
    aper_diff_gal : arr
        annuli mean for circular aperture  in galaxy image
    aper_diff_dark : arr
         annuli mean for circular aperture in dark images
    dark_radii : float
         radii to mask the galaxy (if using circular mask)
    rad1 : arr 
        list of radii
    galaxy : 2d arr
        galaxy imahe
    dark : 2d arr
        dark image
    masks_annulus: 2d arr
        annuli mask pixels


    Returns
    -------
    min_value : minimizer value

    """
    aper_diff_ann = np.array(aper_diff_gal) - a* np.array(aper_diff_dark) - k
    min_value = np.nansum(abs(aper_diff_ann[dark_radii:len(rad1)-1]))
    return min_value
def dark_sub(file_name_flt, file_mask, dark_FLT, \
    dark_RAW, temp, params, params_gal, gal_num, m1,\
     rad1, rad_annulus, masks, masks_annulus, table_dark):
    """
    main module for dark subtraction

    Parameters
    ----------
    file_name_flt : str
    file_mask
    dark_FLT
    dark_RAW
    temp
    params
    params_gal
    gal_num
    m1
    rad1 
    rad_annulus
    masks
    masks_annulus
    table_dark

    Returns
    -------
    minimizer 
    verification
    variance_fits
    dark_data 
    """
    gal_name = params_gal['name']
    primary_dir = params['data_dir']+gal_name+ '/'


    file_name_no_dir = file_name_flt.split("/")
    file_name_no_dir = file_name_no_dir[-1]
    hdulist = fits.open(os.path.dirname(file_name_flt)+'/UV_FLT/FITS/'+ os.path.basename(file_name_flt))
    exp_time = hdulist[0].header["EXPTIME"]
    data = hdulist[1].data

    fits_mask = fits.open(file_mask)
    data_mask  = fits_mask[0].data
    mask = np.where(data_mask!=0)
    data[mask] = 0.0

    hdu = fits.PrimaryHDU(data=data)
    header = hdu.header
    header.add_history("masked image with the mask created by sextractor usng file %s"%(file_mask))
    file_masked_flt = os.path.dirname(file_name_flt) +'/UV_DARK_SUB/FITS/'+ os.path.basename(file_name_flt)
    hdu.writeto(file_masked_flt, overwrite= True)
    


    temp_gal = (float(hdulist[1].header["MDECODT1"]) + float(hdulist[1].header["MDECODT2"]))/2.
    dark_radii=int(params_gal["dark_radii"])-m1
        
    A_lim = float(params["a_lim"])
    K_lim = float(params["k_lim"])
    del_A = float(params["del_a"])
    del_K = float(params["del_k"])
    

    A = np.arange(0, A_lim, del_A)
    K = np.arange(0, K_lim, del_K)
    A_min = np.zeros(len(dark_RAW))
    K_min = np.zeros(len(dark_RAW))
    minimum_var = np.zeros(len(dark_RAW))

    diff_ar = np.zeros((len(A), len(K)))


    for i in range(len(dark_RAW)):
        fits_dark_raw = fits.open(dark_RAW[i])
        ## scaling for difference between dark and galaxy exposure time difffference

        data_dark = fits_dark_raw[1].data*exp_time/float(params["exp_dark"]) ### exp time for darks is 1000 secs
    
        ### looking for DQ array from dark FLT files ###

        fits_dark_FLT =  fits.open(dark_FLT[i])
        data_dark_DQ = fits_dark_FLT[3].data
        

        DQ = fits_dark_FLT[3].data

        DQ[int(hot_pix_x[0]), int(hot_pix_y[0])] = float(params['dq'])
        DQ[int(hot_pix_x[1]), int(hot_pix_y[1])] = float(params['dq'])
        data_dark[mask] = 0.0
        
        data_dark[DQ!=0] = float(params['dq_bad'])
        data[DQ!=0] = float(params['dq_bad'])

        data[DQ!=0] = float(params['dq_bad'])
        aper_diff_gal =  [(np.mean(data[masks_annulus[k]])) for k in range(len(rad1)-1) ] 

        print ("performing minimization  for dark = %s" %( i+1))
        #### minimization steps###

        aper_diff_dark =  [(np.mean(data_dark[masks_annulus[k]])) for k in range(len(rad1)-1) ] 

        diff_ar = [function_dark_sub(a,k, aper_diff_gal, aper_diff_dark, dark_radii, rad1, data, data_dark, masks_annulus)  for a in A  for k in K]
        c1 = [(a,k)  for a in A  for k in K]

        c = (np.unravel_index(np.array(diff_ar).argmin(), np.array(diff_ar).shape))

        par = (np.array(c1)[c])
        print ("loop done")
        scale_factor = par[0]#A[c[0]]
        sky = par[1]#[c[1]]
        dark_final = scale_factor*data_dark +sky
        fits_dark_name = file_name_no_dir.replace("flt.fits", "dark_%s.fits"%(i+1))

      
        hdu = fits.PrimaryHDU(data=dark_final)
        header = hdu.header
        header["RAWNAME"] = (file_name_flt, "FLT file for dark subtraction")
        header["MASK"] = (file_mask, "mask file used for galaxy")
        header["AMIN"] = (par[0], " A value that minimizes the spatial variation" )
        header["KMIN"] = (par[1], "K value that minimizes the spatial variation")
        header["DARKFILE"] = ( dark_RAW[i], "dark file")
        header["ALIM"] = A_lim
        header["KLIM"] = K_lim
        header["ADEL"] = del_A
        header["KDEL"] = del_K
        header["MINVAR"] = np.min(diff_ar)
        header.add_history(" dark file dark subtraction method using G_subtracted = Galaxy - A*Dark - K. THis file has (A*Dark+K) ")
        hdu.writeto("%sUV_DARK_SUB/FITS/%s"%(primary_dir, fits_dark_name), overwrite= True)
        
        print ("minimization done")
        c2 = np.reshape(diff_ar, (len(A),len(K)))

        fits_variance_name = file_name_no_dir.replace("flt.fits", "dark_%s_diff.fits"%(i+1))
        fits.writeto("%sUV_DARK_SUB/FITS/%s"%(primary_dir, fits_variance_name), data = c2, header = header, overwrite= True )

        minimum_var[i] = np.min(diff_ar)
        A_min[i] = par[0]#A[c[0]]
        K_min[i] = par[1]#K[c[1]]
        print ("dark %s  done"%(i+1))

   
    ind = np.argmin(minimum_var)
    print ("dark minimum index ", ind)
    table_dark.add_row(( file_name_no_dir, m1, dark_radii, ind, A_min[ind], K_min[ind], np.min(diff_ar), exp_time ))
    table_name_dark = '%sUV_DARK_SUB/TXT/FLT_drk_sub_gal_%s.txt'%(primary_dir, gal_num+1)
    table_dark.write(table_name_dark, format='ascii', overwrite = True)

    print ("plotting the minimization values now\n")
    fig, ax= plt.subplots(1,1 , figsize =(8, 8))
    ax.plot(temp, minimum_var, "o", markersize =8 )
    ax.set_xlabel(r" Dark Temp[$^o$ C]")    
    ax.set_ylabel("Minimizer |G[r]| [counts]")
    ax.axvline(x = temp_gal, color ="g", label = "galaxy temperature")
    ax.axvline(x = temp[ind], color = "r", label = "minimum value A= %.1e\n, K = %.1e,\n\
     index = %.1e\n"%(A_min[ind], K_min[ind], ind))
    ax.legend(loc='lower left')
    


    minimizer_png = file_name_no_dir.replace("flt.fits", "minimizer.png")
    
    fig.savefig("%sUV_DARK_SUB/PNG/%s"%(primary_dir, minimizer_png), dvi = 400, bbox_inches = 'tight')

    print ("plotting the difference image now ..... with galaxy\n")

    hdulist = fits.open(os.path.dirname(file_name_flt)+'/UV_FLT/FITS/'+os.path.basename(file_name_flt))
    exp_time = hdulist[0].header["EXPTIME"]
    data_gal = hdulist[1].data 

    hdu_dark_selected = fits.open(dark_RAW[ind])
        ## scaling for difference between dark and galaxy exposure time difference
    data_dark = (A_min[ind]*hdu_dark_selected[1].data*exp_time/float(params["exp_dark"]) + K_min[ind])/exp_time### exp time for darks is 1000 secs
      #data_dark[DQ!=0] = 0.0
    data_gal[DQ!=0] = float(params['dq_bad'])
    data_dark[DQ!=0] = float(params['dq_bad'])

    data_sub = data_gal - data_dark*exp_time
    data_sub[DQ!=0] = float(params['dq_bad'])

    hdulist[1].data = data_sub

    hdu_dark_selected[1].data = data_dark

    ################## main outputs#####################
    dark_selected = file_name_flt.replace("flt.fits", "dark_%s.fits"%(ind+1)) 
    dark_selected = os.path.dirname(dark_selected) +'/UV_DARK_SUB/'+ os.path.basename(dark_selected)
    hdu_dark_selected.writeto(dark_selected, overwrite = True)
    hdu_dark_selected.close()
    

    sub_name = file_name_flt.replace("flt.fits", "drk_flt.fits" )
    sub_name = os.path.dirname(sub_name) +'/UV_DARK_SUB/'+ os.path.basename(sub_name)
    hdulist.writeto (sub_name, overwrite = True, output_verify="ignore")
    hdulist.close()
    
    aper_before =  [np.mean(data_gal[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_dark =  [np.mean(data_dark[masks_annulus[k]]*exp_time) for k in range(len(rad1)-1) ]  
    aper_subtracted =  [np.mean(data_sub[masks_annulus[k]]) for k in range(len(rad1)-1) ]  


    fig, (ax1) = plt.subplots(1,1 , figsize =(8, 8))
    ax1.plot (rad_annulus, aper_before, color = "orange", label = "before")
    ax1.plot (rad_annulus, aper_dark, color = "green", label = "dark")
    ax1.plot (rad_annulus, aper_subtracted, color = "blue", label = "subtracted")
    ax1.set_xlabel("pixels")
    ax1.set_ylabel("annuli mean counts")
    ax1.set_title(" Removed Mask for ULIRG %s exposure %s" %(gal_num+1, file_name_no_dir), fontsize = 16)
    ax1.axhline(y=0, color = 'k')
    y1 = 8.11e-6*exp_time
    ax1.axhline(y= y1, linestyle = '--', color = "k", label = "constant dark ISR y = %.1e"%(y1))
    ax1.legend()
   
    verification_png = file_name_no_dir.replace("flt.fits", "verification")

    fig.savefig("%s/UV_DARK_SUB/PNG/%s"%(primary_dir, verification_png), dvi = 400, bbox_inches = 'tight')

    plt.close(fig)

def main(config_file):
    print ( 'the config file used for this analysis ===>', config_file)
    dark_RAW, dark_FLT,  temp = dark_from_ISR(config_file, 'basic')           

    for i in range (5):
        section_gal = 'NULIRG%s' %(int(i+1))
        print (config_file)
        sky_dark_sub(i,config_file, 'basic', section_gal, 'True', dark_RAW, dark_FLT, temp, dark_perform= True)
if __name__ == '__main__': 
    main(config_file)
########<<<<<<<<<< RAMON COPY RAW FILES >>>>>>>>>>>>>>>>>>>>>
### for copying files from RAMON


'''

import ULIRG_params as param
from astropy.table import Table

tab = Table.read("/home/sourabh/ULIRG_v2/scripts/cont_flt_v6_msk_100.txt",format='ascii')
a = list(tab["file_name"])
from subprocess import call
gal_id_all = param.gal_id_all
for i in range (len(tab)):
    for j in range(len(gal_id_all)):
        if tab["galaxy"][i] == gal_id_all[j] :
            cmd = "scp ramon4:/data/highzgal/sourabh/ULIRG/%s /home/sourabh/ULIRG_v2/%s/"%(a[i],gal_id_all[j])
            call(cmd.split(" "))

'''