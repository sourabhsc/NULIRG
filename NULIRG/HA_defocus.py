"""
Finds the HST focus value (in microns) for a given fits file (or list of
files). Uses a reference table.


"""
#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import argparse
from scipy import interpolate
import numpy as np
from astropy.io import fits as fio
from astropy.io import fits
from utils import basic_params
from utils import FLC_centers
from utils import UV_centers

#-----------------------------------------------------------------------------
# Main
#-----------------------------------------------------------------------------
filt = ['775', '782']
from collections import OrderedDict
import pandas as pd
import os


def findfocus(imlst, refdata, noupd, hdrkey):
    """
    This function takes a list of fits frames, finds the appropriate defocus
    value from the supplied reference data file and writes the info to the
    header.
    """
    # Set instrument offsets
    offsets = {'UVIS': -0.24, 'SBC': -0.25, 'WFC': 0.24}

    # Read data from files
    with open(imlst) as ff:
        images = ff.read().split()
    defoctab = np.loadtxt(refdata)

    # Build interpolating function, NB an error will be raised if a date
    # outside of the time range given in the reference file  is requested.
    finterp = interpolate.interp1d(defoctab[:, 0], defoctab[:, 1],
                                   bounds_error=True)

    defsum = 0.
    first = True
    for image in images:
        hdul = fio.open(image, mode='update')
        tjd = (hdul[0].header['EXPSTART'] + hdul[0].header['EXPEND']) / 2.
        imoffs = offsets[hdul[0].header['DETECTOR']]
        defoc = imoffs + finterp(tjd)
        defsum = defsum + defoc
        if noupd:
            if first:
                print('Image                   Defocus')
            first = False
            print(image + '   ' + str(defoc))
        else:
            hdul[0].header[hdrkey] = (defoc, 'model defocus value in microns')
            hdul[0].header['HISTORY'] = 'Defocus information added to header. ' +\
                'Based on model annual summary file: ' + refdata + ', and using ' +\
                'observation time given in the header (average of STARTEXP ' +\
                'and ENDEXP).'
        hdul.close()
    return defsum / len(images)


#-----------------------------------------------------------------------------
# Script I/O
#-----------------------------------------------------------------------------


def main(config_file):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     fromfile_prefix_chars='@')
    # parser.add_argument('imlst', type=str)
    # parser.add_argument('refdata', type=str)
    parser.add_argument('--noupd', action='store_true', help='do not update ' +
                        'header(s), print to std output instead')
    parser.add_argument('--hdrkey', type=str, default='DFOC_MOD',
                        help='header keyword')
    args = parser.parse_args()
    keys = ['galaxy', 'filter', 'exposures', 'defocus']
    dict_defoc = OrderedDict.fromkeys(keys)
    for x in keys:
        dict_defoc[x] = []
    for i in range(5):
        section_gal = 'NULIRG%s' % (int(i + 1))
        params, params_gal = basic_params(config_file, 'basic', section_gal)
        refdata = params['psf_defocus']
        gal_name = params_gal['name']
        primary_dir = params['data_dir'] + gal_name + '/'
        txt_dir = primary_dir + 'HA_PSF/TXT/'

        for j in range(2):

            imlist_file = primary_dir + 'HA_PSF/TXT/' + 'imlist%s_gal%s.txt' % (filt[j], i + 1)
            #loclist = txt_dir + 'loclist%s_gal%s.txt' % (filt[j], i + 1)
            simpos = txt_dir + 'simpos%s_gal%s.txt' % (filt[j], i + 1)

            images = params_gal['f%s' % (filt[j])].split(',')
            ar = []
            simpos_ra = []
            simpos_dec = []
            #loc_x = []
            #loc_y = []
            for x in images:
                name = primary_dir + 'HA_FLC/FITS/jcmc%s_flc.fits' % (x)
                ar.append(name)
                cent_x, cent_y, cent_ra, cent_dec = UV_centers(params, params_gal, i)
                # cent_x is UV center
                simpos_ra.append(cent_ra)
                simpos_dec.append(cent_dec)
                #pixx, pixy = FLC_centers(cent_ra, cent_dec, ar[1])
                # loc_x.append(pixx)
                # loc_y.append(pixy)

            with open(simpos, 'w') as f:
                [f.write("%s \t %s \n" % (simpos_ra[m], simpos_dec[m])) for m in range(len(ar))]
            with open(imlist_file, 'w') as f:
                [f.write("%s\n" % item) for item in ar]
            meandef = findfocus(imlist_file, refdata, args.noupd, args.hdrkey)

            print('Mean defocus of frames filter %s galaxy %s===> %s' % (filt[j], i + 1, str(meandef)))
            dict_defoc['galaxy'].append('gal_%s' % (i + 1))
            dict_defoc['filter'].append('f' + filt[j])
            dict_defoc['exposures'].append(images[0] + '_' + images[1])
            dict_defoc['defocus'].append(meandef)
    df = pd.DataFrame(dict_defoc)
    df.to_csv(os.path.dirname(config_file) + 'defocus_ULIRG.txt', index=False)


if __name__ == '__main__':
    main(config_file)
