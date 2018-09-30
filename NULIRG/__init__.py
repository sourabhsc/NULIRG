import argparse
import os
import sys
import glob
from astropy.io import fits
import numpy as np
from collections import OrderedDict
import pandas as pd
from operator import itemgetter, attrgetter
import shutil

import configparser
from configparser import ExtendedInterpolation

# sys.path.insert(0, '/home/sourabh/ULIRG_v2/scripts/')
# from util import *
from utils import basic_params
from utils import mkdirp
from astropy.io import ascii
#    print("Call your main application code here")

__version__ = '0.0.3'


class make_directories():
    ''' A class to make sub directories given the directory where data is downloaded

    Attributes:
                    data_dir : directory where ULIRG data is downloaded
    '''

    def __init__(self, data_dir, subdir, subsubdir, sub_filetypes):
        ''' function that creatd the sub directories'''

        self.data_dir = data_dir
        self.subdir = subdir
        self.subsubdir = subsubdir
        self.sub_filetypes = sub_filetypes

        [mkdirp(str(args.data_dir + galaxy)) for galaxy in subdir]
        [mkdirp(str(args.data_dir + galaxy + '/' + subsub)) for galaxy in subdir
         for subsub in subsubdir]
        [mkdirp(str(args.data_dir + galaxy + '/' + subsub + '/' + type))
         for galaxy in subdir for subsub in subsubdir for type in sub_filetypes]

    def list_fits(self, text_file, flat_list, dict_entries, keys):
        ''' function to create txt file with all FLT files'''

        dict_files = OrderedDict.fromkeys(dict_entries)
        hdu_list = [fits.getheader(file) + fits.getheader(file, 'sci') for file in flat_list]
        i = 0
        for entry in dict_entries:

            dict_files[entry] = [hdu[keys[i]] for hdu in hdu_list]
            i = i + 1
        dict_files['filtername'] = [dict_files['filtername'][i].
                                    replace('CLEAR1L', 'FR782N') for i in range(len(hdu_list))]

        df = pd.DataFrame(dict_files)
        df.to_csv(text_file, index=False)
        d1 = ascii.read(text_file)
        ascii.write(d1, 'test.txt', format='rst')
        return dict_files

    def copy_files(self, dict_files, galaxy_name, keep_original):
        for k in range(len(dict_files['filename'])):
            if dict_files['detector'][k] == 'SBC' and dict_files['galaxy'][k] == galaxy_name:
                if keep_original:

                    shutil.copy(self.data_dir + '/' + dict_files['filename'][k], self.data_dir + '%s/%s/' % (galaxy_name, 'UV_FLT/FITS/'))
                else:

                    shutil.move(self.data_dir + '/' + dict_files['filename'][k], self.data_dir + '%s/%s/' % (galaxy_name, 'UV_FLT/FITS/'))
            if dict_files['detector'][k] == 'WFC' and dict_files['galaxy'][k] == galaxy_name:
                if keep_original:

                    shutil.copy(self.data_dir + '/' + dict_files['filename'][k], self.data_dir + '%s/%s/' % (galaxy_name, 'UV_FLT/FITS/'))
                else:
                    shutil.move(self.data_dir + '/' + dict_files['filename'][k], self.data_dir + '%s/%s/' % (galaxy_name, 'HA_FLC/FITS/'))

    # def dark_fits(self, )


def dark_selection(mean_obsdate, dark_dir, data_dir, dark_type):
    dark_FLT = sorted(glob.glob(dark_dir + '/*%s.fits' % (dark_type)))
    obs_date_dark = [fits.getheader(x)['DATE-OBS'] for x in dark_FLT]
    date = np.unique(obs_date_dark)
    dark_date = date.astype('datetime64')
    galaxies_date = mean_obsdate.astype('datetime64')
    selected_darks = np.argmin(abs(dark_date - galaxies_date))
    dark_chosen = []
    for i in range(len(dark_FLT)):
        if obs_date_dark[i] == date[selected_darks]:
            dark_chosen.append(dark_FLT[i])
    [shutil.copy(x, data_dir + 'dark/' + os.path.basename(x)) for x in dark_chosen]
    return data_dir + 'dark/'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creating directory structure for ULIRGs ')
    package_dir = os.getcwd().replace('NULIRG', '')
    script_dir = os.getcwd() + '/'
    parser.add_argument('--data_dir', default=os.getcwd().replace('NULIRG', 'data/'),
                        help='Directory where all FLT files have been saved ')

    config_dir = package_dir + 'config/'
    parser.add_argument('--config_file', default=config_dir + 'params_default.cfg',
                        help='config file with ULIRG parameters specified')

    dark_dir = os.getcwd().replace('NULIRG', 'data_dark/')
    parser.add_argument('--dark_dir', default=dark_dir,
                        help='Directory where all dark files are saved')
    args = parser.parse_args()
    print('***data_directory ==> ', args.data_dir)
    print('***config file used for analysis ==> ', args.config_file)

    print ('1) creating sub directories\n')
    file_types = ['flt', 'flc']
    flt_files = sorted([glob.glob(args.data_dir + '*%s*' % (type)) for type in file_types])
    flat_list = [item for sublist in flt_files for item in sublist]
    flat_list = sorted(flat_list)
    galaxies = []
    for f in flat_list:
        hdr = fits.getheader(f)
        galaxies.append(hdr['TARGNAME'])
    subdir = sorted(np.unique(galaxies))
    print ('\n...number of galaxies ', len(np.unique(galaxies)))
    print ('\n...', sorted(np.unique(galaxies)))

    subsubdir = ['UV_RAW', 'UV_FLT', 'UV_DARK_SUB', 'UV_ALIGN', 'UV_DRZ', 'UV_PSF',
                 'HA_FLC', 'HA_INTER', 'HA_DRC', 'HA_PSF']
    sub_filetypes = ['FITS', 'PNG', 'TXT']
    dict_entries = ['filename', 'instrument', 'detector', 'galaxy', 'filtername', 'obsdate', 'orientation']
    keys = ['FILENAME', 'INSTRUME', 'DETECTOR', 'TARGNAME', 'FILTER1', 'DATE-OBS', 'ORIENTAT']

    data_dir = args.data_dir
    text_file = config_dir + 'FLT_list.txt'

    dir_instance = make_directories(data_dir, subdir, subsubdir, sub_filetypes)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ('2) creating text file with names of galaxies and filters \n ')
    print ('Text file with flt files data ===>', text_file, )

    dict_files = dir_instance.list_fits(text_file, flat_list, dict_entries, keys)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ('3) selecting the dark subtraction files\n')
    print ('Directory with all the dark observations ===>', args.dark_dir)
    mean_obsdate = (np.array(dict_files['obsdate'], dtype='datetime64[s]')
                    .view('i8')
                    .mean()
                    .astype('datetime64[s]'))
    print ('\n mean observation date of galaxies ===>', mean_obsdate)

    dark_dir_chosen = dark_selection(mean_obsdate, dark_dir, data_dir, dark_type='flt')
    dark_dir_chosen = dark_selection(mean_obsdate, dark_dir, data_dir, dark_type='raw')

    print ('Directory of selected 20 dark exposures\n closest\
     to obs date of galaxies ===>', dark_dir_chosen)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("4) Populating the config file \n")
    print ("<<<<< USING copy of 'ULIRG_params.cfg' params_default.cfg>>>>>> \n")
    print ('config file is ===> ', config_dir + args.config_file)
    default_config = args.config_file

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("5) Moving input FLT and FLC files\n")
    print (np.unique(galaxies))
    [dir_instance.copy_files(dict_files, y, keep_original=True) for y in np.unique(galaxies)]
    print ("\n6) dark subtraction \n")

    # import UV_dark_main
    # UV_dark_main.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("\n 6) xregister aligning and rotating UV images \n")
    # import UV_xreg_rotate
    # UV_xreg_rotate.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("\n 7) combining all extensions from output of aligning images \n")
    # import UV_header_rename
    # UV_header_rename.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("\n 8) drizzling the aligned images using astrodrizzle \n")
    # import UV_astrodrizzle
    # UV_astrodrizzle.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("\n 9) PSFmatching for UV images\n")
    # import UV_psfmatch
    # UV_psfmatch.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

##########<<<<<<<<<<<<<<<<< HA codes now >>>>>>>>>>>>>>>>>>>>>>> #########################
    print ("\n 10) HA image align and cut\n")
    #import HA_align
    # HA_align.main(default_config)
    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')

    print ("\n 11) HA image find defocus and PSFs \n")
    #import HA_defocus
    # HA_defocus.main(default_config)
    #import HA_mkpsf
    # HA_mkpsf.main(default_config)
    import HA_psf_cut
    HA_psf_cut.main(default_config)
    import HA_psfmatch
    HA_psfmatch.main(default_config)

    print ('\n----------------------------------------------------------\n')
    print ('----------------------------------------------------------\n')


#def setup(app):
#    app.add_html_theme('dask_sphinx_theme', path.abspath(path.dirname(__file__)))

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
