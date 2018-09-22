*******************
Configuration codes
*******************




Installation
============



Data Download
-------------
    * Create the following directory structure.::

        subdir = ['data', 'data_dark', 'config', 'NULIRG']

    * Copy all the codes in 'NULIRG' folder. 'config' folder will have all the config files as listed in the next section. Change the 'data_dir; in 'params_default.cfg'.

    * **dark data** :- Download from .... put it in 'data_dark' folder inside your directory
    * **ULIRG data** :- Download from .... put it in 'data' folder inside your directory
    * **Running __init__.py** :- This code with default config file for ULIRGs (params_default.cfg) creates the appropriate directory structure.::


        subsubdir = ['UV_RAW', 'UV_FLT', 'UV_DARK_SUB', 'UV_ALIGN', 'UV_DRZ', 'UV_PSF',
                 'HA_FLC', 'HA_INTER', 'HA_DRC', 'HA_PSF']

        subsubsubdir = ['FITS','TXT', 'PNG' ]

Config files
------------

These files are located in `/DIR/config/` .
    
    * astrodrizzle_SBC_conf.cfg 
    * astrodrizzle_WFC_conf.cfg
    * defocus.dat
    * params_default.cfg
    * sex_default.conf
    * tinybase775_gal1.txt
    * tinybase775_gal2.txt
    * tinybase775_gal3.txt
    * tinybase775_gal4.txt
    * tinybase775_gal5.txt
    * tinybase782_gal1.txt
    * tinybase782_gal2.txt
    * tinybase782_gal3.txt
    * tinybase782_gal4.txt
    * tinybase782_gal5.txt

Requirements
============
    * Pyraf
    * astropy
    * tinytim
    * configparser
    * argparser
    * numpy 
    * astropy
    * drizzlepac
    * scipy 
    * zap




utilities function
==================


.. automodule:: utils
   :members: 
