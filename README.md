# Data reduction steps for Lyman alpha imaging of in NULIRGs
This pipeline can be used to identify a line emission image from 4 low pass HST SBC filters.

## Getting Started


### Requirements
Here is a list of softwares that needs to be installed before running the code:- 

```
sudo apt-get install sextractor

conda install configparser
conda install astropy
conda install matplotlib
conda install scipy

conda install acstools
conda install pyfiglet
conda install termcolor

```

### Data reduction steps

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. Here are the steps to run the pipeline:-
1. [Required] Download the FLT and FLC files for SBC, WFC filters for different ULIRGs in the working directory **DIR**.
2. [Required] Download dark files for the time of observation in a sub folder **data_dark**.
3. [Required] Download config files and place all of them in sub directory **config**.
4. [Optional] Update **params_default.cfg** text file by changing **package_dir** and **data_dir** keywords. (This step is not required if this repository has been downloaded from here.)

At the end of above step we would have following subdirectories:-
```
subdir = ['data', 'data_dark', 'config', 'NULIRG']
```
Now we are ready to run the code:-

* Go to **package_dir/NULIRG** and run 
```
python __init__.py
```
* Step 1. This will create all sub directories for all the ULIRGs

```
subsubdir = ['UV_RAW', 'UV_FLT', 'UV_DARK_SUB', 'UV_ALIGN', 'UV_DRZ', 'UV_PSF',
         'HA_FLC', 'HA_INTER', 'HA_DRC', 'HA_PSF']

subsubsubdir = ['FITS','TXT', 'PNG' ]

```
* Step 2. Make a text file with list of FLC and FLT files. 
* Step 3. Selecting all the exposures that require dark subtraction
* Step 4. [optional] Populate the parameters in config file. (Adding this feature later on. For now, we are taking default config for this run)
* Step 5. [optional] Make a backup copy of FLT files. (Not using for now)
* Step 6. dark subtraction---> uncomment this
```
    # import UV_dark_main
    # UV_dark_main.main(default_config)
```
* Step 6. **xregister** aligning and rotating UV images --> uncomment this
```
    # import UV_xreg_rotate
    # UV_xreg_rotate.main(default_config)
```
* Step 7. combining all extensions from output of aligning images --> uncomment this
```
    # import UV_header_rename
    # UV_header_rename.main(default_config)
```
* Step 8. drizzling the aligned images using astrodrizzle --> uncomment this
```
    # import UV_astrodrizzle
    # UV_astrodrizzle.main(default_config)
```
* Step 9. PSFmatching for UV images --> uncomment this
```
    # import UV_psfmatch
    # UV_psfmatch.main(default_config)
```
* Step 10. HA image align and cut ---> uncomment this
```
    #import HA_align
    # HA_align.main(default_config)

```
* Step 11. HA image find defocus and PSF matching --> uncomment this
```
    import HA_psf_cut
    HA_psf_cut.main(default_config)
    import HA_psfmatch
    HA_psfmatch.main(default_config)
```
----------------
Now we have all the reduced UV and HA FITS files in directories **UV_PSF/FITS** and **HA_PSF/FITS/** to start fitting for Lyman alpha flux and obtaining segmentation map data. 
----------------



## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/sourabhsc/NULIRG/tags). 

## Authors

* **Sourabh Chauhan** - *Initial work* - [sourabhsc](https://github.com/sourabhsc)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.




