# Data reduction steps for lyman alpha imaging of in NULIRGs


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Requirements
Here is a list of softwares that needs to be installed before running the code:- 

```
conda install configparser
sudo apt-get install sextractor
conda install astropy
conda install matplotlib
conda install scipy
```


### Structure and set up

All the input parameters can be specified in the 'configuration file':- 'dark_sub.cfg'
In the coding directory one needs to have the following set of files :-

 * main code --> 'dark_sub_SBC.py'
 * configuration file --> 'dark_sub.cfg'
 * sextractor files
 ```
 default.conv
 default.param
 sex_default.config
```

* create a 'work_dir' where all the input files are present. 
* within 'work_dir', subdirectories will be created while running the code to store intermediate files
* all the output files will also be created in 'work_dir'
* 'dark_dir':- a directory where all the dark exposres are stored which are closest to the observation date of your observation.

## Inputs in 'dark_sub.cfg' files :-
* name :- In the current form, the configuration file is created for a given galaxy. For separate galaxies, one needs to create separate configuration files
* input_frames :- list out all the input exoposures separated by comma and no spaces. Name of input files should be in the standard form of SBC images e.g. 'jcmc11e6q_flt.fits'
* temp_threshold :- Temperature cutoff to decide when one needs dark subtraction and when only sky subtraction is enough for removing background. For our case, we are taking T>22C for hot frames that need dark subtraction
* chisq :- True or False , if chi square minimization method is selected for selecting accurate dark exposure.
    
     old method:-- method =='|G(r)|'
     chi sq methd :- method =='chisq'  
     
* mask_radii :- i) 0 , for the case when only sextractor segmentation maps are chosen for masking the galaxy
        ii) R_95 (in pixels),  non-zero radius selected to cover 95% flux of the full image. for the example case, I took R =250. Here galaxy is masked by large circular region instead of sextractor maps.
* sky_min, sky_max :- For cold frames minimimum and maximum radius of circular aperture for sky. Mean of annulus mean is              taken as sky value where annulus go from 'sky_min' to 'sky_max'.
* dq, dq_bad:- (1000, 0.0) for differentiating between additional hot pixels found on SBC detector and pixels with [DQ!=0]

* sex_config :- sex_default.conf

* smooth_scale :- 3 , gaussian smoothing scale which is used to smooth the galaxy image and segementation maps on the   '''         smoothed image are found to mask the galaxy.
```
* Sextractor parameters
 i) detect_thresh :- 8.0
 ii) detect_minarea :- 3.5
 iii) analysis_thresh :- 3.5
```

Some additional parameters are listed below:-

```
###  parameters for annulus
cent_x = 512
cent_y = 512
width = 1
aper_lim = 500 
###grid search parameters ###
a_lim = 5
k_lim = 0.02
del_a = 0.1
del_k = 0.001
exp_dark = 1000
nx = 1024 
ny = 1024 
### additional hot pixels found on SBC detector #####
hot_pix_x=698,669
hot_pix_y=420,410
```


## Output files :--
Each run of the code produces output dark subtracted files for both minimization method 'G(r)' and 'chisq'. Sextractor maps can be used if one 'mask_radii =0' in configuration files. 

   Final output is the dark subtracted frame with name 
   ```
   dark_sub_name = input_flt.replace("flt", "drk_flt_%s_dr_%s" %(method, dark_radii)) 
   ```
   method = 'Gr' or 'chisq', based on the method of minization
   dark_radii = '0' or 'R95', based on sextraction segmentaion maps or not. 
   
   e.g. For 1)  jctd43pqq_drk_flt_chisq_dr_250.fits, method = chisq, dr =250(no sextractor maps) 
            2)  jctd43pqq_drk_flt_Gr_dr_0.fits, method = G(r), dr = 0 (sextractor maps are used)
   
   The dark frame for the given dark subtraction is saved as 
   
   ```
   corrected_dark_frame = input_flt.replace("flt.fits", "dark_%s_%s_dr_%s.fits"%(index of selected dark, method, dark_radii )
   jctd43pqq_dark_18_chisq_dr_0.fits
   ```
   
   * Sky subtracted output :- 
   ```
    sky_sub_name = input_flt.replace("flt", "sky_flt")
   ```
   
I am still updaing the documentaion for intermediate files.

b) Intemediate files :---
 i) FITS:-- 

    fits files with the minimizer values will be saved as 

    ----->>>> input_hot_suffix + "_dark_" + "i_" + method +'.fits'
    
    Subtracted frame for each dark index is saved as :-

    ----->>>>> input_hot_suffix + "_dark_" + "i_" + 'diff_' + method +'.fits' 

    where i goes from 0-20 corresponding to the 20 dark frames while  method takes two values depending upon minimization process:--
 ii) PNG:--- 
    Intermediatiate PNG files for photometry step, verification step and value of minimizer as a function of dark are also saved.
 iii) TXT:---
   The text file stores the information of parameters , A, K and ind (refer to documentation)


Briefly these are all the steps  that are performed:--
1) creating directories for intermediate files ... 
2) creating apertures for circular annuli
3) performing photometry on input cold and hot frames
4) creating table to write the output for dark subtraction process
5) Performing dark subtraction
6) plotting the minimization values now
7) plotting the difference image now ..... with galaxy



## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/sourabhsc/NULIRG/tags). 

## Authors

* **Sourabh Chauhan** - *Initial work* - [sourabhsc](https://github.com/sourabhsc)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.




