# Data reduction steps for lyman alpha imaging of in NULIRGs
This pipeline can be used to identify a line emission image from 4 low pass HST SBC filters.

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




