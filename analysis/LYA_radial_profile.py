
# coding: utf-8

# ### Here is the list:-
# * all ULIRGs check size and run all ...
#     * test for S/N cut
#     * get error on lya fluxes
#     * goodness of fit
#     * radial profiles and comparison
#     * host the interface
# * start new draft 
#     * add dark subtraction
#     * add sdss analysis, ppxf fitting and make it work
#     * LARS, GOAL , CM15 data
#     
# * additional research
# 
#     * add SB99 or cloudy analysis to match a few numbers
#     * Look for clues IR herschel and molecular lines methods
#     * Look for cosmic ray feedback methods etc
#     * 
# * additional 
#     * improve web interface 
#     * update website ULIRGs
#     
# * NN interface learning 
#     * finish that COURSERA course
#     * finish paper on CNN and mahcine learning 
#     * 
#     
#     
# 
# 

# In[1]:


import sys
sys.path.append("/home/sourabh/ULIRG_package/NULIRG/") # go to parent dir
config_file = "/home/sourabh/ULIRG_package/config/ULIRG_params.cfg"
from astropy.io import fits
from matplotlib import pyplot as plt
import utils
from utils import UV_centers_drz
from utils import masks_circular
from utils import basic_params
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from astropy.visualization import (MinMaxInterval, AsinhStretch,
								   ImageNormalize)
from scipy.ndimage.filters import gaussian_filter

import glob
import numpy as np
from matplotlib import pylab
plt.rcdefaults()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 7),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)
from utils import circle
import matplotlib
from collections import OrderedDict
matplotlib.__version__
import pandas as pd


# In[2]:


min_contour = ([0.05,0.02, 0.01, 0.01, 0.02])

def plot_rectangle(x,y, szx, szy, color, ax):
    utils.rectangle(x, y, szx, szy, color, ax)
def sum_ext(i, ext ):
    files_all = glob.glob('./ULIRG%s_test/*%s*'%(i+1, ext ))
    s = [fits.getdata(f) for f in files_all]
    sum_all = np.zeros((data.shape[0], data.shape[1]))
    for s1 in s:
        sum_all = sum_all + s1
    return sum_all
scale = [2.750, 2.750,2.750,2.750, 2.349 ]
#sz = np.array([80, 120, 180, 180, 150])
sz = np.array([500, 500, 500, 500, 500])

petro_all = np.array([8.95, 18.9, 7.37, 16.3, 9.53])
dl = np.array([754.9, 754.9, 749.7, 760.2, 615.3] )# lum distance in mpc using 
dl_all = dl*1e6*3.08568e18

color =['r', 'g', 'b', 'c', 'm', 'g', 'y', 'orange', 'r', 'g', 'b', 'c', 'm', 'g', 'y', 'orange']

szx = 200
szy = 200


def non_zero(data):
    c1 = data[data!=0]
    return len (c1)



def 
for i in range(5):
    file125 = 'gal%s_UV_F125_scale_04_psfmatch.fits'%(i+1)
    hd = fits.open(file125)
    photflam = float(hd[1].header['PHOTFLAM'])
    data_125 = hd[1].data
    section_gal = 'NULIRG%s' % (int(i + 1))
    params, params_gal = basic_params(config_file, 'basic', section_gal)
    cent_x, cent_y = UV_centers_drz(file125)


    print ('=============================================\n')
    print ('<<<<<<<<<<<<<<<<<<<<<ULIRG %s>>>>>>>>>>>>>>>>>>\n'%(i+1))
    print ('=============================================\n')

    print ('Original shape', data.shape)
    print ('UV 125 centers --> (%s, %s)'% (cent_x, cent_y))

    

    
    print (' ----> Making annuli masks\n')
    width = 4.0
    aper_lim = 300

    if i ==4:
        aper_lim = 300*0.158/0.131
    nx = data_125.shape[0]
    ny = data_125.shape[1]
    rad1, rad_annulus, masks, masks_annulus = masks_circular(cent_x, cent_y, width, aper_lim, nx, ny)

    ################## computing center boxes ################
    if i ==1:
        x = [150, 150, 150, 350, 350, 350, 550, 550, 550]
        y = [150, 350, 550, 150, 350, 550, 150, 350, 550]
    elif i ==2:
        x = [150, 150, 150, 150, 350, 350, 350, 350, 550, 550, 550, 550, 750, 750, 750, 750]
        y = [150, 350, 550, 750, 150, 350, 550, 750, 150, 350, 550, 750, 150, 350, 550, 750]
    else:
        x = [350, 350, 350, 550, 550, 550, 750, 750, 750]
        y = [350, 550, 750, 350, 550, 750, 350, 550, 750]
      
        
    fig,((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize = (30,22))
    
    [plot_rectangle(x[i],y[i], szx, szy, color[i], ax1) for i in range(len(x))]
    norm = ImageNormalize(data_125, stretch=AsinhStretch())
    img = ax1.imshow(data_125, origin = 'lower', vmin =0.0005, vmax = 0.001, norm = norm, cmap = plt.cm.Greys_r, zorder =0)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.0)
    cbar = plt.colorbar(img,cax=cax, orientation='vertical')
    cbar.minorticks_on
    

    data_lya = sum_ext(i, 'A_lmfit_x')
    err_lya = sum_ext(i, 'A_lmfit_err')
    ### positive flux cond
    
    #data_all[data_all<0] = 0.0
    cond_positive = data_lya < 0.0
    cond_negative = data_lya > 0.0

    #### SN >1 condition
    
    SN_all = abs(data_lya/err_lya)
    

    #SN_add = data_lya + err_lya
    #SN_diff = data_lya - err_lya

    #cond_sn = np.isnan(SN_all)
    #SN_all[cond_sn] = 0.0
    #SN_all[cond_sn_tot] = 0.0
    b_all = sum_ext(i, 'b_lmfit_x')
    cond_b = abs(abs(b_all)-3) <1e-5

    status_all = sum_ext(i, 'status_lmfit_x')
    cond_status = status_all == 0
    
    data_neg = np.copy(data_all)
    data_test = np.copy(data_all)
    aper_total =  [len(data_neg[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    
    aper_nonzero_total =  [non_zero(data_test[masks_annulus[k]]) for k in range(len(rad1)-1) ] 
    cond_test = cond_sn_tot + cond_b + cond_status
    data_test[cond_test] = 0.0
    aper_nonzero_total_cutoff = [non_zero(data_test[masks_annulus[k]]) for k in range(len(rad1)-1) ] 

    
    cond_all = cond_sn_tot + cond_b + cond_status + cond_positive   
    data_all[cond_all] = 0.0
    SN_add[cond_all] = 0.0
    SN_diff[cond_all] = 0.0
    
    cond_neg = cond_sn_tot + cond_b + cond_status + cond_negative   

    data_neg[cond_neg] = 0.0

    #data_neg = abs(data_neg)
    #img = ax1.imshow(data_all, origin = 'lower', vmin =-1e-15, vmax = 1e-20, norm = norm, cmap = plt.cm.Greys_r, zorder =0)

    
    fits.writeto('./ULIRG%s_test/data%s.fits'%(i+1, i+1), data = data_test, overwrite = True)
    fits.writeto('./ULIRG%s_test/SN%s.fits'%(i+1, i+1), data = SN_all, overwrite = True)
    fits.writeto('./ULIRG%s_test/b%s.fits'%(i+1, i+1), data = b_all, overwrite = True)
    fits.writeto('./ULIRG%s_test/status%s.fits'%(i+1, i+1), data = status_all, overwrite = True)

    

    aper_SN_low =  [np.mean(SN_add[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    aper_SN_high =  [np.mean(SN_diff[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    
    
    sum_SN_low =  [np.sum(SN_add[masks[k]]) for k in range(len(rad1)) ]  
    sum_SN_high =  [np.sum(SN_diff[masks[k]]) for k in range(len(rad1)) ]  

    
    sum_125 =  [np.sum(data_all[masks[k]]) for k in range(len(rad1)) ]  
    aper_125 =  [np.mean(data_all[masks_annulus[k]]) for k in range(len(rad1)-1) ]  
    
    sum_125_neg =  [np.sum(abs(data_neg[masks[k]])) for k in range(len(rad1)) ]  
    aper_125_neg =  [np.mean(abs(data_neg[masks_annulus[k]])) for k in range(len(rad1)-1) ]  
    
    sum_UV =  [np.nansum(data[masks[k]]*photflam*150) for k in range(len(rad1)) ]  
    aper_UV =  [np.nanmean(data[masks_annulus[k]]*photflam*150) for k in range(len(rad1)-1) ] 
    
    aper_nonzero =  [non_zero(data_all[masks_annulus[k]]) for k in range(len(rad1)-1) ] 
    aper_nonzero_neg =  [non_zero(data_neg[masks_annulus[k]]) for k in range(len(rad1)-1) ]  

    rad_kpc = [rad*0.04*scale[i] for rad in rad1]
    rad_annulus_kpc = [rad*0.04*scale[i] for rad in rad_annulus]
    diff = [(-aper_125_neg[k] + aper_125[k]) for k in range(len (rad_annulus))]
    
    diff_sum = [(-sum_125_neg[k] + sum_125[k]) for k in range(len (rad1))]

    #norm_sum = [x1/sum_125[-1] for x1 in sum_125]

    #low = [aper_125[i] -aper_125[i] /aper_SN[i] for i in range(len (rad_annulus))]
    #high = [aper_125[i] +aper_125[i]/aper_SN[i] for i in range(len (rad_annulus))]

    ax4.plot(rad_kpc, sum_125, 'r', label = r'$Ly\alpha$ positive')
    ax4.fill_between(rad_kpc, sum_SN_low, sum_SN_high, alpha =0.1, color = 'b', zorder =1)
    
    ax4.plot(rad_kpc, sum_125_neg, 'g', label = r'abs($Ly\alpha$ negative)')
    ax4.plot(rad_kpc, diff_sum, 'b', label = r'($Ly\alpha$ total)')
    ax4.plot(rad_kpc, sum_UV, 'k', label = 'F125 flux')

    
    
    ax3.plot(rad_annulus_kpc, aper_125, 'r', label = r'$Ly\alpha$ positive')

    ax3.fill_between(rad_annulus_kpc, aper_SN_low, aper_SN_high, alpha =0.1, color = 'b', zorder =1)

    ax3.plot(rad_annulus_kpc, aper_125_neg, 'g', label = r'abs($Ly\alpha$ negative)')
    ax3.plot(rad_annulus_kpc, diff, 'b', label = r'($Ly\alpha$ total)')
    ax3.plot(rad_annulus_kpc, aper_UV, 'k', lw = 2.0, label = 'F125 flux')
    
    
    theoretical = [2*np.pi*k*width for k in rad_annulus]
    ax2.plot(rad_annulus_kpc, theoretical, color = 'k', zorder =2,lw = 2.0, label = 'Expected')
    
    ax2.fill_between(rad_annulus_kpc, aper_total, alpha =0.1, color = 'b', zorder =1, label = 'Total')
    ax2.fill_between(rad_annulus_kpc, aper_nonzero_total, alpha =0.2, color = 'c', zorder =1, label = 'Total non-zero')
    ax2.fill_between(rad_annulus_kpc, aper_nonzero_total_cutoff, alpha =0.3, color = 'm', zorder =1, label = 'Quality fit cutoff')

    ax2.fill_between(rad_annulus_kpc, aper_nonzero, alpha =0.3, color = 'r', zorder =2, label ='Positive')
    ax2.fill_between(rad_annulus_kpc, aper_nonzero_neg, alpha =0.6, color = 'g', zorder =3, label = 'Negative')


    ''' #testing the total counts
    test = [aper_nonzero[k] + aper_nonzero_neg[k]- aper_nonzero_total_cutoff[k] for k in range(len (rad_annulus))]
    ax1.fill_between(rad_annulus_kpc, test, alpha =0.5, color = 'g', zorder =2, label = 'negative')
    '''
    #ax2.plot(rad_annulus_kpc, aper_UV, )
    
    
    
    ax3.set_xlim(0, 32)
    ax4.set_xlim(0, 32)
    ax2.set_xlim(0, 32)
    ax2.set_yscale('log')
    ax2.set_ylabel('Number of points')
    ax2.set_xlabel('Distance from UV center [kpc]')
    
    ax3.set_ylabel(r'Aperture mean Flux[in ergs-$s^{-1}$-$cm^{-2}$-arcsec$^{-2}$]')
    ax3.set_xlabel('Distance from UV center [kpc]')
    
    
    ax3.set_yscale('log')
    
    ax3.axvline(x = petro_all[i], color = 'm', linestyle='--', label = r'$2\times R_p$')
    ax4.axvline(x = petro_all[i], color = 'm', linestyle = '--', label = r'$2\times R_p$')
    
    ax2.legend(loc ='lower left')
    ax3.legend(loc ='upper right')
    
    ax1.set_xlim(0,1000)
    ax1.set_ylim(0,1000)
    ax1.plot(cent_x, cent_y, marker ='x', color = 'b', mew=5, ms=20)
    
    circle(cent_x, cent_y, aper_lim, 'c', ax1)
    circle(cent_x, cent_y, petro_all[i]/(0.05*scale[i]), 'm', ax1)
    
    new_label_x = np.linspace((cent_x-sz[i]), (cent_x+sz[i]), 10)#*0.05*phys )
    new_label_y = np.linspace((cent_y-sz[i]), (cent_y+sz[i]), 10)#*0.05*phys )
    ar_x = np.round((new_label_x- cent_x)*0.04*scale[i] )
    ar_y = np.round((new_label_y - cent_y)*0.04*scale[i] )

    ax1.set_yticks(new_label_y)
    ax1.set_yticklabels(int(x) for x in (ar_y))
    ax1.set_xticks(new_label_x)
    ax1.set_xticklabels((int(x) for x in (ar_x)))
    ax1.set_xlabel('kpc')
    ax1.set_ylabel('kpc')
    ax1.set_aspect(1.)
    ax4.legend(loc = 'upper left')
    
    
    ax4.set_ylabel(r'Cumulative Flux[in ergs-$s^{-1}$-$cm^{-2}$-arcsec$^{-2}$]')
    ax4.set_xlabel('Distance from UV center [kpc]')
    
    
    keys = ['radius', 'Lya_pos', 'UV', 'Lya_neg', 'Lya_total', ]
    
    dict_annuli = OrderedDict.fromkeys(keys)
    dict_annuli['radius'] = rad_annulus_kpc
    dict_annuli['Lya_pos'] = aper_125
    dict_annuli['UV'] = aper_UV
    dict_annuli['Lya_neg'] = aper_125_neg
    dict_annuli['Lya_total'] = diff
    
    
    df = pd.DataFrame.from_dict(dict_annuli)
    df.to_csv('ULIRG%s_annuli.csv'%(i+1), columns = keys, index = True)
    
    
    dict_sum = OrderedDict.fromkeys(keys)
    dict_sum['radius'] = rad_kpc
    dict_sum['Lya_pos'] = sum_125
    dict_sum['UV'] = sum_UV
    dict_sum['Lya_neg'] = sum_125_neg
    dict_sum['Lya_total'] = diff_sum
    
    
    df = pd.DataFrame.from_dict(dict_sum)
    df.to_csv('ULIRG%s_sum.csv'%(i+1), columns = keys, index = True)
    
    
    
    ax1.set_title('ULIRG %s F125 image'%(i+1), fontsize = 15)
    fig.savefig('ULIRG%s_radial.png'%(i+1), dvi = 400)
    plt.show()

