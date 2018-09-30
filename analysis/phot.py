import photutils
from photutils import create_matching_kernel
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import scipy
from photutils import TopHatWindow
from photutils import CosineBellWindow

window = CosineBellWindow(alpha=0.35)
from scipy.signal import convolve as scipy_convolve
window = TopHatWindow(0.35)

dir1 = '/home/sourabh/ULIRG_package/data/OPTICAL_PSF/'
data_ref = fits.getdata(dir1 + 'f165psf.fits')
data_psf = fits.getdata(dir1 + 'PSF_775_gal4_rotate_cut.fits')
kernel = create_matching_kernel(data_psf, data_ref)  # , window = window )
fits.writeto('ker.fits', data=kernel, overwrite=True)
plt.imshow(kernel, cmap='Greys_r', origin='lower')
filename = '/home/sourabh/ULIRG_package/data/IRASF10594+3818/gal1_HA.fits'
fileout = '/home/sourabh/ULIRG_package/data/IRASF10594+3818/gal1_HA_psfmatch.fits'
ker = 'ker.fits'
#ker_shift = np.pad(kernel, ((0, 1), (0, 1)), mode='constant')
data1 = scipy_convolve(data_psf, kernel, mode='same')
fits.writeto('test2.fits', data=data1, overwrite=True)
data3 = data1 - data_ref
fits.writeto('test3.fits', data=data3, overwrite=True)


def psf_match(filename, fileout, ker):
    hdulist = fits.open(filename)
    data = hdulist[0].data

    hdu_ker = fits.open(ker)
    ker_data = hdu_ker[0].data
    ker_shift = np.pad(ker_data, ((0, 1), (0, 1)), mode='constant')

    data_out = scipy_convolve(data, ker_shift, mode='same', method='fft')  # convolve

    hdulist[0].data = data_out
    fits.writeto('test.fits', data=data_out - data, overwrite=True)
    hdulist.writeto(fileout, overwrite=True, output_verify="ignore")


psf_match(filename, fileout, ker)


plt.colorbar()
plt.show()


# In[36]:


from matplotlib import pylab

params = {'legend.fontsize': 10,
          'figure.figsize': (15, 7),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'xx-large',
          'ytick.labelsize': 'xx-large'}
import pylab as plot
plot.rcParams.update(params)
pylab.rcParams.update(params)


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


cent_x = 65
cent_y = 65
width = 3
aper_lim = 120
nx = 130
ny = 130
rad1, rad_annulus, masks, masks_annulus = masks_circular(cent_x, cent_y, width, aper_lim, nx, ny)


aper_ref = [(np.mean(data_ref[masks_annulus[k]])) for k in range(len(rad1) - 1)]
aper_psf = [(np.mean(data_psf[masks_annulus[k]])) for k in range(len(rad1) - 1)]
#plt.plot(rad_annulus, aper_ref/aper_ref[0], label ='f165')
#plt.plot(rad_annulus, aper_psf/aper_psf[0], label = 'f775')
aper_ref_sum = [(np.sum(data_ref[masks[k]])) for k in range(len(rad1))]
aper_psf_sum = [(np.sum(data_psf[masks[k]])) for k in range(len(rad1))]


plt.plot(rad1, aper_ref_sum, label='f165')
plt.plot(rad1, aper_psf_sum, label='f775')

#plt.plot(rad_annulus, aper_ref/aper_ref[0], label ='f165')
#plt.plot(rad_annulus, aper_psf/aper_psf[0], label = 'f775')

plt.legend()
plt.show()
