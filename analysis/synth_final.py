
# coding: utf-8

# ## A test for synthetic filter counts as a function of input spectrum parameters
#
# This notebook summarizes a few tests to study the variation in counts as a function of input spectrum of the form
# $$ f_{\lambda} =  {norm}*(m*\lambda^{\beta} + A*\frac{1}{\sqrt{2\pi\sigma}}exp^{-\frac{(\lambda-\lambda_0)^2}{2\sigma^2}})$$
#
# where $$norm = 3 \times 10^{-6}$$
#
#
# counts per seconds are given by
# $$\frac{C(P)}{\Delta t}  = \frac{a}{hc}\int f_{\lambda} P(\lambda) \lambda d\lambda $$
# I am taking :-
# $$ A =1$$
# $$ m =1$$
# $$ \beta = -1 $$
# $$\sigma = 1$$
#
# $$ \frac{a}{hc} = 2.277380717335305\times 10^{20} cm^{-1}{ergs}^{-1} $$
#
# I am trying to calculate the counts in synthetic and original filters

# ### Following code attempts to determine the minimization method to be used for this problem

# ### <<<<<<<<<<<<<<<< IMPORTS >>>>>>>>

# In[ ]:


import numpy as np
import pysynphot as S
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d
import scipy.optimize as op
from matplotlib import pylab
from matplotlib import pyplot as plt
import inspect
from scipy.optimize import curve_fit


params = {'legend.fontsize': 10,
          'figure.figsize': (15, 7),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'xx-large',
          'ytick.labelsize': 'xx-large'}
import pylab as plot
plot.rcParams.update(params)
pylab.rcParams.update(params)


from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import scipy.optimize as op
from collections import OrderedDict


from astropy import constants as c
from astropy import units as u
A = S.refs.PRIMARY_AREA
fact = A / (c.h.cgs * c.c.cgs)

dir1 = '/home/sourabh/ULIRG_package/data/IRASF10594+3818/'

red_shift = np.array([0.158, 0.158, 0.157, 0.159, 0.131])
filters = ['F125', 'F140', 'F150', 'F165']
filter_name = np.array(['f125lp', 'f140lp', 'f150lp', 'f165lp'])
fil_num = np.array([125, 140, 150, 165])

m_true = 1
b_true = -1
A_true = 5
sigma = 1


theta = [m_true, b_true, A_true]


def original_filt(filt):
    bp = S.ObsBandpass('acs,sbc,%s' % (filt))
    through = bp.throughput
    wave = bp.wave
    return wave, through


def gauss(x, x0, sig):
    return 1 / np.sqrt(2 * np.pi * sig**2) * np.exp(-0.5 * (x - x0)**2 / sig**2)


def spec_norm(m, b, A, wave):
    spec = ((m * wave**b + A * gauss(wave, lya, sigma)))
    # spec = spec * 3e-4 / 100  # /(0.0001400168)*1e-8
    return spec


def spec_calc(m, b, A, norm):
    spec = ((m * wave**b + A * gauss(wave, lya, sigma)))
    spec = spec * norm
    # spec = spec * 3e-4 / 100  # /(0.0001400168)*1e-8
    return spec


def counts_syn(x, m, b, A):

    spec = spec_calc(m, b, A, norm)
    sp = S.ArraySpectrum(wave, spec, waveunits='angstrom', fluxunits='photlam', name='MySource', keepneg=True)
    counts2 = np.zeros((4))
    for i in range((4)):
        obs = S.Observation(sp, S.ObsBandpass('acs,sbc,%s' % (filter_name[i])))
        counts2[i] = obs.countrate()  # effstim('counts')
    return counts2


def counts_syn_no_norm(x, m, b, A, wave):

    spec = spec_norm(m, b, A, wave)
    sp = S.ArraySpectrum(wave, spec, waveunits='angstrom', fluxunits='photlam', name='MySource', keepneg=True)
    counts2 = np.zeros((4))
    for i in range((4)):
        obs = S.Observation(sp, S.ObsBandpass('acs,sbc,%s' % (filter_name[i])))
        counts2[i] = obs.countrate()  # effstim('counts')
    return counts2

# likelihood functions


def lnlike(theta, x, y, yerr):
    m, b, A = theta
    model = counts_syn(x, m, b, A)
    inv_sigma2 = 1.0 / (yerr**2 + model**2)
    return -0.5 * (np.sum((y - model)**2 * inv_sigma2 - np.log(inv_sigma2)))


def lnprior(theta):
    m, b, A = theta
    if 0.1 < m < 10 and -2.0 < b < 2.0 and 0.1 < A < 100.0:  # 3, 0.1
        return 0.0
    return -np.inf


def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


# likelihood functions


def obtain_data(dir1, ind, filters):

    dict_data = OrderedDict.fromkeys(filters)
    dict_err = OrderedDict.fromkeys(filters)
    for f in filters:
        hdu = fits.open(dir1 + 'gal%s_UV_%s_psfmatch.fits' % (ind + 1, f))
        dat = hdu[1].data
        err = hdu[2].data
        err = 1 / np.sqrt(err)
        dict_data[f] = dat
        dict_err[f] = err
    return dict_data, dict_err


def fit_method(method, x, y, yerr, m_true, b_true, A_true):
    nll = lambda *args: -lnlike(*args)

    result = op.minimize(nll, [m_true, b_true, A_true],
                         args=(x, y, yerr), method=method[0], jac='2-point',
                         hess='2-point', options={'disp': True, 'maxfev': 1000, 'xtol': 1e-5, 'ftol': 1e-5})
    return result


# from joblib import Parallel, delayed

# results = Parallel(n_jobs=-1, verbose=51)(delayed(run_me)(i, j)
                                          for i in range(min_x, max_x) for j in range(min_x, max_x))


def new():
    global wave, x, lya, filters, min_x, max_x, norm
    x=np.array([1438.19187645, 1527.99574111, 1612.22929754, 1762.54619064])
    color=['r', 'g', 'b', 'c']

    # this is just to get wavelengths
    wave, through=original_filt('f125lp')
    method=['Powell']
    ind=0
    lya=1215.67 * (1 + red_shift[ind])
    print (lya)
    counts_no_norm=counts_syn_no_norm(x, m_true, b_true, A_true, wave)
    print ('INPUT:-counts from synphot calculation:-\n', counts_no_norm)

    dict_data, dict_err=obtain_data(dir1, ind, filters)
    nx=dict_data[filters[0]].shape[0]
    ny=dict_data[filters[0]].shape[1]
    min_x=621
    max_x=631
    data_m=np.zeros((1024, 1024))
    data_b=np.zeros((1024, 1024))
    data_A=np.zeros((1024, 1024))
    data_status=np.zeros((1024, 1024))

    for i in range(min_x, max_x):
        for j in range(min_x, max_x):
            yerr=np.zeros((4))
            y=np.zeros((4))
            for k in range(4):
                y[k]=dict_data[filters[k]][i, j]
                yerr[k]=dict_err[filters[k]][i, j]
            norm=y[0] / counts_no_norm[0]
            result=fit_method(method, x, y, yerr, m_true, b_true, A_true)
            m_ml, b_ml, A_ml=result["x"]
            data_m[i, j]=m_ml
            data_b[i, j]=b_ml
            data_A[i, j]=A_ml
            # if result["success"] == True:
            #    data_status[i, j] = 1
            # else:
            #    data_status[i, j] = 0
            data_status[i, j]=result["status"]
            print (i, j, m_ml, b_ml, A_ml, result["status"])
            print ("\n-----\n")
            print (result)
            print ("\n------\n")

    fits.writeto('m_full.fits', data = data_m, overwrite = True)
    fits.writeto('b_full.fits', data = data_b, overwrite = True)
    fits.writeto('A_full.fits', data = data_A, overwrite = True)
    fits.writeto('status_full.fits', data = data_status, overwrite = True)


# if __name__ == '__main__':
new()
