# !/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx---------------------------PLOT 1-D SPECTRA---------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import uncertainties as unc
from astropy.io import fits
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
light_speed = 2.99792458e5  # In km/s
norm_factor = 1e-15
file_name = 'z_NGC2276_HostHII.dat'
DIR_SPEC = '/home/avinash/Supernovae_Data/2016gfy/NGC2276/HCTSpec/Iter6/'

dict_17hcc = {'Halpha': (6562.81, r'H$\rm \alpha\ 6563\ \AA$', [1.4, 6562.81, 60, 2.2, 6562.81, 20, 0], '6475,6640'),
              'Hbeta': (4861.63, r'H$\rm \beta\ 4861\ \AA$', [1.1, 4861.63, 40, 1.4, 4861.63, 10, 0], '4800,4920')}
dict_12ab = {'Halpha': (6562.81, r'H$\rm \alpha\ 6563\ \AA$', [1.8, 6562.81, 60, 1.6, 6562.81, 20, 0], '6450,6660'),
             'Hbeta': (4861.63, r'H$\rm \beta\ 4861\ \AA$', [1.1, 4861.63, 40, 1.4, 4861.63, 10, 0], '4800,4920')}
dict_16gfy = {'Halpha': (6562.81, r'H$\rm \alpha\ 6563\ \AA$', [0.4, 6548, 5, 2.3, 6564, 5, 0.7, 6584, 5, 1],
              '6500,6620')}

feature = 'Halpha'
wave_feature, title, guess_fit, wave_range = dict_16gfy[feature]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Confidence Intervals
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xarr, fcolor='orange'):
    """
    Plots 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xarr    : Array of X-Values over which confidence intervals are to be plotted
        fcolor  : Fill color for the confidence intervals
    Returns:
        None
    """
    h1, c1, w1, h2, c2, w2, h3, c3, w3, offset = unc.correlated_values(optpar, covpar)
    func = h1 * unp.exp(-(xarr - c1) ** 2 / (2 * w1 ** 2)) + h2 * unp.exp(-(xarr - c2) ** 2 / (2 * w2 ** 2)) + h3 * unp.exp(-(xarr - c3) ** 2 / (2 * w3 ** 2)) + offset
    fit = unp.nominal_values(func)
    sigma = unp.std_devs(func)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr, fitlow, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xarr, fithigh, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.fill_between(xarr, fitlow, fithigh, facecolor=fcolor, alpha=0.3)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting
# ------------------------------------------------------------------------------------------------------------------- #

def gaussian(x, height, center, width, offset):
    return height * np.exp(-(x - center) ** 2 / (2 * width ** 2)) + offset


def two_gaussians(x, h1, c1, w1, h2, c2, w2, offset):
    return gaussian(x, h1, c1, w1, offset=0) + gaussian(x, h2, c2, w2, offset=0) + offset


def three_gaussians(x, h1, c1, w1, h2, c2, w2, h3, c3, w3, offset):
    return (gaussian(x, h1, c1, w1, offset=0) + gaussian(x, h2, c2, w2, offset=0) +
            gaussian(x, h3, c3, w3, offset=0) + offset)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Data From The Normalised Spectra
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_SPEC + file_name, sep='\s+', header=None, names=['Wavelength', 'Flux'], comment='#')
wav_data = np.array(data_df['Wavelength'])
flux_data = np.array(data_df['Flux']) / norm_factor
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Choose Region For Fitting & Fit The Spectral Feature
# ------------------------------------------------------------------------------------------------------------------- #
# wave_lim = raw_input("Enter The Wavelength Range For The Spectral Line: {0}\r".format(wave_range)) or wave_range
wave_lim = wave_range.split(',')
lower_lim = int(wave_lim[0])
upper_lim = int(wave_lim[1])

usable_flux = flux_data[np.where((wav_data > lower_lim) & (wav_data < upper_lim))]
usable_wav = wav_data[np.where((wav_data > lower_lim) & (wav_data < upper_lim))]

fit_wav = np.arange(min(usable_wav), max(usable_wav), 1)
fit_vel = np.array([((wave - wave_feature) / (1000 * wave_feature)) * light_speed for wave in fit_wav])
usable_vel = np.array([((wave - wave_feature) / (1000 * wave_feature)) * light_speed for wave in usable_wav])

opt, cov = curve_fit(three_gaussians, usable_wav, usable_flux, p0=[guess_fit])
err = np.sqrt(np.diag(cov))

optn1 = [abs(opt[index]) for index in [0, 1, 2, 9]]
covn1 = [abs(err[index]) for index in [0, 1, 2, 9]]
optha = [abs(opt[index]) for index in [3, 4, 5, 9]]
covha = [abs(err[index]) for index in [3, 4, 5, 9]]
optn2 = [abs(opt[index]) for index in [6, 7, 8, 9]]
covn2 = [abs(err[index]) for index in [6, 7, 8, 9]]

print optn1
print optha
print optn2
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plotting The Fit To The Spectral Feature
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)

ax.scatter(usable_wav, usable_flux, c='k', marker='x', s=40, label='Original Spectra')
ax.plot(fit_wav, three_gaussians(fit_wav, *opt), lw=2.5, ls='-', c='r', label='Triple Gaussian Fit')
ax.plot(fit_wav, gaussian(fit_wav, *optn1), lw=1.5, ls='--', c='b', label='_nolegend_')
ax.plot(fit_wav, gaussian(fit_wav, *optha), lw=1.5, ls='--', c='darkorange', label='_nolegend_')
ax.plot(fit_wav, gaussian(fit_wav, *optn2), lw=1.5, ls='--', c='g', label='_nolegend_')

plot_confintervals(ax, opt, cov, fit_wav, fcolor='grey')

ax.set_yticklabels([])
ax.set_xlim(6520, 6610)
ax.legend(markerscale=2, fontsize=14, frameon=False)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=16)
ax.tick_params(which='minor', direction='in', width=0.7, length=4, labelsize=16)

ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=18)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$]', fontsize=18)

ax.text(0.29, 0.37, r"$\rm [N\,II]\ 6548\ \AA$", fontsize=16, color='b', rotation=90, transform=ax.transAxes)
ax.text(0.49, 0.6, r"$\rm H\alpha$", fontsize=18, color='darkorange', rotation=90, transform=ax.transAxes)
ax.text(0.71, 0.49, r"$\rm [N\,II]\ 6584\ \AA$", fontsize=16, color='g', rotation=90, transform=ax.transAxes)

fig.savefig('PLOT_SpecHostHalpha.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
