# !/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx---------------------------PLOT 1-D SPECTRA---------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
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
# file_name = 'Norm_z2012ab_26d.fits'
file_name = 'Norm_z2017hcc_27d.fits'
DIR_SPEC = '/home/avinash/Supernovae_Data/2017hcc/'

dict_17hcc = {'Halpha': (6562.81, r'H$\rm \alpha\ 6563\ \AA$', [1.4, 6562.81, 60, 2.2, 6562.81, 20, 0], '6475,6640'),
              'Hbeta': (4861.63, r'H$\rm \beta\ 4861\ \AA$', [1.1, 4861.63, 40, 1.4, 4861.63, 10, 0], '4800,4920')}
dict_12ab = {'Halpha': (6562.81, r'H$\rm \alpha\ 6563\ \AA$', [1.8, 6562.81, 60, 1.6, 6562.81, 20, 0], '6450,6660'),
             'Hbeta': (4861.63, r'H$\rm \beta\ 4861\ \AA$', [1.1, 4861.63, 40, 1.4, 4861.63, 10, 0], '4800,4920')}

feature = 'Hbeta'
wave_feature, title, guess_fit, wave_range = dict_17hcc[feature]
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
    h1, c1, w1, h2, c2, w2, offset = unc.correlated_values(optpar, covpar)
    func = h1 * unp.exp(-(xarr - c1) ** 2 / (2 * w1 ** 2)) + h2 * unp.exp(-(xarr - c2) ** 2 / (2 * w2 ** 2)) + offset
    fit = unp.nominal_values(func)
    sigma = unp.std_devs(func)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    velarr = np.array([((wave - wave_feature) / (1000 * wave_feature)) * light_speed for wave in xarr])
    ax_obj.plot(velarr, fitlow, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.plot(velarr, fithigh, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.fill_between(velarr, fitlow, fithigh, facecolor=fcolor, alpha=0.5)

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
read_data = spec.read_fits_spectrum1d(DIR_SPEC + file_name)

with fits.open(str(file_name)) as hdulist:
    image_data = hdulist[0].data

wav_data = read_data.dispersion
flux_data = image_data
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Setup The Plot For Fitting The Spectral Feature
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)

ax.plot(wav_data, flux_data, ls='-', lw=1.5, c='k')
ax.set_xlim(min(wav_data), max(wav_data))
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
num_data = np.linspace(1, len(wav_data), len(ax.get_xticks()))

ax2 = ax.twiny()
ax2.set_xticks(num_data)

plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Choose Region For Fitting & Fit The Spectral Feature
# ------------------------------------------------------------------------------------------------------------------- #
# wave_lim = raw_input("Enter The Wavelength Range For The Spectral Line: {0}\r".format(wave_range)) or wave_range
wave_lim = wave_range
wave_lim = wave_lim.split(',')

lower_lim = int(wave_lim[0])
upper_lim = int(wave_lim[1])

usable_flux = flux_data[np.where((wav_data > lower_lim) & (wav_data < upper_lim))]
usable_wav = wav_data[np.where((wav_data > lower_lim) & (wav_data < upper_lim))]

fit_wav = np.arange(min(usable_wav), max(usable_wav), 1)
fit_vel = np.array([((wave - wave_feature) / (1000 * wave_feature)) * light_speed for wave in fit_wav])
usable_vel = np.array([((wave - wave_feature) / (1000 * wave_feature)) * light_speed for wave in usable_wav])

opt, cov = curve_fit(two_gaussians, usable_wav, usable_flux, p0=[guess_fit])
err = np.sqrt(np.diag(cov))
print (opt)

optbroad = [abs(opt[index]) for index in [0, 1, 2, 6]]
covbroad = [abs(err[index]) for index in [0, 1, 2, 6]]
optnar = [abs(opt[index]) for index in [3, 4, 5, 6]]
covnar = [abs(err[index]) for index in [3, 4, 5, 6]]

velbroad = 2.3546 * optbroad[2] * light_speed / wave_feature
errbroad = 2.3546 * covbroad[2] * light_speed / wave_feature
velnar = 2.3546 * optnar[2] * light_speed / wave_feature
errnar = 2.3546 * covnar[2] * light_speed / wave_feature
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plotting The Fit To The Spectral Feature
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111)

ax.scatter(usable_vel, usable_flux, c='k', marker='o', label='Normalised Spectra')
ax.plot(fit_vel, two_gaussians(fit_wav, *opt), lw=1.5, ls='--', c='g', label='Combined Fit')
ax.plot(fit_vel, gaussian(fit_wav, *optbroad), lw=1.5, ls='-', c='r', label='Broad Gaussian')
ax.plot(fit_vel, gaussian(fit_wav, *optnar), lw=1.5, ls=':', c='b', label='Narrow Gaussian')
plot_confintervals(ax, opt, cov, fit_wav)

# ax.set_yticklabels([])
ax.set_xlim(-3.25, 3.25)
ax.set_ylim(min(usable_flux) / 1.03, max(usable_flux) * 1.02)
ax.legend(markerscale=2, fontsize=14, frameon=False)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(1.0))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=16)
ax.tick_params(which='minor', direction='in', width=0.7, length=4, labelsize=16)

ax.set_title(title, fontsize=16)
ax.set_xlabel(r'Doppler Velocity [$\rm \times 10^3\ km\ s^{-1}$]', fontsize=16)
ax.set_ylabel(r'Normalised Flux, $\rm F_{\lambda}$', fontsize=16)

ax.text(0.03, 0.6, r"$\rm FWHM_{Narrow}$ : " + r"{0:.0f}$\pm${1:.0f} km/s".format(velnar, errnar),
        fontsize=14, color='b', transform=ax.transAxes)
ax.text(0.63, 0.6, s=r"$\rm FWHM_{Broad}$ : " + r"{0:.0f}$\pm${1:.0f} km/s".format(velbroad, errbroad),
        fontsize=14, color='r', transform=ax.transAxes)

fig.savefig('PLOT_2017hccFit{0}.pdf'.format(feature), format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
