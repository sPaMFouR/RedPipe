#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxx----------------FIT THE OLIVARES MODEL TO TYPE IIP SN LIGHT CURVES---------------xxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import uncertainties as unc
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
data_fmt = "{0:.3f}"
DIR_PHOT = "/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/"

file_name = '2016gfy_HCT.dat'
date_explosion = 2457641.80
date_maximum = 2457650.0
fit_epoch = 250
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting And Conversions
# ------------------------------------------------------------------------------------------------------------------- #

def fermifunc(t, a0, tPT, w0):
    return -a0 / (1 + np.exp((t - tPT) / w0))


def linefunc(t, p0, tPT, m0):
    return p0 * (t - tPT) + m0


def gaussfunc(t, P, Q, R):
    return -P * np.exp(-((t - Q) / R) ** 2)


def olifunc(t, a0, tPT, w0, p0, m0, P, Q, R):
    return -a0 / (1 + np.exp((t - tPT) / w0)) + p0 * (t - tPT) + m0 - P * np.exp(-((t - Q) / R) ** 2)


def redchisq(ydata, ymod, sd, n=8):
    """
    Args:
        ydata    : Observed data
        ymod     : Model data
          sd     : Uncertainties in ydata
          n      : Number of free parameters in the model
    Returns:
        RedChiSq : Reduced Chi-Square
    """
    ydata = np.array(ydata)
    ymod = np.array(ymod)
    sd = np.array(sd)

    if not sd.any():
        chisq = np.sum((ydata - ymod) ** 2)
    else:
        chisq = np.sum(((ydata - ymod) / sd) ** 2)

    nu = ydata.size - 1 - n
    redchisq = chisq / nu

    return redchisq

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Confidence Intervals
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xarr, fcolor='grey'):
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
    a0, tPT, w0, p0, m0, P, Q, R = unc.correlated_values(optpar, covpar)
    func = -a0 / (1 + unp.exp((xarr - tPT) / w0)) + p0 * (xarr - tPT) + m0 - P * unp.exp(-((xarr - Q) / R) ** 2)
    fit = unp.nominal_values(func)
    sigma = unp.std_devs(func)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr, fitlow, ls='-.', c='k', lw=0.7, alpha=0.3, label='_nolegend_')
    ax_obj.plot(xarr, fithigh, ls='-.', c='k', lw=0.7, alpha=0.3, label='_nolegend_')
    ax_obj.fill_between(xarr, fitlow, fithigh, facecolor=fcolor, alpha=0.3)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Apparent Magnitude LC & Lay Out The Guess Parameters For The Fit
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_PHOT + file_name, sep='\s+', comment='#').drop('Date', axis=1)
data_df = data_df.replace('INDEF', np.nan).astype('float64')[data_df['Phase'] < fit_epoch]

list_bands = ['B', 'V', 'R', 'I']
dict_guess = {'B': [1.970, 99.247, 7.575, 0.0078, 19.7122, 0.8606, -5.0377, 28.1691],
              'V': [1.2449, 110.2644, 2.8179, 0.00989, 18.2976, 0.4402, 78, 29.1894],
              'R': [1.392, 104.2053, 3.590, 0.0087, 17.3458, -0.2, 10, 30],
              'I': [1.220, 105.5630, 3.598, 0.0110, 16.8611, -0.9668, 0.2952, 51.6011]}

# dict_guess = {'B': [1.970, 99.247, 7.575, 0.0078, 19.7122, 0.8606, -5.0377, 28.1691],
#               'V': [1.2449, 110.2644, 2.8179, 0.00989, 18.2976, 0.4402, 78.5989, 29.1894],
#               'R': [1.392, 104.2053, 3.590, 0.0087, 17.3458, -0.2, 10, 20],
#               'I': [1.220, 105.5630, 3.598, 0.0110, 16.8611, -0.9668, 0.2952, 51.6011]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Olivares Model & Determine Light Curve Parameters For SN In Study (Technique - 1)
# ------------------------------------------------------------------------------------------------------------------- #
fig, axes = plt.subplots(nrows=2, ncols=4, gridspec_kw={'height_ratios': [4, 1]}, figsize=(30, 10), sharex=True)

for j in range(axes.shape[1]):
    band = list_bands[j]
    band_df = data_df[['JD', 'Phase', band, band + 'Err']].copy().dropna()
    band_df = band_df[band_df['JD'] >= date_maximum]

    opt, cov = curve_fit(olifunc, band_df['Phase'], band_df[band], sigma=band_df[band + 'Err'], p0=dict_guess[band])

    print (np.round(np.array(opt).astype('float64'), 4))
    print (np.round(np.sqrt(np.diagonal(cov).astype('float64')), 4))
    print (np.round(redchisq(band_df[band], olifunc(band_df['Phase'], *opt), sd=band_df[band + 'Err']), 2))

    jdarr = np.round(np.arange(band_df['Phase'].min(), band_df['Phase'].max(), 0.1), 1)

    axes[0, j].scatter(band_df['Phase'], band_df[band], marker='*', s=40, c='k', label='Observed Data')
    axes[0, j].plot(jdarr, olifunc(jdarr, *opt), c='r', lw=1.2, label='Best Fit')
    axes[0, j].plot(jdarr, fermifunc(jdarr, *opt[0:3]) + opt[4], c='orange', lw=1.2, ls='--', label='Fermi Dirac')
    axes[0, j].plot(jdarr, linefunc(jdarr, *[opt[i] for i in [3, 1, 4]]), c='b', lw=1.2, ls='-.', label='Linear Decay')
    axes[0, j].plot(jdarr, gaussfunc(jdarr, *opt[5:]) + opt[4], c='g', lw=1.2, ls=':', label='Gaussian Peak')
    axes[0, j].axvline(opt[1], ls='--', lw=0.8, c='k')

    axes[0, j].set_ylim(20.5, 15)
    axes[0, j].set_xlim(-2, 242)
    axes[0, j].yaxis.set_ticks_position('both')
    axes[0, j].xaxis.set_ticks_position('both')
    axes[0, j].yaxis.set_major_locator(MultipleLocator(1))
    axes[0, j].yaxis.set_minor_locator(MultipleLocator(0.1))
    axes[0, j].xaxis.set_major_locator(MultipleLocator(50))
    axes[0, j].xaxis.set_minor_locator(MultipleLocator(5))
    axes[0, j].tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
    axes[0, j].tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)
    axes[0, j].text(190, 15.5, '${0}$-Band'.format(band), fontsize=18)

    axes[1, j].scatter(band_df['Phase'], band_df[band] - olifunc(band_df['Phase'], *opt), marker='^', c='k', s=30)
    axes[1, j].axvline(opt[1], ls='--', lw=0.8, c='k')

    axes[1, j].set_ylim(-0.25, 0.25)
    axes[1, j].yaxis.set_ticks_position('both')
    axes[1, j].xaxis.set_ticks_position('both')
    axes[1, j].yaxis.set_major_locator(MultipleLocator(0.2))
    axes[1, j].yaxis.set_minor_locator(MultipleLocator(0.04))
    axes[1, j].xaxis.set_major_locator(MultipleLocator(50))
    axes[1, j].xaxis.set_minor_locator(MultipleLocator(5))
    axes[1, j].tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
    axes[1, j].tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)
    axes[1, j].set_xlabel('Time Since Explosion [Days]', fontsize=16)

    if j != 0:
        axes[0, j].set_yticklabels([])
        axes[1, j].set_yticklabels([])

    plot_confintervals(axes[0, j], opt, cov, jdarr, fcolor='dimgrey')
    axes[0, j].text(opt[1] + 1.5, axes[0, j].get_ylim()[-1] + 0.3, r'Transition Phase', rotation=90, fontsize=12)

axes[0, 3].legend(markerscale=3, fontsize=14, loc='lower right')
axes[0, 0].set_ylabel('Apparent Magnitude', fontsize=18)
axes[1, 0].set_ylabel(r'Residuals', fontsize=18)

fig.subplots_adjust(hspace=0.01, wspace=0.01)
fig.savefig('PLOT_FitOlivares.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Olivares Model & Determine Light Curve Parameters For SN In Study (Technique - 2)
# ------------------------------------------------------------------------------------------------------------------- #

def create_plot(ax1, band_df, opt):
    divider = make_axes_locatable(ax1)
    ax2 = divider.append_axes('bottom', size='20%', pad=0)
    ax1.figure.add_axes(ax2)

    ax1.set_xlim(-2, 242)
    ax2.set_xlim(-2, 242)
    ax1.set_xticklabels([])

    jdarr = np.round(np.arange(band_df['Phase'].min(), band_df['Phase'].max(), 0.1), 1)

    ax1.plot(band_df['Phase'], band_df[band], marker='*', ls='', ms=12, c='k', markerfacecolor='grey',
             label='Observed Data')
    ax1.plot(jdarr, olifunc(jdarr, *opt), c='r', lw=1.2, label='Best Fit')
    ax1.plot(jdarr, fermifunc(jdarr, *opt[0:3]) + opt[4], c='darkgoldenrod', lw=1.4, ls='--', label='Fermi Dirac')
    ax1.plot(jdarr, linefunc(jdarr, *[opt[i] for i in [3, 1, 4]]), c='b', lw=1.4, ls='-.', label='Linear Decay')
    ax1.plot(jdarr, gaussfunc(jdarr, *opt[5:]) + opt[4], c='g', lw=1.4, ls=':', label='Gaussian Peak')
    ax1.axvline(opt[1], ls='--', lw=0.8, c='k')

#     ax1.yaxis.tick_right()
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.yaxis.set_label_position('right')
    ax1.yaxis.set_major_locator(MultipleLocator(1))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.xaxis.set_minor_locator(MultipleLocator(5))
    ax1.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=18, pad=8)
    ax1.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=18)

    ax2.scatter(band_df['Phase'], band_df[band] - olifunc(band_df['Phase'], *opt), marker='^', c='r', s=30)
    ax2.axvline(opt[1], ls='--', lw=0.8, c='k')

    ax2.set_ylim(-0.25, 0.25)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_label_position('right')
    ax2.yaxis.set_major_locator(MultipleLocator(0.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.04))
    ax2.xaxis.set_major_locator(MultipleLocator(50))
    ax2.xaxis.set_minor_locator(MultipleLocator(5))

    ax2.tick_params(axis='y', which='major', direction='in', length=8, width=1.4, labelsize=16)
    ax2.tick_params(axis='x', which='major', direction='in', length=8, width=1.4, labelsize=18)
    ax2.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=18)

    if ax1 == axes[0, 0] or ax1 == axes[1, 0]:
        ax1.set_ylabel('Apparent Magnitude', fontsize=18)
        ax2.set_ylabel('Residuals', color='r', fontsize=18)
    else:
        ax1.set_yticklabels([])
        ax2.set_yticklabels([])

    if ax1 == axes[0, 0] or ax1 == axes[0, 1]:
        ax2.set_xticklabels([])
        ax1.set_ylim(20.5, 15.8)
        ax1.text(200, 16.2, '${0}$-Band'.format(band), color='orangered', fontsize=18)
    else:
        ax1.set_ylim(18.7, 15.2)
        ax1.text(200, 15.5, '${0}$-Band'.format(band), color='orangered', fontsize=18)
        ax2.set_xlabel('Time Since Explosion [Days]', fontsize=18)

    plot_confintervals(ax1, opt, cov, jdarr, fcolor='grey')
    ax1.text(opt[1] + 1.5, ax1.get_ylim()[-1] + 0.2, r'Transition Phase', rotation=90, fontsize=13)


fig, axes = plt.subplots(nrows=len(list_bands) / 2, ncols=2, figsize=(18, 18))

for i in range(axes.shape[0]):
    for j in range(axes.shape[1]):
        band = list_bands[i * axes.shape[0] + j]
        band_df = data_df[['JD', 'Phase', band, band + 'Err']].copy().dropna()
        band_df = band_df[band_df['JD'] >= date_maximum]

        opt, cov = curve_fit(olifunc, band_df['Phase'], band_df[band], sigma=band_df[band + 'Err'],
                             p0=dict_guess[band])
        create_plot(axes[i, j], band_df, opt)

axes[0, 1].legend(markerscale=2.5, fontsize=18, frameon=False, loc='lower left')
fig.subplots_adjust(hspace=0.01, wspace=0.05)
fig.savefig('PLOT_FitOlivares2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
