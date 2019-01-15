#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx------------FIT COWEN MODEL TO THE EARLY PHASE LIGHT CURVE IN TYPE IIP SNe------------xxxxxxxxxxxxxxx #
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
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
data_fmt = "{0:.3f}"
file_name = '2016gfy_HCT.dat'
DIR_PHOT = '/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/'
DIR_CODE = '/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/'
jdfit_epoch = 2457670.0
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of SN In Study
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
EBV_mag = 0.21
EBV_err = 0.11
date_discovery = 2457644.60
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting And Conversions
# ------------------------------------------------------------------------------------------------------------------- #

def earlyrisefunc(t, A, B, t0, alpha, beta):
    return A * ((t - t0) ** beta) / (np.exp(B * ((t - t0) ** alpha) - 1))


def cowenfunc(t, a1, a2, t0, a3):
    return a1 * ((t - t0) ** 1.6) / (np.exp(a2 * ((t - t0) ** 0.5) - 1)) + a3 * ((t - t0) ** 2)


def calc_magflux(data, band):
    mag = data[band]
    err = data[band + 'Err']
    rlambda = filter_df.loc[band, 'RLambda']
    zp = filter_df.loc[band, 'ZeroPoint']

    flux = 10 ** (-0.4 * (mag - rlambda * EBV_mag + zp + 21.100))
    fluxerr = abs(flux - 10 ** (-0.4 * (mag + err - rlambda * EBV_mag + zp + 21.100)))

    return flux, fluxerr


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

def plot_confintervals(ax_obj, optpar, covpar, xarr, offset, fcolor='grey'):
    """
    Plots 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xarr    : Array of X-Values over which confidence intervals are to be plotted
        offset  : Offset epoch in JD
        fcolor  : Fill color for the confidence intervals
    Returns:
        None
    """
    a1, a2, t0, a3 = unc.correlated_values(optpar, covpar)
    func = a1 * ((xarr - t0) ** 1.6) / (unp.exp(a2 * ((xarr - t0) ** 0.5) - 1)) + a3 * ((xarr - t0) ** 2)
    fit = unp.nominal_values(func)
    sigma = unp.std_devs(func)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr + offset, fitlow, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xarr + offset, fithigh, ls='-.', c='k', lw=0.7, alpha=0.5, label='_nolegend_')
    ax_obj.fill_between(xarr + offset, fitlow, fithigh, facecolor=fcolor, alpha=0.2)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS and Other Type II SNe
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
list_filters = filter_df.index.tolist()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Apparent Magnitude LC & Lay Out The Guess Parameters For The Fit
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_PHOT + file_name, sep='\s+', comment='#').drop('Date', axis=1)
data_df = data_df.replace('INDEF', np.nan).astype('float64')[data_df['JD'] < jdfit_epoch]

list_bands = ['U', 'B', 'V', 'R', 'I']
dict_bands = {'U': [[3.2, 1.6, date_discovery - 3, 0.0001], '*', 'navy'],
              'B': [[5, 2, date_discovery - 3, 0.0001], 'o', 'b'],
              'V': [[2, 2, date_discovery - 3, 0.005], 'P', 'g'],
              'R': [[2, 2, date_discovery - 3, 0.005], 'x', 'r'],
              'I': [[1, 2, date_discovery - 3, 0.005], 's', 'darkorange']}

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Cowen Model & Determine Rise Time Parameters For SN In Study
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 9))
ax = fig.add_subplot(111)

maxepochs = []
expepochs = []
expepochs_err = []
for band in list_bands:
    band_df = data_df[['JD', 'Phase', band, band + 'Err']].copy().dropna()
    band_df = band_df.reindex(band_df.columns.tolist() + ['Flux', 'FluxErr'], axis=1)

    for index, row in band_df.iterrows():
        data_magflux = calc_magflux(data=row, band=band)
        row['Flux'] = data_magflux[0] / 1e-15
        row['FluxErr'] = data_magflux[1] / 1e-15

    opt, cov = curve_fit(cowenfunc, band_df['JD'], band_df['Flux'], sigma=band_df['FluxErr'], p0=dict_bands[band][0])
    offset = -date_discovery

    print np.round(opt, 4)
    print np.sqrt(np.diagonal(cov))
#     print redchisq(band_df['Flux'], cowenfunc(band_df['JD'], *opt), sd=band_df['FluxErr'])

    jdarr = np.arange(opt[2], band_df['JD'].max(), 0.1)
    jdmax = jdarr[np.where(cowenfunc(jdarr, *opt) == max(cowenfunc(jdarr, *opt)))] + offset

    maxepochs.append(jdmax)
    expepochs.append(opt[2])
    expepochs_err.append(np.sqrt((np.diagonal(cov)[2])))

    ax.axvline(jdmax, c=dict_bands[band][2], ls='--', lw=1.2)
    ax.text(jdmax + 0.1, 0.5, color=dict_bands[band][2], s=band + r'$\rm _{Max}$', rotation=90, fontsize=12)
    ax.plot(jdarr + offset, cowenfunc(jdarr, *opt), c='k', lw=1.2, alpha=0.8, label='_nolegend_')
    ax.scatter(band_df['JD'] + offset, band_df['Flux'], marker=dict_bands[band][1], s=100, c=dict_bands[band][2],
               label=band + r' ,$\rm t_{Rise}$ = ' + '{0:.1f} d'.format(float(jdmax - opt[2] + date_discovery)))
    ax.errorbar(band_df['JD'] + offset, band_df['Flux'], yerr=band_df['FluxErr'], c='k', lw=1.2,
                ls='', capsize=4, capthick=2, alpha=0.7, label='_nolegend_')

    plot_confintervals(ax, opt, cov, jdarr, offset=offset, fcolor=dict_bands[band][2])

ax.set_ylim(0, 7.1)
ax.set_xlim(-6, 26)
ax.legend(markerscale=1.6, frameon=False, fontsize=16)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)

ax.set_xlabel('Time Since Discovery [Days]', fontsize=16)
ax.set_ylabel('Flux [Arbitrary Units]', fontsize=16)

exp_mean = np.mean(expepochs)
exp_err = (np.std(expepochs) ** 2 + (np.sum(np.array(expepochs_err) ** 2)) / len(expepochs) ** 2) ** 0.5
ax.set_title('Mean Explosion Epoch = {0:.1f}$\pm${1:.1f}'.format(exp_mean, exp_err), fontsize=16)

print maxepochs - exp_mean + date_discovery

fig.subplots_adjust(hspace=0.01, wspace=0.01)
# fig.savefig('PLOT_FitRiseTime.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
