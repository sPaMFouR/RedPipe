#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx--------------FIT PIECEWISE LINEAR POLYNOMIAL TO TYPE IIP SN LIGHT CURVES-------------xxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import pwlf
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
file_name = '2016gfy_HCT.dat'
DIR_PHOT = '/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/'
date_maximum = 2457650.00
fit_epoch = 88
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting
# ------------------------------------------------------------------------------------------------------------------- #

def redchisq(ydata, ymod, sd=np.empty(0), n=8):
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
# Read The Apparent Magnitude LC & Lay Out The Guess Parameters For The Fit
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_PHOT + file_name, sep='\s+', comment='#').drop('Date', axis=1)
data_df = data_df.replace('INDEF', np.nan).astype('float64')[data_df['Phase'] < fit_epoch]

list_bands = ['U', 'B', 'V', 'R', 'I']
dict_bands = {'U': [[7.0, 40.0], 1.0, 0.1, 15.2, 19.2], 'B': [[8.5, 40.0], 0.5, 0.05, 16.0, 18.0],
              'V': [[9.7, 39.0], 0.1, 0.01, 15.90, 16.38], 'R': [[9.0, 40.0], 0.1, 0.01, 15.6, 16.14],
              'I': [[9.0, 38.0], 0.1, 0.01, 15.30, 15.97]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Olivares Model & Determine Light Curve Parameters For SN In Study (Technique - 1)
# ------------------------------------------------------------------------------------------------------------------- #
fig, axes = plt.subplots(nrows=2, ncols=5, gridspec_kw={'height_ratios': [4, 1]}, figsize=(32, 9), sharex=True)

for j, band in enumerate(list_bands):
    band_df = data_df[['JD', 'Phase', band, band + 'Err']].copy().dropna()
    x = band_df['Phase'].values
    y = band_df[band].values
    yerr = band_df[band + 'Err'].values

    number_segments = 3
    myPWLF = pwlf.PiecewiseLinFit(x, y, sorted_data=True)
    myPWLF.fit(number_segments)

    xarr = np.linspace(min(x), max(x), num=1000)
    xguess = np.zeros(number_segments - 1)
    xguess[0:] = dict_bands[band][0]

    res = minimize(myPWLF.fit_with_breaks_opt, xguess)

    fit2 = myPWLF.fit_with_breaks([np.min(x)] + dict_bands[band][0] + [np.max(x)])
    yfit2 = myPWLF.predict(xarr)
    fitvar2 = myPWLF.prediction_variance(xarr)
    fitsigma2 = np.sqrt(fitvar2)
    slope2 = myPWLF.slopes * 100

    xbreaks = np.zeros(number_segments + 1)
    xbreaks[0] = np.min(x)
    xbreaks[-1] = np.max(x)
    xbreaks[1:-1] = res.x

    fit = myPWLF.fit_with_breaks(xbreaks)
    yfit = myPWLF.predict(xarr)
    fitvar = myPWLF.prediction_variance(xarr)
    fitsigma = np.sqrt(fitvar)

    print "{0}-Band".format(band)
    print "Number of Parameters: {0}".format(myPWLF.n_parameters)
    print "Best Fit: ", myPWLF.slopes * 100
    print "Manual Fit: ", slope2
    print redchisq(y, myPWLF.predict(x), n=myPWLF.n_parameters)

    axes[0, j].plot(x, y, marker='*', ms=10, ls='', c='k', label='Observed Data')
    axes[0, j].plot(xarr, yfit, ls='-', lw=1.5, c='blue', label='Optimised Fit')
    axes[0, j].plot(xarr, yfit2, ls='--', lw=1.5, c='darkgreen', label='Manual BreakPoints Fit')
    axes[0, j].fill_between(xarr, yfit - (fitsigma * 3), yfit + (fitsigma * 3), alpha=0.4, color='red')
    axes[0, j].fill_between(xarr, yfit2 - (fitsigma2 * 3), yfit2 + (fitsigma2 * 3), alpha=0.2, color='blue')

    axes[0, j].set_ylim(dict_bands[band][4], dict_bands[band][3])
    axes[0, j].yaxis.set_ticks_position('both')
    axes[0, j].xaxis.set_ticks_position('both')
    axes[0, j].yaxis.set_major_locator(MultipleLocator(dict_bands[band][1]))
    axes[0, j].yaxis.set_minor_locator(MultipleLocator(dict_bands[band][2]))
    axes[0, j].xaxis.set_major_locator(MultipleLocator(20))
    axes[0, j].xaxis.set_minor_locator(MultipleLocator(2))
    axes[0, j].tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
    axes[0, j].tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
    axes[0, j].set_title(r'${0}$-Band'.format(band), fontsize=16)

    s1 = (myPWLF.slopes)[1]
    s2 = (myPWLF.slopes)[2]
    yticks = axes[0, j].get_yticks(minor=True)

    for xbreak in xguess:
        axes[0, j].axvline(xbreak, ls='--', lw=0.8, c='k')
        axes[1, j].axvline(xbreak, ls='--', lw=0.8, c='k')
        axes[0, j].text(xbreak + 0.6, yticks[-7], '{:.1f} d'.format(xbreak), rotation=90, color='b', fontsize=13)

    axes[0, j].text(xbreaks[2] + 5, np.mean(yticks), '$s_1$={:.2f}'.format(s1 * 100) +
                    r' mag $\rm (100 d)^{-1}$' + '\n$s_2$={:.2f}'.format(s2 * 100) + r' mag $\rm (100 d)^{-1}$',
                    color='r', fontsize=14)

    axes[1, j].scatter(x, y - myPWLF.predict(x), marker='^', c='k', label=None)

    axes[1, j].set_ylim(-0.13, 0.13)
    axes[1, j].yaxis.set_ticks_position('both')
    axes[1, j].xaxis.set_ticks_position('both')
    axes[1, j].yaxis.set_major_locator(MultipleLocator(0.1))
    axes[1, j].yaxis.set_minor_locator(MultipleLocator(0.02))
    axes[1, j].xaxis.set_major_locator(MultipleLocator(20))
    axes[1, j].xaxis.set_minor_locator(MultipleLocator(2))
    axes[1, j].tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
    axes[1, j].tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
    axes[1, j].set_xlabel('Time Since Explosion [Days]', fontsize=16)

axes[0, 0].set_ylabel('Apparent Magnitude', fontsize=16)
axes[1, 0].set_ylabel(r'Residuals', fontsize=16)
axes[0, 4].legend(markerscale=2, fontsize=12)

fig.subplots_adjust(hspace=0.01, wspace=0.1)
fig.savefig('PLOT_FitPieceWise.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Fit Olivares Model & Determine Light Curve Parameters For SN In Study (Technique - 2)
# # ------------------------------------------------------------------------------------------------------------------- #
# fig, axes = plt.subplots(nrows=2, ncols=4, gridspec_kw={'height_ratios': [4, 1]}, figsize=(30, 10), sharex=True)

# for band in list_bands:
#     band_df = data_df[['JD', 'Phase', band, band + 'Err']].copy().dropna()
# #     band_df = band_df[band_df['JD'] >= date_maximum]

#     x = band_df['Phase'].values
#     y = band_df[band].values
#     yerr = band_df[band + 'Err'].values

#     number_segments = 3
#     myPWLF = pwlf.PiecewiseLinFit(x, y, sorted_data=True)
#     myPWLF.fit(number_segments)

#     xarr = np.linspace(min(x), max(x), num=1000)
#     xguess = np.zeros(number_segments - 1)
#     xguess[0:] = dict_bands[band][0]

#     res = minimize(myPWLF.fit_with_breaks_opt, xguess)

#     fit2 = myPWLF.fit_with_breaks([np.min(x)] + dict_bands[band][0] + [np.max(x)])
#     yfit2 = myPWLF.predict(xarr)
#     fitvar2 = myPWLF.prediction_variance(xarr)
#     fitsigma2 = np.sqrt(fitvar2)
#     slope2 = myPWLF.slopes * 100

#     xbreaks = np.zeros(number_segments + 1)
#     xbreaks[0] = np.min(x)
#     xbreaks[-1] = np.max(x)
#     xbreaks[1:-1] = res.x

#     fit = myPWLF.fit_with_breaks(xbreaks)
#     yfit = myPWLF.predict(xarr)
#     fitvar = myPWLF.prediction_variance(xarr)
#     fitsigma = np.sqrt(fitvar)

#     print "{0}-Band".format(band)
#     print "Number of Parameters: {0}".format(myPWLF.n_parameters)
#     print "Best Fit: ", myPWLF.slopes * 100
#     print "Manual Fit: ", slope2
#     print redchisq(y, myPWLF.predict(x), n=myPWLF.n_parameters)

#     fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 9), sharex=True,
#                                    gridspec_kw={'height_ratios': [4, 1]})

#     ax1.plot(x, y, marker='*', ms=10, ls='', c='k', label='Observed Data')
#     ax1.plot(xarr, yfit, ls='-', lw=1.5, c='blue', label='Optimised Fit')
#     ax1.plot(xarr, yfit2, ls='--', lw=1.5, c='darkgreen', label='Manual BreakPoints Fit')
#     ax1.fill_between(xarr, yfit - (fitsigma * 3), yfit + (fitsigma * 3), alpha=0.4, color='red')
#     ax1.fill_between(xarr, yfit2 - (fitsigma2 * 3), yfit2 + (fitsigma2 * 3), alpha=0.2, color='blue')

#     ax1.set_ylim(dict_bands[band][4], dict_bands[band][3])
#     ax1.legend(markerscale=2, fontsize=12)
#     ax1.yaxis.set_ticks_position('both')
#     ax1.xaxis.set_ticks_position('both')
#     ax1.yaxis.set_major_locator(MultipleLocator(dict_bands[band][1]))
#     ax1.yaxis.set_minor_locator(MultipleLocator(dict_bands[band][2]))
#     ax1.xaxis.set_major_locator(MultipleLocator(20))
#     ax1.xaxis.set_minor_locator(MultipleLocator(2))
#     ax1.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
#     ax1.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
#     ax1.set_title(r'${0}$-Band'.format(band), fontsize=16)
#     ax1.set_ylabel('Apparent Magnitude', fontsize=16)

#     s1 = (myPWLF.slopes)[1]
#     s2 = (myPWLF.slopes)[2]
#     yticks = ax1.get_yticks(minor=True)

#     for xbreak in xbreaks:
#         ax1.axvline(xbreak, ls='--', lw=0.8, c='k')
#         ax2.axvline(xbreak, ls='--', lw=0.8, c='k')
#         ax1.text(xbreak + 0.5, yticks[-7], '{:.1f} d'.format(xbreak), rotation=90, color='brown', fontsize=12)

#     ax1.text((xbreaks[1] + xbreaks[2]) / 2, yticks[0], '$s_1$={:.2f}'.format(s1 * 100), color='r', fontsize=14)
#     ax1.text((xbreaks[2] + xbreaks[3]) / 2, yticks[0], '$s_2$={:.2f}'.format(s2 * 100), color='r', fontsize=14)

#     ax2.scatter(x, y - myPWLF.predict(x), marker='^', c='k', label=None)

#     ax2.set_ylim(-0.13, 0.13)
#     ax2.yaxis.set_ticks_position('both')
#     ax2.xaxis.set_ticks_position('both')
#     ax2.yaxis.set_major_locator(MultipleLocator(0.1))
#     ax2.yaxis.set_minor_locator(MultipleLocator(0.02))
#     ax2.xaxis.set_major_locator(MultipleLocator(20))
#     ax2.xaxis.set_minor_locator(MultipleLocator(2))
#     ax2.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
#     ax2.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
#     ax2.set_xlabel('Time Since Explosion [Days]', fontsize=16)
#     ax2.set_ylabel(r'Residuals', fontsize=16)

#     fig.subplots_adjust(hspace=0.01)
#     fig.savefig('PLOT_FitPieceWise{0}.pdf'.format(band), format='pdf', dpi=2000, bbox_inches='tight')
#     plt.show()
#     plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #
