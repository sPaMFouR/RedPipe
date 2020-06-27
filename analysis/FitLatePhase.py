#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx----------FIT RADIOACTIVE DECAY PHASE TO ACCOUNT FOR THE GAMMA-RAY TRAPPING-----------xxxxxxxxxxxxxxx #
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
from scipy.interpolate import CubicSpline
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
file_applc = '2016gfy_HCT.dat'
DIR_PHOT = '/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/'
DIR_CODE = '/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/'
kgamma = 0.033
tCo = 111.4
tNi = 8.8

name_SN = '2016gfy'
date_explosion = 2457641.4
EBV_mag = 0.21
EBV_err = 0.05
dist_val = 29.64
dist_err = 2.65
fit_epoch = 130
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting And Conversions
# ------------------------------------------------------------------------------------------------------------------- #

def calc_Ek(Menv, v):
    return (3 / 10) * Menv * (v ** 2)


def calc_timescale(Menv, v):
    global kgamma
    return (9 * kgamma * (Menv ** 2) / (40 * np.pi * calc_Ek(Menv, v))) ** 0.5


def line(x, m, c):
    return m * x + c


def latephasefunc(t, mNi, tc):
    global tCo, tNi
    return 1.41e43 * mNi * (np.exp(-t / tCo) - np.exp(-t / tNi)) * (1 - np.exp(-(tc / t) ** 2))


def calc_lum(flux):
    val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
    lum = val * dist_val ** 2
    lumerr = val * ((dist_val + dist_err) ** 2 - dist_val ** 2)

    return lum, lumerr


def calc_magflux(mag, err, band):
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

def plot_confintervals(ax_obj, fit, sigma, xarr, fcolor='black'):
    """
    Plots 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        fit     : Functional fit for which the confidence interval is to be plotted
        sigma   : Deviation of the fit for which the confidence interval is to be plotted
        xarr    : Array of X-Values over which confidence intervals are to be plotted
        fcolor  : Fill color for the confidence intervals
    Returns:
        None
    """
    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr, fitlow, ls='-.', c='k', lw=0.8, alpha=0.8, label='_nolegend_')
    ax_obj.plot(xarr, fithigh, ls='-.', c='k', lw=0.8, alpha=0.8, label='_nolegend_')
    ax_obj.fill_between(xarr, fitlow, fithigh, facecolor=fcolor, alpha=0.2)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS and Other Type II SNe
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Apparent Magnitude Optical Light Curve Data
# Convert The Data From A Column-Wise (Based on FILTERs) To A Row-Wise Representation
# ------------------------------------------------------------------------------------------------------------------- #
data_app = pd.read_csv(file_applc, sep='\s+', engine='python').drop(['Date', 'Phase'], axis=1)
data_app = data_app.replace('INDEF', np.nan).astype('float64')[data_app['JD'] > date_explosion + fit_epoch]
data_app = data_app.interpolate(method='linear', limit=2).round(2).set_index('JD')
list_filters = [band for band in data_app.columns.values if 'Err' not in band]

data_arr = data_app.values
size = data_arr.shape

list_jd = np.repeat(data_app.index.values, (len(list_filters)))
data_arr = np.reshape(data_arr, [size[0] * len(list_filters), 2])

master_df = pd.DataFrame(data_arr, index=list_jd, columns=['FMAG', 'FERR'], dtype='float64')
master_df.index.name = 'JD'
master_df = master_df.reset_index(drop=False)
master_df['FILTER'] = list_filters * size[0]
master_df = master_df.reindex(columns=master_df.columns.tolist() + ['Flux', 'FluxErr'])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate The Apparent Flux In Each Photometric Band
# Compute The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
for index, band in master_df['FILTER'].items():
    magflux = calc_magflux(mag=master_df.loc[index, 'FMAG'], err=master_df.loc[index, 'FERR'], band=band)
    master_df.loc[index, ['Flux', 'FluxErr']] = magflux

master_df = master_df.set_index('JD')

dict_val = {}
for index, row in master_df.iterrows():
    if index not in dict_val.keys():
        dict_val[index] = {}
    if row['FILTER'] not in dict_val[index]:
        dict_val[index][row['FILTER']] = []
    dict_val[index][row['FILTER']].append(row['Flux'])

for (day, dict_date) in dict_val.items():
    for (band, list_flux) in dict_date.items():
        dict_val[day][band] = float(np.mean(list_flux))

flux_df = pd.DataFrame(dict_val)
flux_df.index = [filter_df.loc[band, 'CentreWave'] for band in flux_df.index.values]
flux_df = flux_df.sort_index()

dict_flux = {}
for caljd in flux_df.columns.values:
    series = flux_df[caljd].dropna().apply(lambda x: float(x))

    if caljd > 110:
        wave_data = np.linspace(4000, 9200, 1000)
    else:
        wave_data = np.linspace(3100, 9200, 1000)

    spline = CubicSpline(series.index.values.tolist(), series.values.tolist(), bc_type='natural', extrapolate=True)

    flux_data = spline(wave_data)
    flux_data[flux_data < 0] = 0
    netflux = np.trapz(flux_data, wave_data)

    dict_flux[caljd] = {}
    dict_flux[caljd]['Flux'] = netflux
    dict_flux[caljd]['Lum'] = calc_lum(netflux)[0]
    dict_flux[caljd]['LumErr'] = calc_lum(netflux)[1]

#     print caljd
#     print caljd - date_explosion
#     fig_temp = plt.figure(figsize=(10, 8))
#     ax = fig_temp.add_subplot(111)
#     ax.plot(series.index.values, series.values, 'o', label='Data Points')
#     ax.plot(wave_data, flux_data, 'r-', label='CubicSpline')
#     ax.legend()
#     ax.grid()
#     plt.show()
#     plt.close(fig_temp)

data_df = pd.DataFrame(dict_flux).T
data_df.index.name = 'JD'
data_df = data_df.reset_index()
data_df['Phase'] = (data_df['JD'] - date_explosion).round(2)
data_df = data_df[['JD', 'Phase', 'Lum', 'LumErr']].dropna()
data2_df = data_df[data_df['Phase'] < 250].dropna()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Late Phase Model & Determine Nickel Mass For SN In Study
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

jdarr = np.arange(data_df['Phase'].min(), data_df['Phase'].max(), 0.1)
jdarr2 = np.arange(data2_df['Phase'].min(), data2_df['Phase'].max(), 0.1)

opt, cov = curve_fit(latephasefunc, data_df['Phase'], data_df['Lum'], sigma=data_df['LumErr'], p0=[0.045, 100])
opt2, cov2 = curve_fit(line, data2_df['Phase'], data2_df['Lum'], sigma=data2_df['LumErr'], p0=[1e39, 1e41])

mNi, tc = unc.correlated_values(opt, cov)
latefunc = 1.41e43 * mNi * (unp.exp(-jdarr / tCo) - np.exp(-jdarr / tNi)) * (1 - unp.exp(-(tc / jdarr) ** 2)) / 1e40
fit = unp.nominal_values(latefunc)
sigma = unp.std_devs(latefunc)

m, c = unc.correlated_values(opt2, cov2)
linefunc = line(jdarr2, m, c) / 1e40
fit2 = unp.nominal_values(linefunc)
sigma2 = unp.std_devs(linefunc)

print (np.round(opt, 4))
print (np.sqrt(np.diagonal(cov)))
print (redchisq(data_df['Lum'], latephasefunc(data_df['Phase'], *opt), sd=data_df['LumErr']))

print (np.round(opt2, 4))
print (np.sqrt(np.diagonal(cov2)))
print (redchisq(data2_df['Lum'], line(data2_df['Phase'], *opt2), sd=data2_df['LumErr']))

ax.plot(data_df['Phase'], data_df['Lum'] / 1e40, marker='*', ls='', ms=12, markerfacecolor='dodgerblue', c='k',
        label='Observed Data')
ax.plot(jdarr, latephasefunc(jdarr, *opt) / 1e40, c='k', lw=1.5, ls='--', alpha=0.8, label='Late Phase Fit')
ax.plot(jdarr, line(jdarr, *opt2) / 1e40, c='orangered', lw=2.0, ls='-', label='Linear Fit')
ax.errorbar(data_df['Phase'], data_df['Lum'] / 1e40, yerr=data_df['LumErr'] / 1e40, c='k', lw=1.2, ls='', capsize=4,
            capthick=1, alpha=0.7, label='_nolegend_')

ax.text(140, 0.5, s='\nDecline Rate = {0:.2e}'.format(100 * opt2[0]) + r' $\rm erg\ s^{-1}\ (100\ d)^{-1}$',
        color='mediumblue', fontsize=15)
ax.text(300, 9, color='k', fontsize=16, s=r'$\rm\ \ t_c\ \sim$' + ' {0:.0f} d\n'.format(opt[1]) +
        r'$\rm M_{Ni}$' + ' = {0:.3f} '.format(opt[0]) + r'$\rm M_{\odot}$')

plot_confintervals(ax, fit, sigma, jdarr)
# plot_confintervals(ax, fit2, sigma2, jdarr2, fcolor='blue')

ax.set_ylim(0, 14)
ax.set_xlim(135, 390)
ax.legend(markerscale=1.6, frameon=False, fontsize=16)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(3))
ax.yaxis.set_minor_locator(MultipleLocator(0.3))
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)

ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax.set_ylabel(r'Pseudo-Bolometric Luminosity [x $\rm 10^{40}\ erg\ s^{-1}$]', fontsize=16)

fig.subplots_adjust(hspace=0.01, wspace=0.01)
fig.savefig('PLOT_FitLatePhase.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
