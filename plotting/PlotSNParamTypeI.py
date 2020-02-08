#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxx--------------------PLOT CORRELATION BETWEEN TYPE I SNe PARAMETERS-------------------xxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import container
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables & Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/Ib_Data/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def line(x, m, c):
    return x * m + c


def calc_sigma(x, covmc):
    return np.sqrt((x ** 2) * covmc[0, 0] + covmc[1, 1] + 2 * x * covmc[0, 1])

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plotting Functions
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xlim=(-15, -20)):
    """
    Plots 1-Sigma and 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xlim    : Limit on the X-Values
    Returns:
        None
    """
    xdata = np.linspace(xlim[0], xlim[1], 1000)
    fit = line(xdata, *optpar)
    sigma = calc_sigma(xdata, covmc=covpar)

    fitlow = fit - sigma
    fithigh = fit + sigma
    fitlow2 = fit - 3 * sigma
    fithigh2 = fit + 3 * sigma

    ax_obj.plot(xdata, fit, linestyle='--', color='k', linewidth=1, label='Our Fit')
    ax_obj.plot(xdata, fitlow, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fithigh, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fitlow2, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')

    ax_obj.fill_between(xdata, fitlow, fithigh, facecolor='grey', alpha=0.7)
    ax_obj.fill_between(xdata, fitlow, fitlow2, facecolor='lightgrey', alpha=0.5)
    ax_obj.fill_between(xdata, fithigh, fithigh2, facecolor='lightgrey', alpha=0.5)

    
def plot_logconfintervals(ax_obj, optpar, covpar, xlim=(-15, -20)):
    """
    Plots 1-Sigma and 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xlim    : Limit on the X-Values
    Returns:
        None
    """
    def exponential(x, m, c):
        return 10 ** line(x, m, c)
    
    xdata = np.linspace(xlim[0], xlim[1], 1000)
    fit = exponential(xdata, *optpar)
    sigma = exponential(xdata, *optpar) * calc_sigma(xdata, covmc=covpar)
    
    fitlow = fit - sigma
    fithigh = fit + sigma
    fitlow2 = fit - 3 * sigma
    fithigh2 = fit + 3 * sigma

    ax_obj.plot(xdata, fit, linestyle='--', color='k', linewidth=1, label='Our Fit')
    ax_obj.plot(xdata, fitlow, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fithigh, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fitlow2, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, linestyle='-.', color='k', linewidth=0.7, alpha=0.7, label='_nolegend_')

    ax_obj.fill_between(xdata, fitlow, fithigh, facecolor='grey', alpha=0.7)
    ax_obj.fill_between(xdata, fitlow, fitlow2, facecolor='lightgrey', alpha=0.6)
    ax_obj.fill_between(xdata, fithigh, fithigh2, facecolor='lightgrey', alpha=0.6)

    
def plot_log3sigmaconfintervals(ax_obj, optpar, covpar, xlim=(-15, -20)):
    """
    Plots 1-Sigma and 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xlim    : Limit on the X-Values
    Returns:
        None
    """
    def exponential(x, m, c):
        return 10 ** line(x, m, c)
    
    xdata = np.linspace(xlim[0], xlim[1], 1000)
    fit = exponential(xdata, *optpar)
    sigma = exponential(xdata, *optpar) * calc_sigma(xdata, covmc=covpar)
    
    fitlow = fit - sigma
    fithigh = fit + sigma
    fitlow2 = fit - 3 * sigma
    fithigh2 = fit + 3 * sigma

    ax_obj.plot(xdata, fit, linestyle='-', color='k', linewidth=1.2, label='Our Fit')
    ax_obj.plot(xdata, fitlow2, linestyle='-.', color='k', linewidth=0.9, alpha=0.7, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, linestyle='-.', color='k', linewidth=0.9, alpha=0.7, label='_nolegend_')

    ax_obj.fill_between(xdata, fit, fitlow2, facecolor='grey', alpha=0.5)
    ax_obj.fill_between(xdata, fit, fithigh2, facecolor='grey', alpha=0.5)
    
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Parameter Data For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_sn = pd.read_csv(DIR_SNe + 'Phot_Data/SESNe.dat', sep='\s+', comment='#', engine='python')

data_sn = data_sn.replace('INDEF', np.nan).set_index('Name', drop=True).dropna(axis=1, how='any').dropna()
data_sn[['BolMag', 'MNi', 'MNiErr']] = data_sn[['BolMag', 'MNi', 'MNiErr']].astype('float64')
data_sn['LogMNi'] = data_sn['MNi'].apply(lambda x: np.log10(x))
# data_sn['LogMNiErr'] = data_sn['MNiErr'] / data_sn['MNi']
data_sn['LogMNiErr+'] = (data_sn['MNi'] + data_sn['MNiErr']).apply(lambda x: np.log10(x)) - data_sn['LogMNi']
data_sn['LogMNiErr-'] = data_sn['LogMNi'] - (data_sn['MNi'] - data_sn['MNiErr']).apply(lambda x: np.log10(x))

data_Ic = data_sn[(data_sn['Type'] == 'Ic') | (data_sn['Type'] == 'Ic-tran')].copy()
data_IcBL = data_sn[data_sn['Type'] == 'Ic-BL'].copy()
data_Ib = data_sn[(data_sn['Type'] == 'Ib') | (data_sn['Type'] == 'Ib-pec')].copy()
data_IIb = data_sn[data_sn['Type'] == 'IIb'].copy()

bolmagmax = data_sn['BolMag'].max()
bolmagmin = data_sn['BolMag'].min()
nimassmax = data_sn['LogMNi'].max()
nimassmin = data_sn['LogMNi'].min()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Absolute Bolometric Magnitude Vs Nickel Mass In Linear Scale
# ------------------------------------------------------------------------------------------------------------------- #
opt, cov = curve_fit(line, data_sn['BolMag'], data_sn['LogMNi'], p0=[-0.5, -7.0])
print(opt, cov)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

plot_confintervals(ax, opt, cov, xlim=(bolmagmax + 1, bolmagmin - 1))

ax.errorbar(-17.3, np.log10(0.09), yerr=0.027 / 0.09, xerr=0.2, color='k', fmt='*', ms=18, capsize=6, label=name_SN)
ax.errorbar(data_Ic['BolMag'], data_Ic['LogMNi'], yerr=data_Ic['LogMNiErr+'], color='red', fmt='D', capsize=5,
            markersize=10, capthick=0.5, elinewidth=0.5, label='Ic')
ax.errorbar(data_IcBL['BolMag'], data_IcBL['LogMNi'], yerr=data_IcBL['LogMNiErr+'], color='darkorange', fmt='D', capsize=5,
            markersize=10, capthick=0.5, elinewidth=0.5, markerfacecolor='None', markeredgewidth=2, label='Ic-BL')
ax.errorbar(data_Ib['BolMag'], data_Ib['LogMNi'], yerr=data_Ib['LogMNiErr+'], color='green', fmt='o', capsize=5,
            markersize=10, capthick=0.5, elinewidth=0.5, label='Ib')
ax.errorbar(data_IIb['BolMag'], data_IIb['LogMNi'], yerr=data_IIb['LogMNiErr+'], color='blue', fmt='o', capsize=5,
            markersize=10, capthick=0.5, elinewidth=0.5, markerfacecolor='None', markeredgewidth=2, label='IIb')

ax.set_ylim(nimassmin - 0.08, nimassmax + 0.08)
ax.set_xlim(bolmagmax + 0.5, bolmagmin - 0.5)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.2))
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels, fontsize=14, markerscale=1.1, frameon=False)

corr = pearsonr(data_sn['BolMag'], data_sn['MNi'])
ax.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_sn.shape[0], corr[0], corr[1]))
ax.set_xlabel(r'Bolometric Magnitude, $\rm M_{peak}$ [mag]', fontsize=16)
ax.set_ylabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=16)

fig.savefig('PLOT_NiVsBolMag.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Absolute Bolometric Magnitude Vs Nickel Mass In Log Scale
# ------------------------------------------------------------------------------------------------------------------- #
opt, cov = curve_fit(line, data_sn['BolMag'], data_sn['LogMNi'], sigma=data_sn['LogMNiErr'], p0=[-0.4, -7])
print(opt, cov)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

plot_log3sigmaconfintervals(ax, opt, cov, xlim=(bolmagmax + 1, bolmagmin - 1))

ax.errorbar(-17.3, 0.09, yerr=0.027, xerr=0.2, color='k', fmt='*', ms=14, capsize=0, elinewidth=1, label=name_SN)
ax.errorbar(data_Ic['BolMag'], data_Ic['MNi'], yerr=data_Ic['MNiErr'], color='red', fmt='D', capsize=0,
            markersize=8, capthick=0.5, elinewidth=0.5, label='Ic')
ax.errorbar(data_IcBL['BolMag'], data_IcBL['MNi'], yerr=data_IcBL['MNiErr'], color='darkorange', fmt='D', capsize=0,
            markersize=8, capthick=0.5, elinewidth=0.5, markerfacecolor='None', markeredgewidth=2, label='Ic-BL')
ax.errorbar(data_Ib['BolMag'], data_Ib['MNi'], yerr=data_Ib['MNiErr'], color='green', fmt='o', capsize=0,
            markersize=8, capthick=0.5, elinewidth=0.5, label='Ib')
ax.errorbar(data_IIb['BolMag'], data_IIb['MNi'], yerr=data_IIb['MNiErr'], color='blue', fmt='o', capsize=0,
            markersize=8, capthick=0.5, elinewidth=0.5, markerfacecolor='None', markeredgewidth=2, label='IIb')

ax.set_yscale('log')
ax.set_ylim(2e-2, 1)
ax.set_xlim(bolmagmax + 0.4, bolmagmin - 0.4)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels, fontsize=14, markerscale=1.4, frameon=False)

corr = pearsonr(data_sn['BolMag'], data_sn['MNi'])
ax.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_sn.shape[0], corr[0], corr[1]))
ax.set_xlabel(r'Bolometric Magnitude, $\rm M_{peak}$ [mag]', fontsize=16)
ax.set_ylabel(r'$\rm M_{Ni}\ [M_{\odot}]$', fontsize=16)

fig.savefig('PLOT_NiVsBolMagSemiLogScale.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #

