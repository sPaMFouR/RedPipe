#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx------------------------PLOT THE TYPE II SNE PARAMETERS-------------------------xxxxxxxxxxxxxxxxxx #
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
# Details Of SN In Study & Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
OPTd = 90
OPTdErr = 5
s1 = 0.94
s1Err = 0.02
s2 = 0.12
s2Err = 0.02
Mv50 = -16.74
Mv50Err = 0.24
MMax = -17.06
MMaxErr = 0.24
MNi = 0.033
MNiErr = 0.003
V50 = 3913
V50Err = 320

# All stars with Rc-Ic < 1.15 [Jordi et al. 2006]
#         Transformation                RMS residual
#     u-g    =    1.28*(U-B)   + 1.13      0.06
#     g-r    =    1.02*(B-V)   - 0.22      0.04
#     r-i    =    0.91*(Rc-Ic) - 0.20      0.03
#     r-z    =    1.72*(Rc-Ic) - 0.41      0.03
#     g      =    V + 0.60*(B-V) - 0.12    0.02
#     r      =    V - 0.42*(B-V) + 0.11    0.03

UB30 = 0.00
UB30Err = 0.09
ug30 = 1.28 * UB30 + 1.13
ug30Err = (0.06 ** 2 + (1.28 * UB30Err) ** 2) ** 0.5

BV70 = 1.35
BV70Err = 0.055
gr70 = 1.02 * BV70 - 0.22
gr70Err = (0.04 ** 2 + (1.02 * BV70Err) ** 2) ** 0.5

BV15 = 0.125
BV15Err = 0.055
gr15 = 1.02 * BV15 - 0.22
gr15Err = (0.04 ** 2 + (1.02 * BV15Err) ** 2) ** 0.5

VHa70 = 7000
VHa70Err = 160

name_hostgal = 'NGC 2276'
# V 11.30+-0.13
# B 11.75+-0.13
# U 11.61+-0.14
distmod = 32.36
distmodErr = 0.18
hostMv = 11.30 - distmod
hostMvErr = 0.20
hostZ = 8.61
hostZErr = 0.18

DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def line(x, m, c):
    return x * m + c


def calc_sigma(x, covmc):
    return np.sqrt(covmc[0, 0] * (x ** 2) + covmc[1, 1] + 2 * x * covmc[0, 1])

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Confidence Intervals And Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xlim=(-3.5, 0.0), alpha=0.6, fcolor='dimgrey'):
    """
    Plots 1-Sigma and 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xlim    : Limit on the X-Values
        alpha   : Transparency for the color used for the confidence intervals
        fcolor  : Fill color for the confidence intervals
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

    ax_obj.plot(xdata, fit, ls='--', c='k', label='Our Fit')
    ax_obj.plot(xdata, fitlow, ls='-.', c='k', lw=0.6, alpha=alpha, label='_nolegend_')
    ax_obj.plot(xdata, fithigh, ls='-.', c='k', lw=0.6, alpha=alpha, label='_nolegend_')
    ax_obj.plot(xdata, fitlow2, ls='-.', c='k', lw=0.6, alpha=alpha, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, ls='-.', c='k', lw=0.6, alpha=alpha, label='_nolegend_')

    ax_obj.fill_between(xdata, fitlow, fithigh, facecolor=fcolor, alpha=alpha)
    ax_obj.fill_between(xdata, fitlow, fitlow2, facecolor=fcolor, alpha=alpha / 3)
    ax_obj.fill_between(xdata, fithigh, fithigh2, facecolor=fcolor, alpha=alpha / 3)


def set_plotparams(ax_obj, legend=True, markerscale=1.5, markerfont=13):
    """
    Sets plot parameters to the axes object 'ax_obj'.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        legend      : Boolean stating whether legend is to be displayed in the plot
        markerscale : Scaled marker size to be used in the legend
        markerfont  : Font of the labels in the legend
    Returns:
        None
    """
    if legend:
        handles, labels = ax_obj.get_legend_handles_labels()
        handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
        ax_obj.legend(handles, labels, fontsize=markerfont, markerscale=markerscale, frameon=False)

    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
    ax_obj.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Parameter Data For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_sn = pd.read_csv(DIR_SNe + 'Param_Data/MNiSample.dat', sep='\s+', comment='#', engine='python')
data_sn = data_sn.replace('INDEF', np.nan).set_index('Name', drop=True).astype('float64')

data_sn['LogMNi'] = data_sn['MNi'].apply(lambda x: np.log10(x))
data_sn['LogMNiErr+'] = data_sn['MNiErr+'] / data_sn['MNi']
data_sn['LogMNiErr-'] = data_sn['MNiErr-'] / data_sn['MNi']
data_sn['LogV50'] = data_sn['V50'].apply(lambda x: np.log10(x))
data_sn['LogV50Err'] = data_sn['V50Err'] / data_sn['V50']

data_spiro = data_sn.iloc[0:15].copy()
data_hamuy = data_sn.iloc[15:38].copy()
data_valenti = data_sn.iloc[38:46].copy()
data_others = data_sn.iloc[46:54].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Light Curve Parameters of Type II SNe (Anderson et al.(2014))
# ------------------------------------------------------------------------------------------------------------------- #
data_and = pd.read_csv(DIR_SNe + 'Param_Data/AndersonSample.dat', sep='\s+', comment='#', engine='python')
data_and = data_and.replace('INDEF', np.nan).set_index('Name', drop=True).astype('float64')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Galaxy Properties Data
# ------------------------------------------------------------------------------------------------------------------- #
data_gal = pd.read_csv(DIR_SNe + 'Param_Data/TremontiGalZSample.dat', sep='\s+', comment='#', engine='python')
data_gal = data_gal[['MB1', 'MB2', 'OH']].sort_values('MB2').replace('INDEF', np.nan).astype('float64').dropna()
data_gal = data_gal[(data_gal['MB2'] > -30) & (data_gal['MB2'] < -12)]

data_sngal = pd.read_csv(DIR_SNe + 'Param_Data/PrietoSNZSample.dat', sep='\s+', comment='#', engine='python')
data_sngal = data_sngal.replace('INDEF', np.nan).astype('float64')[data_sngal.index < 135]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Color Correlation Data For Type II SNe (De Jaeger et al.(2018))
# ------------------------------------------------------------------------------------------------------------------- #
data_Mmaxug30 = pd.read_csv(DIR_SNe + 'Param_Data/DeJaeger_u-g30VsMMax.dat', sep='\s+', comment='#')
data_Mmaxug30 = data_Mmaxug30.replace('INDEF', np.nan).astype('float64')

data_s2gr70 = pd.read_csv(DIR_SNe + 'Param_Data/DeJaeger_g-r70Vss2.dat', sep='\s+', comment='#')
data_s2gr70 = data_s2gr70.replace('INDEF', np.nan).astype('float64')

data_s1gr15 = pd.read_csv(DIR_SNe + 'Param_Data/DeJaeger_g-r15Vss1.dat', sep='\s+', comment='#')
data_s1gr15 = data_s1gr15.replace('INDEF', np.nan).astype('float64')

data_VHagr70 = pd.read_csv(DIR_SNe + 'Param_Data/DeJaeger_g-r70VsVHa70.dat', sep='\s+', comment='#')
data_VHagr70 = data_VHagr70.replace('INDEF', np.nan).astype('float64')
data_VHagr70[['VHa70', 'VHa70Err']] = data_VHagr70[['VHa70', 'VHa70Err']].apply(lambda x: x / 1000.)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Equivalent Width Data From Dessart et al. (2013)
# ------------------------------------------------------------------------------------------------------------------- #
data_ewz2m3 = pd.read_csv(DIR_SNe + 'EW_Data/m15z2m3_FeII.dat', sep='\s+', comment='#')
data_ewz8m3 = pd.read_csv(DIR_SNe + 'EW_Data/m15z8m3_FeII.dat', sep='\s+', comment='#')
data_ewz2m2 = pd.read_csv(DIR_SNe + 'EW_Data/m15z2m2_FeII.dat', sep='\s+', comment='#')
data_ewz4m2 = pd.read_csv(DIR_SNe + 'EW_Data/m15z4m2_FeII.dat', sep='\s+', comment='#')
data_ewmlt1 = pd.read_csv(DIR_SNe + 'EW_Data/m15mlt1.dat', sep='\s+', comment='#')
data_ewmlt3 = pd.read_csv(DIR_SNe + 'EW_Data/m15mlt3.dat', sep='\s+', comment='#')

data_ewsn = pd.read_csv(DIR_SNe + 'EW_Data/2016gfy_FeII.dat', sep='\s+', comment='#')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Absolute V-Band Magnitude Vs Nickel Mass At Mid-Plateau
# ------------------------------------------------------------------------------------------------------------------- #
data_ni = data_sn[['Mv', 'MvErr', 'LogMNi', 'LogMNiErr+', 'LogMNiErr-']].dropna()[:-1]
opt, cov = curve_fit(line, data_ni['LogMNi'], data_ni['Mv'], sigma=data_ni['MvErr'], p0=[-2.5, -20.0])
print (opt, cov)

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

plot_confintervals(ax, opt, cov, xlim=(-3.5, 0.0))
ax.errorbar(data_hamuy['LogMNi'], data_hamuy['Mv'], yerr=data_hamuy['MvErr'], c='darkorange', fmt='^', capsize=5,
            ms=7, capthick=0.5, xerr=np.vstack((data_hamuy['LogMNiErr+'], data_hamuy['LogMNiErr-'])),
            elinewidth=0.8, label='Hamuy (2003)')
ax.errorbar(data_spiro['LogMNi'], data_spiro['Mv'], yerr=data_spiro['MvErr'], c='b', fmt='<', capsize=5,
            ms=7, capthick=0.5, xerr=np.vstack((data_spiro['LogMNiErr+'], data_spiro['LogMNiErr-'])),
            elinewidth=0.8, label='Spiro et al. (2014)')
ax.errorbar(data_valenti['LogMNi'], data_valenti['Mv'], yerr=data_valenti['MvErr'], c='g', fmt='s', capsize=5,
            ms=7, capthick=0.5, xerr=np.vstack((data_valenti['LogMNiErr+'], data_valenti['LogMNiErr-'])),
            elinewidth=0.8, label='Valenti et al. (2015)')
ax.errorbar(data_others['LogMNi'], data_others['Mv'], yerr=data_others['MvErr'], c='k', fmt='o', capsize=5,
            ms=7, capthick=0.5, xerr=np.vstack((data_others['LogMNiErr+'], data_others['LogMNiErr-'])),
            elinewidth=0.8, label='Singh et al. (2018)')
ax.errorbar(np.log10(MNi), Mv50, xerr=MNiErr / MNi, yerr=Mv50Err, c='r', fmt='*', ms=18, capsize=7,
            elinewidth=1.5, label=name_SN)

set_plotparams(ax, markerfont=14, markerscale=1.4)
ax.set_ylim(-12.5, -19)
ax.set_xlim(-3.3, -0.3)
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))

corr = pearsonr(data_ni['LogMNi'], data_ni['Mv'])
ax.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_ni.shape[0], corr[0], corr[1]))
ax.set_ylabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}\ [mag]$', fontsize=18)
ax.set_xlabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=18)

fig.savefig('PLOT_Mv50+56Ni.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Mid-Plateau Velocity Vs Nickel Mass (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
data_V50 = data_sn[['LogV50', 'LogV50Err', 'LogMNi', 'LogMNiErr+', 'LogMNiErr-']].copy().dropna()
opt2, cov2 = curve_fit(line, data_V50['LogMNi'], data_V50['LogV50'], sigma=data_V50['LogV50Err'], p0=[2000, 9])
print (opt2, cov2, pearsonr(data_V50['LogMNi'], data_V50['LogV50']))

fig2 = plt.figure(figsize=(9, 9))
ax2 = fig2.add_subplot(111)

plot_confintervals(ax2, opt2, cov2, xlim=(-3.5, 0.0))
ax2.errorbar(data_hamuy['LogMNi'], data_hamuy['LogV50'], yerr=data_hamuy['LogV50Err'], c='darkorange', fmt='^',
             capsize=5, ms=7, xerr=np.vstack((data_hamuy['LogMNiErr+'], data_hamuy['LogMNiErr-'])),
             capthick=0.5, elinewidth=0.8, label='Hamuy (2003)')
ax2.errorbar(data_spiro['LogMNi'], data_spiro['LogV50'], yerr=data_spiro['LogV50Err'], c='b', fmt='<', capsize=5,
             ms=7, xerr=np.vstack((data_spiro['LogMNiErr+'], data_spiro['LogMNiErr-'])),
             capthick=0.5, elinewidth=0.8, label='Spiro et al. (2014)')
ax2.errorbar(data_valenti['LogMNi'], data_valenti['LogV50'], yerr=data_valenti['LogV50Err'], c='g', fmt='s',
             capsize=5, ms=7, xerr=np.vstack((data_valenti['LogMNiErr+'], data_valenti['LogMNiErr-'])),
             capthick=0.5, elinewidth=0.8, label='Valenti et al. (2015)')
ax2.errorbar(data_others['LogMNi'], data_others['LogV50'], yerr=data_others['LogV50Err'], c='k', fmt='o',
             capsize=5, ms=7, xerr=np.vstack((data_others['LogMNiErr+'], data_others['LogMNiErr-'])),
             capthick=0.5, elinewidth=0.8, label='Singh et al. (2018)')
ax2.errorbar(np.log10(MNi), np.log10(V50), xerr=MNiErr / MNi, yerr=V50Err / V50, c='r', fmt='*', ms=18,
             capsize=7, elinewidth=1.5, label=name_SN)

set_plotparams(ax2, markerfont=14, markerscale=1.4)
ax2.set_ylim(3.02, 3.95)
ax2.set_xlim(-3.3, -0.3)
ax2.yaxis.set_major_locator(MultipleLocator(0.2))
ax2.yaxis.set_minor_locator(MultipleLocator(0.02))
ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.xaxis.set_minor_locator(MultipleLocator(0.05))

corr2 = pearsonr(data_V50['LogMNi'], data_V50['LogV50'])
ax2.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_V50.shape[0], corr2[0], corr2[1]))
ax2.set_ylabel(r'Log [$\rm V_{50}/10^3\ km\ s^{-1}$]', fontsize=18)
ax2.set_xlabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=18)

fig2.savefig('PLOT_V50+56Ni.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Absolute V-Band Magnitude Vs Mid-Plateau Velocity (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
data_mv = data_sn[['LogV50', 'LogV50Err', 'Mv', 'MvErr']].copy().dropna()
opt3, cov3 = curve_fit(line, data_mv['Mv'], data_mv['LogV50'], sigma=data_mv['LogV50Err'], p0=[-0.16, -1.01])
print (opt3, cov3)

fig3 = plt.figure(figsize=(9, 9))
ax3 = fig3.add_subplot(111)

plot_confintervals(ax3, opt3, cov3, xlim=(-12.5, -19.0))
ax3.errorbar(data_hamuy['Mv'], data_hamuy['LogV50'], xerr=data_hamuy['MvErr'], yerr=data_hamuy['LogV50Err'], fmt='^',
             capsize=5, c='darkorange', ms=7, capthick=0.5, elinewidth=0.8, label='Hamuy (2003)')
ax3.errorbar(data_spiro['Mv'], data_spiro['LogV50'], xerr=data_spiro['MvErr'], yerr=data_spiro['LogV50Err'], fmt='<',
             capsize=5, c='b', ms=7, capthick=0.5, elinewidth=0.8, label='Spiro et al. (2014)')
ax3.errorbar(data_valenti['Mv'], data_valenti['LogV50'], xerr=data_valenti['MvErr'], yerr=data_valenti['LogV50Err'],
             fmt='s', capsize=5, c='g', ms=7, capthick=0.5, elinewidth=0.8, label='Valenti et al. (2015)')
ax3.errorbar(data_others['Mv'], data_others['LogV50'], xerr=data_others['MvErr'], yerr=data_others['LogV50Err'],
             fmt='o', c='k', capsize=5, ms=7, capthick=0.5, elinewidth=0.8, label='Singh et al. (2018)')
ax3.errorbar(Mv50, np.log10(V50), xerr=Mv50Err, yerr=V50Err / V50, c='r', fmt='*', ms=18, capsize=7,
             elinewidth=1.5, label=name_SN)

set_plotparams(ax3, markerfont=14, markerscale=1.4)
ax3.set_ylim(3.02, 3.95)
ax3.set_xlim(-19, -12.7)
ax3.yaxis.set_major_locator(MultipleLocator(0.2))
ax3.yaxis.set_minor_locator(MultipleLocator(0.02))
ax3.xaxis.set_major_locator(MultipleLocator(1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.1))

corr3 = pearsonr(data_mv['Mv'], data_mv['LogV50'])
ax3.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_mv.shape[0], corr3[0], corr3[1]))
ax3.set_ylabel(r'Log [$\rm V_{50}/10^3\ km\ s^{-1}$]', fontsize=18)
ax3.set_xlabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}\ [mag]$', fontsize=18)

fig3.savefig('PLOT_V50+Mv50.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig3)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Correlations - Absolute V-Band Magnitude Vs Nickel Mass & Mid-Plateau Velocity (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(18, 9))
ax1 = fig.add_subplot(121)

plot_confintervals(ax1, opt, cov, xlim=(-3.5, 0.0))
ax1.errorbar(data_hamuy['LogMNi'], data_hamuy['Mv'], yerr=data_hamuy['MvErr'], c='darkorange', fmt='^', capsize=5,
             ms=7, capthick=0.5, xerr=np.vstack((data_hamuy['LogMNiErr+'], data_hamuy['LogMNiErr-'])),
             elinewidth=0.8, label='Hamuy (2003)')
ax1.errorbar(data_spiro['LogMNi'], data_spiro['Mv'], yerr=data_spiro['MvErr'], c='b', fmt='<', capsize=5,
             ms=7, capthick=0.5, xerr=np.vstack((data_spiro['LogMNiErr+'], data_spiro['LogMNiErr-'])),
             elinewidth=0.8, label='Spiro et al. (2014)')
ax1.errorbar(data_valenti['LogMNi'], data_valenti['Mv'], yerr=data_valenti['MvErr'], c='g', fmt='s', capsize=5,
             ms=7, capthick=0.5, xerr=np.vstack((data_valenti['LogMNiErr+'], data_valenti['LogMNiErr-'])),
             elinewidth=0.8, label='Valenti et al. (2015)')
ax1.errorbar(data_others['LogMNi'], data_others['Mv'], yerr=data_others['MvErr'], c='k', fmt='o', capsize=5,
             ms=7, capthick=0.5, xerr=np.vstack((data_others['LogMNiErr+'], data_others['LogMNiErr-'])),
             elinewidth=0.8, label='Singh et al. (2018)')
ax1.errorbar(np.log10(MNi), Mv50, xerr=MNiErr / MNi, yerr=Mv50Err, c='r', fmt='*', ms=18, capsize=7,
             elinewidth=1.5, label=name_SN)

set_plotparams(ax1, markerfont=14, markerscale=1.6)
ax1.set_ylim(-12.7, -18.8)
ax1.set_xlim(-3.3, -0.3)
ax1.yaxis.set_major_locator(MultipleLocator(1))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.xaxis.set_major_locator(MultipleLocator(0.5))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))

ax1.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_ni.shape[0], corr[0], corr[1]))
ax1.set_ylabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}\ [mag]$', fontsize=18)
ax1.set_xlabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=18)

ax2 = fig.add_subplot(122)

optt, covt = curve_fit(line, data_mv['LogV50'], data_mv['Mv'], method='dogbox', p0=[-6, -6])
corrt = pearsonr(data_mv['LogV50'], data_mv['Mv'])
print optt, covt

plot_confintervals(ax2, optt, covt, xlim=(3.0, 4.0))
ax2.errorbar(data_hamuy['LogV50'], data_hamuy['Mv'], yerr=data_hamuy['MvErr'], xerr=data_hamuy['LogV50Err'], fmt='^',
             capsize=5, c='darkorange', ms=7, capthick=0.5, elinewidth=0.8, label='Hamuy (2003)')
ax2.errorbar(data_spiro['LogV50'], data_spiro['Mv'], yerr=data_spiro['MvErr'], xerr=data_spiro['LogV50Err'], fmt='<',
             capsize=5, c='b', ms=7, capthick=0.5, elinewidth=0.8, label='Spiro et al. (2014)')
ax2.errorbar(data_valenti['LogV50'], data_valenti['Mv'], yerr=data_valenti['MvErr'], xerr=data_valenti['LogV50Err'],
             fmt='s', capsize=5, c='g', ms=7, capthick=0.5, elinewidth=0.8, label='Valenti et al. (2015)')
ax2.errorbar(data_others['LogV50'], data_others['Mv'], yerr=data_others['MvErr'], xerr=data_others['LogV50Err'],
             fmt='o', c='k', capsize=5, ms=7, capthick=0.5, elinewidth=0.8, label='Singh et al. (2018)')
ax2.errorbar(np.log10(V50), Mv50, yerr=Mv50Err, xerr=V50Err / V50, c='r', fmt='*', ms=18, capsize=7,
             elinewidth=1.5, label=name_SN)

set_plotparams(ax2, legend=False)
ax2.set_xlim(3.02, 3.95)
ax2.set_ylim(-12.7, -18.8)
ax2.set_yticklabels([])
ax2.xaxis.set_major_locator(MultipleLocator(0.2))
ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))

ax2.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_mv.shape[0], corrt[0], corrt[1]))
ax2.set_xlabel(r'Log [$\rm V_{50}/10^3\ km\ s^{-1}$]', fontsize=18)

fig.subplots_adjust(wspace=0.01)
fig.savefig('PLOT_HamuyCorrelation.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Decline Rate (s2) Vs Peak Maximum Magnitude For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_s2max = data_and[['MMax', 'MMaxErr', 's2', 's2Err']].copy().dropna()
opt4, cov4 = curve_fit(line, data_s2max['MMax'], data_s2max['s2'], sigma=data_s2max['s2Err'], p0=[-0.5, -10])
print (opt4, cov4)

fig4 = plt.figure(figsize=(9, 9))
ax4 = fig4.add_subplot(111)

plot_confintervals(ax4, opt4, cov4, xlim=(-13.5, -18.5), fcolor='forestgreen')
ax4.errorbar(data_s2max['MMax'], data_s2max['s2'], xerr=data_s2max['MMaxErr'], yerr=data_s2max['s2Err'], fmt='D',
             c='k', capsize=3, ms=7, capthick=0.5, elinewidth=0.8, alpha=0.6, label='Anderson et al. (2014)')
ax4.errorbar(MMax, s2, yerr=s2Err, xerr=MMaxErr, c='r', fmt='*', ms=20, capsize=5, elinewidth=1.2, label=name_SN)
ax4.errorbar(MMax, s2 + 1, yerr=s2Err, xerr=MMaxErr, c='b', fmt='o', ms=14, capsize=5, markerfacecolor='grey',
             markeredgewidth=2, elinewidth=1.2, label=name_SN + r' [Without $\rm ^{56}Ni$]')

set_plotparams(ax4, markerfont=15, markerscale=1.2)
ax4.set_ylim(-0.8, 3.7)
ax4.set_xlim(-18.5, -13.5)
ax4.yaxis.set_major_locator(MultipleLocator(1.0))
ax4.yaxis.set_minor_locator(MultipleLocator(0.1))
ax4.xaxis.set_major_locator(MultipleLocator(1))
ax4.xaxis.set_minor_locator(MultipleLocator(0.1))

corr4 = pearsonr(data_s2max['MMax'], data_s2max['s2'])
ax4.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_s2max.shape[0], corr4[0], corr4[1]), fontsize=14)
ax4.set_ylabel(r's2 [$\rm mag\ (100\ d)^{-1}$]', fontsize=18)
ax4.set_xlabel(r'$\rm M^{max}_{V}\ [mag]$', fontsize=18)

fig4.savefig('PLOT_MMax+s2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig4)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Decline Rate (s1) Vs Plateau Length (OPTd) For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_opts1 = data_and[['OPTd', 'OPTdErr', 's1', 's1Err']].copy().dropna()
data_opts1 = data_opts1[data_opts1['s1'] < 6]
opt5, cov5 = curve_fit(line, data_opts1['OPTd'], data_opts1['s1'], sigma=data_opts1['s1Err'], p0=[-0.5, -10])
print (opt5, cov5)

fig5 = plt.figure(figsize=(9, 9))
ax5 = fig5.add_subplot(111)

plot_confintervals(ax5, opt5, cov5, xlim=(50, 130), fcolor='forestgreen')
ax5.errorbar(data_opts1['OPTd'], data_opts1['s1'], xerr=data_opts1['OPTdErr'], yerr=data_opts1['s1Err'], fmt='D',
             c='k', capsize=3, ms=7, capthick=0.5, elinewidth=0.8, alpha=0.6, label='Anderson et al. (2014)')
ax5.errorbar(OPTd, s1, yerr=s1Err, xerr=OPTdErr, c='r', fmt='*', ms=20, capsize=5, elinewidth=1.5, label=name_SN)
ax5.errorbar(OPTd, s1 + 1, yerr=s1Err, xerr=OPTdErr, c='b', fmt='o', ms=14, capsize=5, markerfacecolor='grey',
             markeredgewidth=2, elinewidth=1.2, label=name_SN + r' [Without $\rm ^{56}Ni$]')

set_plotparams(ax5, markerfont=15, markerscale=1.2)
ax5.set_ylim(0.4, 5.2)
ax5.set_xlim(50, 130)
ax5.yaxis.set_major_locator(MultipleLocator(1))
ax5.yaxis.set_minor_locator(MultipleLocator(0.1))
ax5.xaxis.set_major_locator(MultipleLocator(20))
ax5.xaxis.set_minor_locator(MultipleLocator(2))

corr5 = pearsonr(data_opts1['OPTd'], data_opts1['s1'])
ax5.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_opts1.shape[0], corr5[0], corr5[1]), fontsize=14)
ax5.set_ylabel(r's1 [$\rm mag\ (100\ d)^{-1}$]', fontsize=18)
ax5.set_xlabel(r'Plateau Length, OPTd [days]', fontsize=18)

fig5.savefig('PLOT_OPTd+s1.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig5)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Plateau Length (OPTd) Vs Peak Maximum Magnitude For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_optmmax = data_and[['OPTd', 'OPTdErr', 'MMax', 'MMaxErr']].copy().dropna()
opt6, cov6 = curve_fit(line, data_optmmax['OPTd'], data_optmmax['MMax'], sigma=data_optmmax['MMaxErr'], p0=[0.01, -20])
print (opt6, cov6)

fig6 = plt.figure(figsize=(9, 9))
ax6 = fig6.add_subplot(111)

plot_confintervals(ax6, opt6, cov6, xlim=(35, 130), fcolor='forestgreen')
ax6.errorbar(data_optmmax['OPTd'], data_optmmax['MMax'], xerr=data_optmmax['OPTdErr'], yerr=data_optmmax['MMaxErr'],
             c='k', fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=0.8, alpha=0.6, label='Anderson et al. (2014)')
ax6.errorbar(OPTd, MMax, xerr=OPTdErr, yerr=MMaxErr, c='r', fmt='*', ms=20, capsize=5, label=name_SN)

set_plotparams(ax6, markerfont=14, markerscale=1.2)
ax6.set_ylim(-18.5, -13.9)
ax6.set_xlim(35, 125)
ax6.yaxis.set_major_locator(MultipleLocator(1))
ax6.yaxis.set_minor_locator(MultipleLocator(0.1))
ax6.xaxis.set_major_locator(MultipleLocator(20))
ax6.xaxis.set_minor_locator(MultipleLocator(2))

corr6 = pearsonr(data_optmmax['OPTd'], data_optmmax['MMax'])
ax6.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_optmmax.shape[0], corr6[0], corr6[1]),
              fontsize=14)
ax6.set_ylabel(r'$\rm M^{max}_{V}\ [mag]$', fontsize=18)
ax6.set_xlabel(r'OPTd [days]', fontsize=18)

fig6.savefig('PLOT_OPTd+MMax.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig6)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Color (u-g)30 Vs MMax For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
optc1, covc1 = curve_fit(line, data_Mmaxug30['(u-g)30'], data_Mmaxug30['MMax'], sigma=data_Mmaxug30['MMaxErr'],
                         p0=[-2.5, -16])
print (optc1, covc1)

figc1 = plt.figure(figsize=(9, 9))
axc1 = figc1.add_subplot(111)

plot_confintervals(axc1, optc1, covc1, xlim=(0.5, 2.5), fcolor='forestgreen')
axc1.errorbar(data_Mmaxug30['(u-g)30'], data_Mmaxug30['MMax'], xerr=data_Mmaxug30['ColErr'],
              yerr=data_Mmaxug30['MMaxErr'], c='k', fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=1,
              alpha=0.6, label='De Jaeger et al. (2018a)')

axc1.errorbar(ug30, MMax, xerr=ug30Err, yerr=MMaxErr, c='r', fmt='*', ms=18, capsize=5, label=name_SN)

set_plotparams(axc1, markerscale=1.2, markerfont=14)
axc1.set_ylim(-14.7, -18.5)
axc1.set_xlim(0.5, 2.5)
axc1.yaxis.set_major_locator(MultipleLocator(1))
axc1.yaxis.set_minor_locator(MultipleLocator(0.1))
axc1.xaxis.set_major_locator(MultipleLocator(0.5))
axc1.xaxis.set_minor_locator(MultipleLocator(0.05))

corrc1 = pearsonr(data_Mmaxug30['(u-g)30'], data_Mmaxug30['MMax'])
axc1.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_Mmaxug30.shape[0], corrc1[0], corrc1[1]))
axc1.set_ylabel(r'$\rm M^{max}_{V}\ [mag]$', fontsize=18)
axc1.set_xlabel(r'$\rm (u-g)_{30\ d}$ [mag]', fontsize=18)

figc1.savefig('PLOT_MMax+u-g30.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(figc1)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Color (g-r)70 Vs s2 For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
optc2, covc2 = curve_fit(line, data_s2gr70['(g-r)70'], data_s2gr70['s2'], sigma=data_s2gr70['s2Err'],
                         p0=[2, -0.5])
print (optc2, covc2)

figc2 = plt.figure(figsize=(9, 9))
axc2 = figc2.add_subplot(111)

plot_confintervals(axc2, optc2, covc2, xlim=(0.1, 1.7), fcolor='forestgreen')
axc2.errorbar(data_s2gr70['(g-r)70'], data_s2gr70['s2'], xerr=data_s2gr70['ColErr'], yerr=data_s2gr70['s2Err'], c='k',
              fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=1, ls='',
              alpha=0.6, label='De Jaeger et al. (2018a)')

axc2.errorbar(gr70, s2, xerr=gr70Err, yerr=s2Err, c='r', fmt='*', ms=18, capsize=5, label=name_SN)

set_plotparams(axc2, markerscale=1.2, markerfont=14)
axc2.set_ylim(-0.25, 3.25)
axc2.set_xlim(0.2, 1.7)
axc2.yaxis.set_major_locator(MultipleLocator(0.5))
axc2.yaxis.set_minor_locator(MultipleLocator(0.1))
axc2.xaxis.set_major_locator(MultipleLocator(0.4))
axc2.xaxis.set_minor_locator(MultipleLocator(0.04))

corrc2 = pearsonr(data_s2gr70['(g-r)70'], data_s2gr70['s2'])
axc2.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_s2gr70.shape[0], corrc2[0], corrc2[1]))
axc2.set_ylabel(r'$\rm s_2\ [mag\ (100\ d)^{-1}]$', fontsize=18)
axc2.set_xlabel(r'$\rm (g-r)_{70\ d}$ [mag]', fontsize=18)

figc2.savefig('PLOT_s2+g-r70.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(figc2)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Color (g-r)70 Vs VHa70 For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
optc3, covc3 = curve_fit(line, data_VHagr70['(g-r)70'], data_VHagr70['VHa70'], sigma=data_VHagr70['VHa70Err'],
                         p0=[-40, 4])
print (optc3, covc3)

figc3 = plt.figure(figsize=(9, 9))
axc3 = figc3.add_subplot(111)

plot_confintervals(axc3, optc3, covc3, xlim=(0.2, 1.7), fcolor='forestgreen')
axc3.errorbar(data_VHagr70['(g-r)70'], data_VHagr70['VHa70'], xerr=data_VHagr70['ColErr'],
              yerr=data_VHagr70['VHa70Err'], c='k', fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=1, ls='',
              alpha=0.6, label='De Jaeger et al. (2018a)')
axc3.errorbar(gr70, VHa70 / 1000, xerr=gr70Err, yerr=VHa70Err / 1000, c='r', fmt='*', ms=18, capsize=5, label=name_SN)

set_plotparams(axc3, markerscale=1.2, markerfont=14)
axc3.set_ylim(3, 10)
axc3.set_xlim(0.2, 1.7)
axc3.yaxis.set_major_locator(MultipleLocator(1.0))
axc3.yaxis.set_minor_locator(MultipleLocator(0.1))
axc3.xaxis.set_major_locator(MultipleLocator(0.4))
axc3.xaxis.set_minor_locator(MultipleLocator(0.04))

corrc3 = pearsonr(data_VHagr70['(g-r)70'], data_VHagr70['VHa70'])
axc3.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_VHagr70.shape[0], corrc3[0], corrc3[1]))
axc3.set_ylabel(r'$\rm V^{H\alpha}_{70\ d}\ [\times\ 10^3\ km\ s^{-1}]$', fontsize=18)
axc3.set_xlabel(r'$\rm (g-r)_{70\ d}$ [mag]', fontsize=18)

figc3.savefig('PLOT_g-r70+VHa70.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(figc3)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Color (g-r)70 Vs VHa70,s2 For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
figcol = plt.figure(figsize=(9, 13))
axcol1 = figcol.add_subplot(211)
axcol2 = figcol.add_subplot(212, sharex=axcol1)

axcol1.errorbar(data_s2gr70['(g-r)70'], data_s2gr70['s2'], xerr=data_s2gr70['ColErr'], yerr=data_s2gr70['s2Err'],
                c='k', fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=1, ls='', alpha=0.6,
                label='De Jaeger et al. (2018a)')
axcol1.errorbar(gr70, s2, xerr=gr70Err, yerr=s2Err, c='r', fmt='*', ms=20, capsize=5, label=name_SN)

axcol2.errorbar(data_VHagr70['(g-r)70'], data_VHagr70['VHa70'], xerr=data_VHagr70['ColErr'],
                yerr=data_VHagr70['VHa70Err'], c='k', fmt='D', capsize=3, ms=7, capthick=0.5, elinewidth=1, ls='',
                alpha=0.6, label='_nolegend_')
axcol2.errorbar(gr70, VHa70 / 1000, xerr=gr70Err, yerr=VHa70Err / 1000, c='r', fmt='*', ms=20, capsize=5,
                label='_nolegend_')

axcol1.set_xlim(0.35, 1.7)
axcol1.set_ylim(-0.25, 3.3)
axcol2.set_ylim(3, 9.5)

plot_confintervals(axcol1, optc2, covc2, xlim=(0.1, 1.7), fcolor='dimgrey')
plot_confintervals(axcol2, optc3, covc3, xlim=(0.2, 1.7), fcolor='dimgrey')
set_plotparams(axcol1, markerscale=1.1, markerfont=14)
set_plotparams(axcol2, legend=False)

axcol1.yaxis.set_major_locator(MultipleLocator(1.0))
axcol1.yaxis.set_minor_locator(MultipleLocator(0.1))
axcol2.yaxis.set_major_locator(MultipleLocator(1.0))
axcol2.yaxis.set_minor_locator(MultipleLocator(0.1))
axcol1.xaxis.set_major_locator(MultipleLocator(0.2))
axcol1.xaxis.set_minor_locator(MultipleLocator(0.02))

axcol1.text(1.0, -0.15, r'$\rm N={0}, r = {1:0.2f}, p = {2:3.2e}$'.format(data_s2gr70.shape[0], corrc2[0], corrc2[1]),
            fontsize=13, color='brown')
axcol2.text(1.0, 3.2, r'$\rm N={0}, r = {1:0.2f}, p = {2:3.2e}$'.format(data_VHagr70.shape[0], corrc3[0], corrc3[1]),
            fontsize=13, color='brown')
axcol1.set_ylabel(r'$\rm s_2\ [mag\ (100\ d)^{-1}]$', fontsize=18)
axcol2.set_ylabel(r'$\rm V_{H\alpha}\ [\times\ 10^3\ km\ s^{-1}]$', fontsize=18)
axcol2.set_xlabel(r'$\rm (g-r)_{70\ d}$ [mag]', fontsize=18)

figcol.subplots_adjust(hspace=0.01)
figcol.savefig('PLOT_g-r70+VHa70|s2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(figcol)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Equivalent Width Evolution Of FeII 5169 Feature (Dessart et al. (2013))
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(111)

ax.plot(data_ewz2m3['Phase'], data_ewz2m3['EW5169'], c='k', ls=':', lw=1.5, marker='o', ms=10, alpha=0.6,
        markerfacecolor='darkorange', markeredgewidth=0.7, label=r'0.1 $Z_{\odot}$ [524 $\rm R_{\odot}$]')
ax.plot(data_ewz8m3['Phase'], data_ewz8m3['EW5169'], c='k', ls=':', lw=1.5, marker='D', ms=9, alpha=0.6,
        markerfacecolor='blue', markeredgewidth=0.7, label=r'0.4 $Z_{\odot}$ [611 $\rm R_{\odot}$]')
ax.plot(data_ewz4m2['Phase'], data_ewz4m2['EW5169'], c='k', ls=':', lw=1.5, marker='s', ms=9, alpha=0.6,
        markerfacecolor='red', markeredgewidth=0.7, label=r'2 $Z_{\odot}$ [804 $\rm R_{\odot}$]')
ax.plot(data_ewz2m2['Phase'], data_ewz2m2['EW5169'], c='k', ls=':', lw=1.5, marker='^', ms=10, alpha=0.6,
        markerfacecolor='green', markeredgewidth=0.7, label=r'1 $Z_{\odot}$ [768 $\rm R_{\odot}$]')
ax.plot(data_ewmlt1['Phase'], data_ewmlt1['EW5169'], c='brown', ls='-', lw=2, alpha=0.4,
        label=r'1 $Z_{\odot}$ [1107 $\rm R_{\odot}$]')
ax.plot(data_ewmlt3['Phase'], data_ewmlt3['EW5169'], c='brown', ls='--', lw=2, alpha=0.4,
        label=r'1 $Z_{\odot}$ [501 $\rm R_{\odot}$]')
ax.plot(data_ewsn['Phase'], data_ewsn['EW'], c='k', ls='', marker='*', ms=15, label=name_SN)

ax.set_ylim(-1, 60)
ax.set_xlim(12, 85)
ax.legend(markerscale=1.7, fontsize=15, frameon=False, handlelength=3)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)
ax.set_ylabel(r'pEW [FeII $\rm \lambda\ 5169$]', fontsize=18)
ax.set_xlabel(r'Time Since Explosion [Days]', fontsize=18)

fig.savefig('PLOT_pEWFeII5169.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plot Color (g-r)70 Vs s2 For Type II SNe
# # ------------------------------------------------------------------------------------------------------------------- #
# optc3, covc3 = curve_fit(line, data_s1gr15['(g-r)15'], data_s1gr15['s1'], sigma=data_s1gr15['s1Err'],
#                          p0=[-40, 4])
# print (optc3, covc3)

# figc3 = plt.figure(figsize=(8, 8))
# axc3 = figc3.add_subplot(111)

# plot_confintervals(axc3, optc3, covc3, xlim=(-0.15, 0.35), fcolor='chocolate,peru')
# axc3.errorbar(data_s1gr15['(g-r)15'], data_s1gr15['s1'], xerr=data_s1gr15['ColErr'], yerr=data_s1gr15['s1Err'],
#               c='k', fmt='D', capsize=3, ms=6, capthick=0.5, elinewidth=1, ls='',
#               alpha=0.6, label='De Jaeger et al. (2014)')

# axc3.errorbar(gr15, s1, xerr=gr15Err, yerr=s1Err, c='r', fmt='*', ms=13, capsize=5, label=name_SN)

# set_plotparams(axc3)
# axc3.set_ylim(0.5, 6.0)
# axc3.set_xlim(-0.2, 0.35)
# axc3.yaxis.set_major_locator(MultipleLocator(1.0))
# axc3.yaxis.set_minor_locator(MultipleLocator(0.1))
# axc3.xaxis.set_major_locator(MultipleLocator(0.1))
# axc3.xaxis.set_minor_locator(MultipleLocator(0.01))

# corrc3 = pearsonr(data_s1gr15['(g-r)15'], data_s1gr15['s1'])
# axc3.set_title(r'$\rm N={0}, r = {1:0.3f}, p = {2:4.3e}$'.format(data_s1gr15.shape[0], corrc3[0], corrc3[1]))
# axc3.set_ylabel(r'$\rm s_1\ [mag\ (100\ d)^{-1}]$', fontsize=16)
# axc3.set_xlabel(r'$\rm (g-r)15$ [mag]', fontsize=16)

# figc3.savefig('PLOT_s1+g-r15.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(figc3)
# # ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plot Oxygen Abundance Vs B-Band Absolute Magnitude
# # ------------------------------------------------------------------------------------------------------------------- #
# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111)

# ax.errorbar(hostMv, hostZ, xerr=hostMvErr, yerr=hostZErr, c='r', marker='*', ls='', capsize=5, ms=13, capthick=0.5,
#             elinewidth=1.5, label=name_hostgal)
# ax.scatter(data_gal['MB2'], data_gal['OH'], marker='o', s=1, c='darkgrey', label='_nolegend_')
# ax.scatter(data_sngal['BMag'], data_sngal['Z'], s=35, marker='^', c='navy', label='Type II SNe Hosts')

# set_plotparams(ax)
# ax.set_ylim(7.6, 9.5)
# ax.set_xlim(-11.5, -24.5)
# ax.yaxis.set_major_locator(MultipleLocator(0.4))
# ax.yaxis.set_minor_locator(MultipleLocator(0.04))
# ax.xaxis.set_major_locator(MultipleLocator(2))
# ax.xaxis.set_minor_locator(MultipleLocator(0.4))
# ax.set_ylabel(r'12 + log(O/H) [dex]', fontsize=16)
# ax.set_xlabel(r'$\rm M_B\ [mag]$', fontsize=16)

# fig.savefig('PLOT_GalZMetal.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #
