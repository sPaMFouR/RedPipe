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
# Global Variables & Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = 'ASASSN-14dq'
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
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

def plot_confintervals(ax_obj, optpar, covpar, xlim=(-3.5, 0.0)):
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

    ax_obj.plot(xdata, fit, linestyle='--', color='k', label='Our Fit')
    ax_obj.plot(xdata, fitlow, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xdata, fithigh, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xdata, fitlow2, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')

    ax_obj.fill_between(xdata, fitlow, fithigh, facecolor='grey', alpha=0.3)
    ax_obj.fill_between(xdata, fitlow, fitlow2, facecolor='lightgrey', alpha=0.3)
    ax_obj.fill_between(xdata, fithigh, fithigh2, facecolor='lightgrey', alpha=0.3)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Parameter Data For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_sn = pd.read_csv(DIR_SNe + 'Param_Data/Log_MNi.dat', sep='\s+', comment='#', engine='python')
data_sn = data_sn.replace('INDEF', np.nan).set_index('Name', drop=True).astype('float64')

data_sn['logMNi'] = data_sn['MNi'].apply(lambda x: np.log10(x))
data_sn['logMNiErr+'] = data_sn['MNiErr+'] / data_sn['MNi']
data_sn['logMNiErr-'] = data_sn['MNiErr-'] / data_sn['MNi']
data_sn['logv50'] = data_sn['v50'].apply(lambda x: np.log10(x))
data_sn['logv50Err'] = data_sn['v50Err'] / data_sn['v50']

data_spiro = data_sn.iloc[0:16].copy()
data_hamuy = data_sn.iloc[16:40].copy()
data_valenti = data_sn.iloc[40:48].copy()
data_others = data_sn.iloc[48:54].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Light Curve Parameters of Type II SNe (Anderson(2014))
# ------------------------------------------------------------------------------------------------------------------- #
data_and = pd.read_csv(DIR_SNe + 'Param_Data/Log_AndersonPar.dat', sep='\s+', comment='#', engine='python')
data_and = data_and.replace('INDEF', np.nan).set_index('Name', drop=True).astype('float64')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Galaxy Properties Data
# ------------------------------------------------------------------------------------------------------------------- #
data_gal = pd.read_csv(DIR_SNe + 'Param_Data/Log_TremontiGalZ.dat', sep='\s+', comment='#', engine='python')
data_gal = data_gal[['MB1', 'MB2', 'OH']].sort_values('MB2').replace('INDEF', np.nan).astype('float64').dropna()
data_gal = data_gal[(data_gal['MB2'] > -30) & (data_gal['MB2'] < -12)]

data_sngal = pd.read_csv(DIR_SNe + 'Param_Data/Log_SNZ.dat', sep='\s+', comment='#', engine='python')
data_sngal = data_sngal.replace('INDEF', np.nan).astype('float64')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Absolute V-Band Magnitude Vs Nickel Mass At Mid-Plateau
# ------------------------------------------------------------------------------------------------------------------- #
data_ni = data_sn[['Mv', 'MvErr', 'MNi', 'MNiErr+', 'MNiErr-', 'logMNi', 'logMNiErr+', 'logMNiErr-']].dropna()[:-1]
opt, cov = curve_fit(line, data_ni['logMNi'], data_ni['Mv'], sigma=data_ni['MvErr'], p0=[-2.5, -20.0])
print(opt, cov)

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

plot_confintervals(ax, opt, cov, xlim=(-3.5, 0.0))

ax.errorbar(np.log10(0.029), -16.9, xerr=0.006 / 0.029, yerr=0.2, color='r', fmt='*', markersize=15, capsize=6,
            label=name_SN)
ax.errorbar(data_hamuy['logMNi'], data_hamuy['Mv'], yerr=data_hamuy['MvErr'], color='orange', fmt='^', capsize=5,
            markersize=6, capthick=0.5, xerr=np.vstack((data_hamuy['logMNiErr+'], data_hamuy['logMNiErr-'])),
            elinewidth=1, label='Hamuy (2003)')
ax.errorbar(data_spiro['logMNi'], data_spiro['Mv'], yerr=data_spiro['MvErr'], color='b', fmt='<', capsize=5,
            markersize=6, capthick=0.5, xerr=np.vstack((data_spiro['logMNiErr+'], data_spiro['logMNiErr-'])),
            elinewidth=1, label='Spiro et al. (2014)')
ax.errorbar(data_valenti['logMNi'], data_valenti['Mv'], yerr=data_valenti['MvErr'], color='g', fmt='s', capsize=5,
            markersize=6, capthick=0.5, xerr=np.vstack((data_valenti['logMNiErr+'], data_valenti['logMNiErr-'])),
            elinewidth=1, label='Valenti et al. (2015)')
ax.errorbar(data_others['logMNi'], data_others['Mv'], yerr=data_others['MvErr'], color='k', fmt='o', capsize=5,
            markersize=6, capthick=0.5, xerr=np.vstack((data_others['logMNiErr+'], data_others['logMNiErr-'])),
            elinewidth=1, label='Table 5')

ax.set_ylim(-12.5, -19)
ax.set_xlim(-3.3, -0.3)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))
ax.xaxis.set_major_locator(MultipleLocator(0.5))
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.tick_params(which='both', direction='in', width=1, labelsize=14)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels, fontsize=12, markerscale=1.5, frameon=False)

corr = pearsonr(data_ni['logMNi'], data_ni['Mv'])
ax.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_ni.shape[0], corr[0], corr[1]))
ax.set_ylabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}\ [mag]$', fontsize=16)
ax.set_xlabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=16)

fig.savefig('OUTPUT_PlotNi.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Mid-Plateau Velocity Vs Nickel Mass (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
data_v50 = data_sn[['logv50', 'logv50Err', 'MNi', 'MNiErr+', 'MNiErr-', 'logMNi', 'logMNiErr+', 'logMNiErr-']].dropna()
opt2, cov2 = curve_fit(line, data_v50['logMNi'], data_v50['logv50'], sigma=data_v50['logv50Err'], p0=[2000, 9])
print(opt2, cov2, pearsonr(data_v50['logMNi'], data_v50['logv50']))

fig2 = plt.figure(figsize=(8, 8))
ax2 = fig2.add_subplot(111)

plot_confintervals(ax2, opt2, cov2, xlim=(-3.5, 0.0))

ax2.errorbar(np.log10(0.029), np.log10(4910), xerr=0.005 / 0.029, yerr=100 / 4910, color='r', fmt='*', markersize=15,
             capsize=3, label=name_SN)
ax2.errorbar(data_hamuy['logMNi'], data_hamuy['logv50'], yerr=data_hamuy['logv50Err'], color='orange', fmt='^',
             capsize=5, markersize=6, xerr=np.vstack((data_hamuy['logMNiErr+'], data_hamuy['logMNiErr-'])),
             capthick=0.5, elinewidth=1, label='Hamuy (2003)')
ax2.errorbar(data_spiro['logMNi'], data_spiro['logv50'], yerr=data_spiro['logv50Err'], color='b', fmt='<', capsize=5,
             markersize=6, xerr=np.vstack((data_spiro['logMNiErr+'], data_spiro['logMNiErr-'])),
             capthick=0.5, elinewidth=1, label='Spiro et al. (2014)')
ax2.errorbar(data_valenti['logMNi'], data_valenti['logv50'], yerr=data_valenti['logv50Err'], color='g', fmt='s',
             capsize=5, markersize=6, xerr=np.vstack((data_valenti['logMNiErr+'], data_valenti['logMNiErr-'])),
             capthick=0.5, elinewidth=1, label='Valenti et al. (2015)')
ax2.errorbar(data_others['logMNi'], data_others['logv50'], yerr=data_others['logv50Err'], color='k', fmt='o',
             capsize=5, markersize=6, xerr=np.vstack((data_others['logMNiErr+'], data_others['logMNiErr-'])),
             capthick=0.5, elinewidth=1, label='Table 5')

ax2.set_ylim(3.02, 3.95)
ax2.set_xlim(-3.3, -0.3)
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.xaxis.set_major_locator(MultipleLocator(0.5))
ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
ax2.tick_params(which='both', direction='in', width=1, labelsize=14)

handles2, labels2 = ax2.get_legend_handles_labels()
handles2 = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles2]
ax2.legend(handles2, labels2, fontsize=12, markerscale=1.5, frameon=False)

corr2 = pearsonr(data_v50['logMNi'], data_v50['logv50'])
ax2.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_v50.shape[0], corr2[0], corr2[1]))
ax2.set_ylabel(r'Log [$\rm V_{50}/10^3\ km\ s^{-1}$]', fontsize=16)
ax2.set_xlabel(r'$\rm Log\ [M_{Ni}/ M_{\odot}]$', fontsize=16)

fig2.savefig('OUTPUT_PlotV50.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Absolute V-Band Magnitude Vs Mid-Plateau Velocity (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
data_mv = data_sn[['v50', 'v50Err', 'logv50', 'logv50Err', 'Mv', 'MvErr']].dropna()
opt3, cov3 = curve_fit(line, data_mv['Mv'], data_mv['logv50'], sigma=data_mv['logv50Err'], p0=[-0.16, -1.01])
print(opt3, cov3)

fig3 = plt.figure(figsize=(8, 8))
ax3 = fig3.add_subplot(111)

plot_confintervals(ax3, opt3, cov3, xlim=(-12.5, -19.0))

ax3.errorbar(-16.9, np.log10(4910), yerr=100 / 4910, xerr=0.2, color='r', fmt='*', markersize=15, capsize=3,
             label=name_SN)
ax3.errorbar(data_hamuy['Mv'], data_hamuy['logv50'], xerr=data_hamuy['MvErr'], color='orange', fmt='^', capsize=5,
             markersize=6, capthick=0.5, yerr=data_hamuy['logv50Err'], elinewidth=1, label='Hamuy (2003)')
ax3.errorbar(data_spiro['Mv'], data_spiro['logv50'], xerr=data_spiro['MvErr'], color='b', fmt='<', capsize=5,
             markersize=6, capthick=0.5, yerr=data_spiro['logv50Err'], elinewidth=1, label='Spiro et al. (2014)')
ax3.errorbar(data_valenti['Mv'], data_valenti['logv50'], xerr=data_valenti['MvErr'], color='g', fmt='s', capsize=5,
             markersize=6, capthick=0.5, yerr=data_valenti['logv50Err'], elinewidth=1, label='Valenti et al. (2015)')
ax3.errorbar(data_others['Mv'], data_others['logv50'], xerr=data_others['MvErr'], color='k', fmt='o', capsize=5,
             markersize=6, capthick=0.5, yerr=data_others['logv50Err'], elinewidth=1, label='Table 5')

handles3, labels3 = ax3.get_legend_handles_labels()
handles3 = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles3]
ax3.legend(handles3, labels3, fontsize=12, markerscale=1.5, frameon=False)

ax3.set_ylim(3.02, 3.95)
ax3.set_xlim(-19, -12.7)
ax3.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
ax3.yaxis.set_major_locator(MultipleLocator(0.1))
ax3.yaxis.set_minor_locator(MultipleLocator(0.05))
ax3.xaxis.set_major_locator(MultipleLocator(1))
ax3.xaxis.set_minor_locator(MultipleLocator(0.25))
ax3.tick_params(which='both', direction='in', width=1, labelsize=14)

corr3 = pearsonr(data_mv['Mv'], data_mv['logv50'])
ax3.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_mv.shape[0], corr3[0], corr3[1]))
ax3.set_ylabel(r'Log [$\rm V_{50}/10^3\ km\ s^{-1}$]', fontsize=16)
ax3.set_xlabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}\ [mag]$', fontsize=16)

fig3.savefig('OUTPUT_PlotMv.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig3)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Decline Rate (s2) Vs Peak Maximum Magnitude For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_s2max = data_and[['MMax', 'MMaxErr', 's2', 's2Err']].dropna()
opt4, cov4 = curve_fit(line, data_s2max['MMax'], data_s2max['s2'], sigma=data_s2max['s2Err'], p0=[-0.5, -10])
print(opt4, cov4)

fig4 = plt.figure(figsize=(8, 8))
ax4 = fig4.add_subplot(111)

plot_confintervals(ax4, opt4, cov4, xlim=(-13.5, -18.5))

ax4.errorbar(-16.9, 1.18, yerr=0.05, xerr=0.2, color='r', fmt='*', markersize=15, capsize=3,
             label=name_SN)
ax4.errorbar(data_s2max['MMax'], data_s2max['s2'], xerr=data_s2max['MMaxErr'], color='darkorange', fmt='^', capsize=5,
             markersize=6, capthick=0.5, yerr=data_s2max['s2Err'], elinewidth=1, label='Anderson et al. (2014)')

handles4, labels4 = ax4.get_legend_handles_labels()
handles4 = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles4]
ax4.legend(handles4, labels4, fontsize=12, markerscale=1.5, frameon=False)

# ax4.set_ylim(3.02, 3.95)
ax4.set_xlim(-18.5, -13.5)
ax4.yaxis.set_ticks_position('both')
ax4.xaxis.set_ticks_position('both')
ax4.yaxis.set_major_locator(MultipleLocator(0.5))
ax4.yaxis.set_minor_locator(MultipleLocator(0.125))
ax4.xaxis.set_major_locator(MultipleLocator(1))
ax4.xaxis.set_minor_locator(MultipleLocator(0.25))
ax4.tick_params(which='both', direction='in', width=1, labelsize=14)

corr4 = pearsonr(data_s2max['MMax'], data_s2max['s2'])
ax4.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_s2max.shape[0], corr4[0], corr4[1]))
ax4.set_ylabel(r's2 [$\rm mag\ /\ 100\ d^{-1}$]', fontsize=16)
ax4.set_xlabel(r'$\rm M^{max}_{V}\ [mag]$', fontsize=16)

fig4.savefig('OUTPUT_PlotMMaxs2.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig4)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Decline Rate (s2) Vs Peak Maximum Magnitude For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data_opts1 = data_and[['OPTd', 'OPTdErr', 's1', 's1Err']].dropna()
data_opts1 = data_opts1[data_opts1['s1'] < 6]
opt5, cov5 = curve_fit(line, data_opts1['OPTd'], data_opts1['s1'], sigma=data_opts1['s1Err'], p0=[-0.5, -10])
print(opt5, cov5)

fig5 = plt.figure(figsize=(8, 8))
ax5 = fig5.add_subplot(111)

plot_confintervals(ax5, opt5, cov5, xlim=(25, 125))

ax5.errorbar(90, 1.80, yerr=0.05, xerr=5, color='r', fmt='*', markersize=15, capsize=3, label=name_SN)
ax5.errorbar(data_opts1['OPTd'], data_opts1['s1'], xerr=data_opts1['OPTdErr'], color='darkorange', fmt='^', capsize=5,
             markersize=6, capthick=0.5, yerr=data_opts1['s1Err'], elinewidth=1, label='Anderson et al. (2014)')

handles5, labels5 = ax5.get_legend_handles_labels()
handles5 = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles5]
ax5.legend(handles5, labels5, fontsize=12, markerscale=1.5, frameon=False)

ax5.set_ylim(0.5, 6)
ax5.set_xlim(35, 130)
ax5.yaxis.set_ticks_position('both')
ax5.xaxis.set_ticks_position('both')
ax5.yaxis.set_major_locator(MultipleLocator(2))
ax5.yaxis.set_minor_locator(MultipleLocator(0.5))
ax5.xaxis.set_major_locator(MultipleLocator(40))
ax5.xaxis.set_minor_locator(MultipleLocator(10))
ax5.tick_params(which='both', direction='in', width=1, labelsize=14)

corr5 = pearsonr(data_opts1['OPTd'], data_opts1['s1'])
ax5.set_title(r'$\rm N={0}, r = {1:0.4f}, p = {2:5.4e}$'.format(data_opts1.shape[0], corr5[0], corr5[1]))
ax5.set_ylabel(r's1 [$\rm mag\ /\ 100\ d^{-1}$]', fontsize=16)
ax5.set_xlabel(r'OPTd [days]', fontsize=16)

fig5.savefig('OUTPUT_PlotOPTds1.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig5)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Oxygen Abundance Vs B-Band Absolute Magnitude
# ------------------------------------------------------------------------------------------------------------------- #
fig6 = plt.figure(figsize=(8, 10))
ax6 = fig6.add_subplot(111)

xmag = np.linspace(-12, -24, 1000)
ax6.errorbar(-17.5, 8.40, xerr=0.52, yerr=0.18, color='red', marker='*', linestyle='', capsize=5, ms=12, capthick=0.5,
             label='UGC 11860')
ax6.scatter(data_gal['MB2'], data_gal['OH'], marker='o', s=6, color='darkgrey', label='_nolegend_')
ax6.scatter(data_sngal['II_BMag'], data_sngal['II_Z'], marker='^', color='b', label='Type II SNe Hosts')
ax6.plot(xmag, line(xmag, -0.185, 5.238), linestyle='--', color='k', label='Tremonti et al. (2004)')
ax6.plot(xmag, line(xmag, -0.11, 6.27), linestyle='-.', color='g', label='Berg et al. (2012)')
# ax6.plot(xmag, line(xmag, -0.223, 4.07), linestyle='-', color='y', label='Lamareille et al. (2004) Fit')

handles6, labels6 = ax6.get_legend_handles_labels()
handles6 = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles6]
ax6.legend(handles6, labels6, fontsize=12, markerscale=1.5, frameon=False)

ax6.set_ylim(7.3, 9.5)
ax6.set_xlim(-11.5, -24.5)
ax6.yaxis.set_ticks_position('both')
ax6.xaxis.set_ticks_position('both')
ax6.yaxis.set_major_locator(MultipleLocator(0.5))
ax6.yaxis.set_minor_locator(MultipleLocator(0.1))
ax6.xaxis.set_major_locator(MultipleLocator(2))
ax6.xaxis.set_minor_locator(MultipleLocator(0.5))
ax6.tick_params(which='both', direction='in', width=1.5, labelsize=14)

ax6.set_ylabel(r'12 + log(O/H) [dex]', fontsize=16)
ax6.set_xlabel(r'$\rm M_B\ [mag]$', fontsize=16)

fig6.savefig('OUTPUT_PlotGalMetal.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig6)
# ------------------------------------------------------------------------------------------------------------------- #
