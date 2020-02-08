#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx---------------DETERMINE THE EMPIRICAL RELATION FOR TYPE IIP SUPERNOVA---------------xxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from scipy.stats import linregress, pearsonr, spearmanr
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
file_name = 'TypeIISNe.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def line(x, m, c):
    return m * x + c


def elmfit(x):
    return -6.2295 * x - 0.8147


def extractnum(pattern):
    number = re.findall('\d+', pattern)[0]
    return len(number) - 2


def get_suffix(name):
    name = name[-5:]
    if name[0] == '-':
        return name[1:]
    else:
        return name[extractnum(name):]


def calc_sigma(x, covmc):
    return np.sqrt(covmc[0, 0] * (x ** 2) + covmc[1, 1] + 2 * x * covmc[0, 1])


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

    fitlow2 = fit - 3 * sigma
    fithigh2 = fit + 3 * sigma

    ax_obj.plot(xdata, fit, linestyle='--', color='k', label='Our Fit')
    ax_obj.plot(xdata, fitlow2, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')
    ax_obj.plot(xdata, fithigh2, linestyle='-.', color='k', linewidth=0.5, alpha=0.5, label='_nolegend_')

    ax_obj.fill_between(xdata, fitlow2, fithigh2, facecolor='khaki', alpha=0.3)


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Sample Of Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
data = pd.read_csv(file_name, sep='\s+', comment='#')
data['Label'] = data['Name'].apply(get_suffix)
data = data.replace('INDEF', np.nan).set_index(['Label']).drop(['Name', 'Marker', 'Color'], axis=1).astype('float64')

data['LogMNi'] = data['MNi'].apply(lambda x: np.log10(x))
data['LogMNiErr+'] = data['MNiErr+'] / data['MNi']
data['LogMNiErr-'] = data['MNiErr-'] / data['MNi']
data['LogMNi2'] = data['MNi2'].apply(lambda x: np.log10(x))
data['LogMNiErr2'] = data['MNiErr2'] / data['MNi2']
data = data[['S', 'LogMNi', 'LogMNiErr+', 'LogMNiErr-', 'LogMNi2', 'LogMNiErr2']].dropna()
data['LabelX'] = data['S'] + 0.005
data['LabelY'] = data['LogMNi']

data.loc['99gi', 'LabelY'] -= 0.06
data.loc['99gi', 'LabelX'] -= 0.02
data.loc['13K', 'LabelY'] += 0.03
data.loc['13K', 'LabelX'] -= 0.025
data.loc['13ej', 'LabelY'] -= 0.06
data.loc['13ej', 'LabelX'] -= 0.02
data.loc['04ej', 'LabelY'] -= 0.04
data.loc['04fx', 'LabelY'] -= 0.02
data.loc['12aw', 'LabelX'] -= 0.05
data.loc['12aw', 'LabelY'] -= 0.02
data.loc['88A', 'LabelX'] -= 0.04
data.loc['69L', 'LabelY'] += 0.03
data.loc['69L', 'LabelX'] -= 0.02
data.loc['87A', 'LabelY'] += 0.03
data.loc['87A', 'LabelX'] -= 0.01
data.loc['07it', 'LabelY'] -= 0.07
data.loc['07it', 'LabelX'] -= 0.02
data.loc['04et', 'LabelY'] -= 0.07
data.loc['04et', 'LabelX'] -= 0.03
data.loc['02hx', 'LabelY'] += 0.02
data.loc['17eaw', 'LabelY'] -= 0.04
data.loc['99em', 'LabelY'] -= 0.02
data.loc['09ib', 'LabelY'] -= 0.02
data.loc['70G', 'LabelY'] -= 0.06
data.loc['70G', 'LabelX'] -= 0.02
data.loc['13dpa', 'LabelY'] -= 0.02
data.loc['08in', 'LabelY'] -= 0.05
data.loc['09md', 'LabelX'] -= 0.02
data.loc['09md', 'LabelY'] -= 0.06
data.loc['13ab', 'LabelY'] -= 0.02
data.loc['05af', 'LabelX'] -= 0.04

# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Performing Linear Fit To The Steepness Parameter and Nickel Mass
# ------------------------------------------------------------------------------------------------------------------- #
xdata = np.linspace(0.0, 0.55, 1000)
opt, cov = curve_fit(line, data['S'], data['LogMNi'], p0=[-3, -1])
a = linregress(data['S'], data['LogMNi'])

print(a, cov)
print(pearsonr(data['S'], data['LogMNi']))
print(spearmanr(data['S'], data['LogMNi']))

print(10 ** (-6.2995 * 0.131 - 0.8147))
print(10 ** (a[0] * 0.131 + a[1]))
print(10 ** (a[0] * 0.121 + a[1]) - 10 ** (a[0] * 0.131 + a[1]))
# data = data.drop('14dq', axis=0)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Nickel Mass Vs Steepness Parameter For Type II SNe
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

plot_confintervals(ax, opt, cov, xlim=(-0.02, 0.57))
ax.scatter(0.131, np.log10(0.029), marker='*', s=100, color='r', label='ASASSN-14dq')
ax.plot(xdata, elmfit(xdata), linestyle='-', color='saddlebrown', label='Elmhamdi et al. (2003) Fit')
ax.scatter(data['S'][0:9], data['LogMNi'][0:9], s=40, marker='x', color='chocolate',
           label='Elmhamdi et al. (2003) Sample')
ax.scatter(data['S'][9:], data['LogMNi'][9:], s=25, marker='^', color='b', label='Our Sample')

for index, row in data.iterrows():
    ax.text(row['LabelX'], row['LabelY'], index, fontsize=10)

ax.set_ylim(-3.0, -0.8)
ax.set_xlim(-0.02, 0.57)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.125))
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.025))
ax.legend(markerscale=2, fontsize=12, loc=3, frameon=False)
ax.tick_params(which='both', direction='in', width=1.5, labelsize=14)

pcorr = pearsonr(data['S'], data['LogMNi'])
scorr = spearmanr(data['S'], data['LogMNi'])

ax.text(0.33, -1.0, 'N=39', fontsize=12)
ax.text(0.25, -1.2, r'$\rm r_p={0:0.4f}, p_p\leq {1:4.2e}$'.format(pcorr[0], pcorr[1]) + '\n' +
        r'$\rm r_s={0:0.4f}, p_s\leq {1:4.2e}$'.format(scorr[0], scorr[1], fontsize=12))
ax.set_ylabel(r'Log M($^{56}$Ni)', fontsize=16)
ax.set_xlabel(r'Steepness Parameter [$\rm mag\ day^{-1}$]', fontsize=16)

fig.savefig('OUTPUT_PlotNiSteepness.eps', dpi=1000, format='eps', bbox_inches='tight')
plt.show()

# ------------------------------------------------------------------------------------------------------------------- #
