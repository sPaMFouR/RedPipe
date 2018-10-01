#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxx-------------------Nagy LIGHT CURVE MODEL-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import InterpolatedUnivariateSpline
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Files And Paths of Directories
# ------------------------------------------------------------------------------------------------------------------- #
core_file = 'core.out'
shell_file = 'shell.out'
comb_file = 'comb.out'
comb_file2 = 'comb2.out'
obsbol_file = 'OUTPUT_DateWiseSNBolFlux'

DIR_CURNT = os.getcwd()
DIR_Model = "/home/avinash/Dropbox/ModelSNe/Nagy_LC2/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
date_explosion = 2456841.50

out_cols = ['Phase', 'Lum', 'LogLum', 'LDiff', 'LRec', 'Xi']
comb_cols = ['Phase', 'LumCore', 'LumShell']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Calculate Chi-Square Statistic
# ------------------------------------------------------------------------------------------------------------------- #
def calc_chisquare(obs_data, model_data):
    """
    Calculates the goodness-of-fit through the chi-square parameter. Returns the chi-square statistic.
    Args:
        obs_data             : Observed Data
        model_data           : Model Data
    Returns:
        fit_chisq.statistic  : The sum of squared-errors between the model and the observed data
    """
    fit_interp = InterpolatedUnivariateSpline(model_data['Phase'], model_data['Lum'], k=1)
    fit_chisq = chisquare(obs_data['Lum'], fit_interp(obs_data['Phase']))

    return fit_chisq.statistic

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Plotting
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj):
    """
    Sets the plot parameters for displaying the Light Curves and Expansion Velocity Evolution.
    Args:
        ax_obj   : Axes object to be used for plotting and setting plot parameters
    Returns:
        None
    """
    ax_obj.set_ylim(4e40, 5e42)
    ax_obj.set_xlim(-10, 210)
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.set_xlabel('Time Since Explosion [Days]', fontsize=14)
    ax_obj.tick_params(which='both', direction='in', width=1, labelsize=14)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Reads I/O Data Files
# ------------------------------------------------------------------------------------------------------------------- #
# os.system('./run')

obsbol_df = pd.read_csv(obsbol_file, sep='\s+', engine='python')
obsbol_df = obsbol_df[obsbol_df['Phase'] < 200]
obsbol_df = obsbol_df[(obsbol_df['Date'] != '2014-11-08') & (obsbol_df['Date'] != '2014-11-14')]
obsbol_df = obsbol_df[(obsbol_df['Date'] != '2014-09-29') & (obsbol_df['Date'] != '2014-12-02')]

comb_df = pd.read_csv(comb_file, sep='\s+', header=None, names=comb_cols, comment='#', engine='python')
comb_df['Lum'] = comb_df['LumCore'] + comb_df['LumShell']

comb2_df = pd.read_csv(comb_file2, sep='\s+', header=None, names=comb_cols, comment='#', engine='python')
comb2_df['Lum'] = comb2_df['LumCore'] + comb2_df['LumShell']

print calc_chisquare(obsbol_df, comb_df)
print calc_chisquare(obsbol_df, comb2_df)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 6))

ax1 = fig.add_subplot(121)
ax1.semilogy(obsbol_df['Phase'], obsbol_df['Lum'], c='b', marker='o', linestyle='', label='Observed LC')
ax1.semilogy(comb_df['Phase'], comb_df['Lum'], c='k', linestyle='-', linewidth=0.7, label='Model Core (a=0)')
ax1.semilogy(comb_df['Phase'], comb_df['LumCore'], c='r', linestyle='--', linewidth=0.7, label='Model Shell (n=2)')
ax1.semilogy(comb_df['Phase'], comb_df['LumShell'], c='g', linestyle='--', linewidth=0.7, label='Model Combined LC')

set_plotparams(ax1)
ax1.set_ylabel(r'Bolometric Luminosity [$\rm erg\ s^{-1}$]', fontsize=14)
ax1.text(120, 6e41, s=r'$\rm \kappa_{shell}\ =\ 0.4\ cm^2\ g^{-1}$', fontsize=12)

ax2 = fig.add_subplot(122, sharey=ax1)
ax2.semilogy(obsbol_df['Phase'], obsbol_df['Lum'], c='b', marker='o', linestyle='', label='Observed LC')
ax2.semilogy(comb2_df['Phase'], comb2_df['Lum'], c='k', linestyle='-', linewidth=0.7, label='Model Core (a=0)')
ax2.semilogy(comb2_df['Phase'], comb2_df['LumCore'], c='r', linestyle='--', linewidth=0.7, label='Model Shell (n=2)')
ax2.semilogy(comb2_df['Phase'], comb2_df['LumShell'], c='g', linestyle='--', linewidth=0.7, label='Model Combined LC')

set_plotparams(ax2)
ax2.legend(fontsize=12, frameon=False)
ax2.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax2.text(120, 6e41, s=r'$\rm \kappa_{shell}\ =\ 0.3\ cm^2\ g^{-1}$', fontsize=12)

fig.subplots_adjust(wspace=0.001, top=0.9, right=0.95)
fig.savefig('OUTPUT_FitNagy.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)

