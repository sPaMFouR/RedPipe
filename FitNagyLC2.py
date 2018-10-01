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
from itertools import product
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
obsbol_file = 'OUTPUT_DateWiseSNBolFlux'
output_file ='OUTPUT_ChiSquare'

DIR_CURNT = os.getcwd()
DIR_Model = "/home/avinash/Dropbox/ModelSNe/Nagy_LC2/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
date_explosion = 2456841.50
light_speed = 299792.458

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
    ax_obj.legend()
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.set_xlabel('Time Since Explosion [Days]', fontsize=12)
    ax_obj.tick_params(which='both', direction='in', width=1, labelsize=12)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Generate Parameter Files With Specified Range To Be Used As Input To The Code
# ------------------------------------------------------------------------------------------------------------------- #

core_radii = 1e12 * np.arange(26, 32, 2)
core_mej = np.arange(8, 12, 1)
core_trec = np.arange(5000, 6000, 200)
core_mni = np.array([0.03])
core_eke = np.arange(1.2, 1.6, 0.1)
core_eth = np.arange(0.3, 0.5, 0.1)
core_exp = np.array([0.0])
core_pow = np.array([0.0])
core_kap = np.array([0.2])
core_mag = np.array([0.0])
core_tmag = np.array([0.0])
core_gamma = 1e5 * np.arange(1, 1000, 200)
core_epoch = np.array([100.0])

shell_radii = 1e12 * np.arange(75, 83, 2)
shell_mej = np.arange(0.5, 0.8, 0.1)
shell_trec = np.array([0.0])
shell_mni = np.array([0.0])
shell_eke = np.arange(1.3, 1.6, 0.1)
shell_eth = np.arange(0.07, 0.16, 0.03)
shell_exp = np.array([0.0])
shell_pow = np.array([2.0])
shell_kap = np.array([0.4])
shell_mag = np.array([0.0])
shell_tmag = np.array([0.0])
shell_gamma = 1e7 * np.arange(1, 1000, 200)
shell_epoch = np.array([100.0])

list_core = [core_radii, core_mej, core_trec, core_mni, core_eke, core_eth, core_exp, 
            core_pow, core_kap, core_mag, core_tmag, core_gamma, core_epoch]

list_shell = [shell_radii, shell_mej, shell_trec, shell_mni, shell_eke, shell_eth, shell_exp, 
              shell_pow, shell_kap, shell_mag, shell_tmag, shell_gamma, shell_epoch]

longlist_core = list(product(*list_core))
longlist_shell = list(product(*list_shell))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Reads I/O Data Files
# ------------------------------------------------------------------------------------------------------------------- #
obsbol_df = pd.read_csv(obsbol_file, sep='\s+', engine='python')
# core_df = pd.read_csv(core_file, sep='\s+', header=None, names=out_cols, comment='#', engine='python')
# shell_df = pd.read_csv(shell_file, sep='\s+', header=None, names=out_cols, comment='#', engine='python')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Perform Nagy Model Fit
# ------------------------------------------------------------------------------------------------------------------- #
print len(longlist_core)
print len(longlist_shell)

chisq_df = pd.DataFrame(columns=['CoreIndex', 'ShellIndex', 'ChiSquare'])

for core_index, core_values in enumerate(longlist_core):
    for shell_index, shell_values in enumerate(longlist_shell):
        core_srs = pd.Series(core_values)
        shell_srs = pd.Series(shell_values)
        core_srs.to_csv("par_core.inp", index=None)
        shell_srs.to_csv("par_shell.inp", index=None)
        
        os.system("./run ")
        
        comb_df = pd.read_csv(comb_file, sep='\s+', header=None, names=comb_cols, comment='#', engine='python')
        comb_df['Lum'] = comb_df['LumCore'] + comb_df['LumShell']
        
        chisq_df = chisq_df.append({'CoreIndex': core_index, 'ShellIndex': shell_index, 
                                    'ChiSquare': calc_chisquare(obsbol_df, comb_df)}, ignore_index=True)

min_index = chisq_df['ChiSquare'].idxmin()
min_row = chisq_df.loc[min_index, :]
print min_row['ChiSquare']
print longlist_core[int(min_row['CoreIndex'])]
print longlist_shell[int(min_row['ShellIndex'])]
                     
with open(output_file, 'w') as fout:
    fout.write("# Best-Fit Values Are: \n")
    fout.write('# Chi-Square : ' + str(min_row['ChiSquare']) + '\n')
    fout.write('# Core Index : ' + str(int(min_row['CoreIndex'])) + '\n')
    fout.write('# Shell Index : ' + str(int(min_row['ShellIndex'])) + '\n')
    
chisq_df.to_csv(output_file, index=None, header=True, sep=' ', mode='a')
# ------------------------------------------------------------------------------------------------------------------- #

