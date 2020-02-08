#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxx---------------Plot The Supernova Photospheric Velocity Evolution----------------xxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
from jdcal import gcal2jd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
date_explosion = 2456841.50
light_speed = 2.99792458e5
dict_snmark = {'1999em': ['o', 'c'], '2012aw': ['D', 'b'], '2013ej': ['s', 'g'], '1999gi': ['^', 'r'], 
               '2004et': ['p', 'm'], '2005cs': ['*', 'y'], '2013ab': ['8', 'coral'], '2009bw': ['h', 'violet'], 
               '2006bp' : ['+', 'lime']}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
DIR_PHOT = "/home/avinash/Supernovae_Data/Photometry/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def group_similar_files(text_list, common_text, exceptions=''):
    """
    Groups similar files based on the string 'common_text'. Writes the similar files
    onto the list 'text_list' (only if this string is not empty) and appends the similar
    files to a list 'python_list'.
    Args:
        text_list   : Name of the output text file with names grouped based on the 'common_text'
        common_text : String containing partial name of the files to be grouped
        exceptions  : String containing the partial name of the files that need to be excluded
    Returns:
        list_files  : Python list containing the names of the grouped files
    """
    list_files = glob.glob(common_text)
    if exceptions != '':
        list_exception = exceptions.split(',')
        for file_name in glob.glob(common_text):
            for text in list_exception:
                test = re.search(text, file_name)
                if test:
                    try:
                        list_files.remove(file_name)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name + '\n')

    return list_files

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Converting Calendar Date To Julian Day
# ------------------------------------------------------------------------------------------------------------------- #

def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split('-')
    jd_tuple = gcal2jd(date_comp[0], date_comp[1], str(int(float(date_comp[2])) + 1))
    julian_day = jd_tuple[0] + jd_tuple[1]

    return julian_day

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj, ylim=13.5):
    """
    This function sets the plot parameters for Photospheric Velocity evolution plots using "ax_obj" of the plot.
    Args:
        ax_obj   : Axes object to which the plot parameters are to be applied
        ylim     : Upper limit of Y-axis
    Returns:
        None
    """
    ax_obj.set_xlim(-2, 190)
    ax_obj.set_ylim(0.5, ylim)
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(40))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.yaxis.set_major_locator(MultipleLocator(2))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax_obj.tick_params(which='both', direction='in', width=1.5, labelsize=14)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Pandas DataFrame With Velocity Information From Different Spectral Lines
# ------------------------------------------------------------------------------------------------------------------- #
wave_df = pd.read_csv('OUTPUT_PhotVel', sep='\s+', engine='python')
wave_df = wave_df.replace('INDEF', np.nan, regex=True)
wave_df['JD'] = wave_df['Date'].apply(lambda x: cald_to_jd(x))
wave_df['Phase'] = wave_df['JD'] - date_explosion
wave_df = wave_df.set_index('Phase', drop=True)

vel_df = pd.DataFrame()
vel_df['Date'] = wave_df['Date'].copy()

dict_markers = {6563: ['o', 'k', r'H$\alpha$'], 4861: ['*', 'teal', r'H$\beta$'], 5169: ['s', 'g', r'$\rm Fe\,II$ 5169'], 
                5018: ['^', 'b', r'$\rm Fe\,II\ 5018$'], 4924: ['p', 'orange', r'$\rm Fe\,II$ 4924'], 
                4340: ['8', 'r', r'H$\gamma$']}

for wave in dict_markers.keys():
    vel_df[str(wave)] = wave_df[str(wave)].apply(lambda x: (wave - float(x)) * light_speed / wave).round(0)
    vel_df[str(wave) + 'Err2'] = wave_df[str(wave) + 'Err'].apply(lambda x: float(x) * light_speed / wave).round(0)
    vel_df[str(wave) + 'Err'] = wave_df[str(wave) + 'Err'].apply(lambda x: (float(x) + float(wave) / 1000) * 
                                                                 light_speed / wave).round(0)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Photospheric Velocity Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

for wave in dict_markers.keys():
    temp_series = vel_df[str(wave)].dropna()
    ax.scatter(temp_series.index.values, temp_series.values / 1000, marker=dict_markers[wave][0], 
               color=dict_markers[wave][1], s=45, alpha=0.8, label=dict_markers[wave][2])
    ax.errorbar(temp_series.index.values, temp_series.values / 1000, yerr=vel_df[str(wave) + 'Err'].dropna() / 1000, 
                linestyle='--', capsize=2, capthick=2, color=dict_markers[wave][1], alpha=0.7, label=None)

handles, labels = ax.get_legend_handles_labels()
handles = [handles[0], handles[5], handles[2], handles[1], handles[3], handles[4]]
labels = [labels[0], labels[5], labels[2], labels[1], labels[3], labels[4]]
ax.legend(handles, labels, fontsize=14, markerscale=2, loc=1, frameon=False)

set_plotparams(ax, ylim=12.5)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax.set_ylabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=16)

fig.savefig('OUTPUT_PlotPhotVel.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Photospheric Velocity Inferred From Halpha Line & FeII Line (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
list_halpha = group_similar_files('', DIR_SNe + 'PhotVel_Data/*Halpha.txt')
list_phot = group_similar_files('', DIR_SNe + 'PhotVel_Data/*Phot.txt')

fig2 = plt.figure(figsize=(8, 16))
ax2 = fig2.add_subplot(211)
ax3 = fig2.add_subplot(212)

for file_name in list_halpha:
    name  = file_name.split('/')[-1].split('.')[0].split('_')[0]
    temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', engine='python', names=['Phase', 'Vel'])
    temp_df = temp_df[temp_df['Phase'] < 180]
    ax2.plot(temp_df['Phase'], temp_df['Vel'] / 1000, linestyle='--', marker=dict_snmark[name][0], ms=10,
             c=dict_snmark[name][1], alpha=0.8, label=name)
    
for file_name in list_phot:
    name  = file_name.split('/')[-1].split('.')[0].split('_')[0]
    temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', engine='python', names=['Phase', 'Vel'])
    temp_df = temp_df[temp_df['Phase'] < 180]
    ax3.plot(temp_df['Phase'], temp_df['Vel'] / 1000, linestyle='--', marker=dict_snmark[name][0], ms=10,
             c=dict_snmark[name][1], alpha=0.8, label=name)

temp_series = vel_df['6563'].dropna()
ax2.plot(temp_series.index.values, temp_series / 1000, marker='*', ms=10, linestyle='-', color='k', label='ASASSN-14dq')
temp_series = vel_df['5169'].dropna()
ax3.plot(temp_series.index.values, temp_series / 1000, marker='*', ms=10, linestyle='-', color='k', label='ASASSN-14dq')

set_plotparams(ax2)
ax2.set_ylim(1.5, 13.5)
ax2.legend(fontsize=14, markerscale=1.5, frameon=False)
ax2.set_xticklabels([])
ax2.set_ylabel(r'$\rm V_{H\alpha}\ [x10^3\ km\ s^{-1}]$', fontsize=16)

set_plotparams(ax3)
ax3.legend(fontsize=14, markerscale=1.5, frameon=False)
ax3.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax3.set_ylabel(r'$\rm V_{FeII}\ [x10^3\ km\ s^{-1}]$', fontsize=16)

plt.subplots_adjust(hspace=0)
fig2.savefig('OUTPUT_PlotCompPhotVel.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig2)

# ------------------------------------------------------------------------------------------------------------------- #
