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
name_SN = '2017iro'
file_photwave = '2017iro_V2.dat'
plot_epoch = 110
host_radvel = 1856                      # In km/s
date_bmax = 2458096.24
light_speed = 2.99792458e5              # In km/s
redshift_corr = False
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/Ib_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2017iro/Spectroscopy/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For File Handling
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
# Useful Functions
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


def calc_meanstd(str_num):
    """
    Calculates mean and standard deviation of a list of values.
    Args:
        str_num : Input string containing numbers separated by a comma
    Returns:
        mean    : Mean of the numbers in the input string
        std     : Standard deviation of the numbers in the input string
    """
    if type(str_num) != float:
        list_val = [float(val) for val in str_num.split(',')]
        return np.mean(list_val), np.std(list_val)
    else:
        return np.nan, np.nan


def wavetovel(wave, linewave):
    return (linewave - float(wave)) * light_speed / linewave

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
    ax_obj.set_ylim(2, ylim)
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(25))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax_obj.yaxis.set_major_locator(MultipleLocator(3))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.3))
    ax_obj.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=14)
    ax_obj.tick_params(which='minor', direction='in', width=0.7, length=4, labelsize=14)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Pandas DataFrame With Velocity Information From Different Spectral Lines
# ------------------------------------------------------------------------------------------------------------------- #
wave_df = pd.read_csv(DIR_SPEC + file_photwave, sep='\s+')
wave_df = wave_df.drop('Date', axis=1).replace('INDEF', np.nan)
wave_df = wave_df.set_index('JD').astype('float64')

vel_df = pd.DataFrame(index=wave_df.index.values)
vel_df.index.name = 'JD'

dict_markers = {5876: ['*', 'k', r'$\rm He\,I$ 5876'], 6678: ['^', 'red', r'$\rm He\,I$ 6678'],
                5169: ['s', 'g', r'$\rm Fe\,II$ 5169'], 6355: ['o', 'orange', r'$\rm Si\,II$ 6355']}

for wave in dict_markers.keys():
    if redshift_corr:
        vel_df[str(wave)] = wave_df[str(wave)].apply(lambda x: wavetovel(x, wave)).round(0)
    else:
        vel_df[str(wave)] = wave_df[str(wave)].apply(lambda x: wavetovel(x, wave) + host_radvel).round(0)

vel_df['Phase'] = vel_df.index - date_bmax
vel_df = vel_df.set_index('Phase')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Photospheric Velocity Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

for wave in dict_markers.keys():
    temp_series = vel_df[str(wave)].dropna()
    ax.scatter(temp_series.index.values, temp_series.values / 1000, marker=dict_markers[wave][0],
               c=dict_markers[wave][1], s=70, alpha=0.7, label=dict_markers[wave][2])

set_plotparams(ax, ylim=12)
ax.set_xlim(-20, plot_epoch)
ax.legend(fontsize=14, markerscale=1.4, loc=1, frameon=False)
ax.set_xlabel('Time Since B-Band Maximum [Days]', fontsize=16)
ax.set_ylabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=16)

fig.savefig('PLOT_PhotVel.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Photospheric Velocity Inferred From Halpha Line & FeII Line (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
dict_snmark = {'2007uy': ['o', 'c'], '2007Y': ['D', 'b'], '2008D': ['s', 'g'], '2009jf': ['^', 'r'],
               'iPTF13bvn': ['P', 'orange'], '2012au': ['h', 'sienna']}

list_feii = group_similar_files('', DIR_SNe + 'PhotVel_Data/*FeII5169.dat')
list_hei = group_similar_files('', DIR_SNe + 'PhotVel_Data/*HeI5876.dat')

fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8), sharey=True)

for wave in dict_markers.keys():
    temp_series = vel_df[str(wave)].dropna()
    ax1.scatter(temp_series.index.values, temp_series.values / 1000, marker=dict_markers[wave][0],
                c=dict_markers[wave][1], s=100, alpha=0.8, label=dict_markers[wave][2])

temp_series = vel_df['5169'].dropna()
ax2.plot(temp_series.index.values, temp_series / 1000, marker='*', ms=15, ls='-', c='k', label=name_SN)

for file_name in list_feii:
    name = file_name.split('/')[-1].split('.')[0].split('_')[0]
    temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', names=['Phase', 'Vel'])
    temp_df = temp_df[temp_df['Phase'] < plot_epoch]
    ax2.plot(temp_df['Phase'], temp_df['Vel'] / 1000, ls='--', marker=dict_snmark[name][0], ms=10,
             c=dict_snmark[name][1], alpha=0.8, label=name)

temp_series = vel_df['5876'].dropna()
ax3.plot(temp_series.index.values, temp_series / 1000, marker='*', ms=15, ls='-', c='k', label=name_SN)

for file_name in list_hei:
    name = file_name.split('/')[-1].split('.')[0].split('_')[0]
    temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', names=['Phase', 'Vel'])
    temp_df = temp_df[temp_df['Phase'] < plot_epoch]
    ax3.plot(temp_df['Phase'], temp_df['Vel'] / 1000, ls='--', marker=dict_snmark[name][0], ms=10,
             c=dict_snmark[name][1], alpha=0.8, label=name)


set_plotparams(ax1)
ax1.text(20, 15, 'SN 2017iro', fontsize=16)
ax1.set_xlim(-20, plot_epoch)
ax1.legend(fontsize=14, markerscale=1.5, loc=1, frameon=False)
ax1.set_xlabel('Time Since B-Band Maximum [Days]', fontsize=16)
ax1.set_ylabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=16)

set_plotparams(ax2)
ax2.set_xlim(-20, plot_epoch)
ax2.legend(fontsize=14, markerscale=1.5, frameon=False)
ax2.text(20, 15, r'$\rm Fe\,II$ 5169', fontsize=16)
ax2.set_xlabel('Time Since B-Band Maximum [Days]', fontsize=16)

set_plotparams(ax3, ylim=16.7)
ax3.set_xlim(-20, plot_epoch)
ax3.legend(fontsize=14, markerscale=1.5, frameon=False)
ax3.text(20, 15, r'$\rm He\,I$ 5876', fontsize=16)
ax3.set_xlabel('Time Since B-Band Maximum [Days]', fontsize=16)

fig2.subplots_adjust(wspace=0.01)
fig2.savefig('PLOT_CompPhotVel.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #

