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

a = np.ar

# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
file_photwave = 'OUTPUT_SpecPhotMinima'
date_explosion = 2457644.60
light_speed = 2.99792458e5  # In Km/s
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
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
    ax_obj.set_xlim(-2, 145)
    ax_obj.set_ylim(1.5, ylim)
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(25))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(2.5))
    ax_obj.yaxis.set_major_locator(MultipleLocator(2))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax_obj.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=14)
    ax_obj.tick_params(which='minor', direction='in', width=0.7, length=4, labelsize=14)


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Pandas DataFrame With Velocity Information From Different Spectral Lines
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_SPEC + file_photwave, sep='\s+')
data_df = data_df.set_index('Date').replace('INDEF', np.nan)
wave_df = pd.DataFrame(index=data_df.index.values)
wave_df.index.name = 'Date'

for column in data_df.columns.values:
    wave_df[column] = data_df[column].apply(lambda x: calc_meanstd(x)[0])
    wave_df[column + 'Err'] = data_df[column].apply(lambda x: calc_meanstd(x)[1])

wave_df.round(2).to_csv('OUTPUT_PhotVel', sep=' ', index=True, header=True, na_rep='INDEF')

wave_df = wave_df.reset_index()
wave_df['JD'] = wave_df['Date'].apply(lambda x: cald_to_jd(x))
wave_df['Phase'] = wave_df['JD'] - date_explosion
wave_df = wave_df.set_index('Phase', drop=True)

dict_markers = {6563: ['*', 'k', r'H$\alpha$'], 4861: ['v', 'red', r'H$\beta$'], 5169: ['s', 'g', r'$\rm Fe\,II$ 5169'],
                5018: ['+', 'b', r'$\rm Fe\,II\ 5018$'], 4924: ['p', 'orange', r'$\rm Fe\,II$ 4924'],
                4340: ['o', 'teal', r'H$\gamma$']}
# 6142: ['o', 'orange', r'Ba$\rm \,II$ 6142'], 6246: ['*', 'teal', r'Sc$\rm \,II$ 6246']

vel_df = pd.DataFrame()
vel_df['Date'] = wave_df['Date'].copy()

for wave in dict_markers.keys():
    vel_df[str(wave)] = wave_df[str(wave)].apply(lambda x: (wave - float(x)) * light_speed / wave).round(0)
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
               c=dict_markers[wave][1], s=80, alpha=0.8, label=dict_markers[wave][2])
    ax.errorbar(temp_series.index.values, temp_series.values / 1000, yerr=vel_df[str(wave) + 'Err'].dropna() / 1000,
                ls='--', lw=1, capsize=2, capthick=2, c=dict_markers[wave][1], label=None)

handles, labels = ax.get_legend_handles_labels()
handles = [handles[0], handles[5], handles[2], handles[1], handles[3], handles[4]]
labels = [labels[0], labels[5], labels[2], labels[1], labels[3], labels[4]]
ax.legend(handles, labels, fontsize=14, markerscale=1.5, loc=1, frameon=False)

set_plotparams(ax, ylim=16)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax.set_ylabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=16)

fig.savefig('PLOT_PhotVel.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
