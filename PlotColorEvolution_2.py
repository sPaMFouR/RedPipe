#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxx------------TABULATE THE SUPERNOVA MAGNITUDES & PLOT COLOR EVOLUTION---------xxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
from datetime import date
import matplotlib.pyplot as plt
from jdcal import jd2gcal, gcal2jd
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 2
epoch = 2400000.5
filters = ['U', 'B', 'V', 'R', 'I']
colors = ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Central Wavelength & Zero Point Correction For Different Photometric Bands
# ------------------------------------------------------------------------------------------------------------------- #
dict_centralwav = {'U': 3700, 'B': 4200, 'V': 5300, 'R': 6000, 'I': 8050,
                   'u': 3655, 'g': 5105, 'r': 6480, 'i': 7105, 'z': 8640, 'Z': 7500}
dict_zp = {'U': -0.152, 'B': -0.602, 'V': 0.000, 'R': 0.555, 'I': 1.271,
           'u': -0.002, 'g': -0.455, 'r': 0.393, 'i': 1.028, 'z': 1.577}
dict_flambda = {'U': 417.5e-11, 'B': 632e-11, 'V': 363.1e-11, 'R': 217.7e-11, 'I': 112.6e-11,
                'u': 363.9e-11, 'g': 552.1e-11, 'r': 252.9e-11, 'i': 140.9e-11, 'z': 85.01e-11}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
dict_R = {'U': 4.95, 'B': 4.27, 'V': 3.15, 'R': 2.65, 'I': 1.72,
          'u': 5.02, 'g': 3.32, 'r': 2.37, 'i': 2.08, 'z': 1.54}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
DIR_PHOT = "/home/avinash/Supernovae_Data/PhotTemp/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of The SNe In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SNe = '2016gfy'
EBV_mag = 0.0865
EBV_err = 0.0018
dist_val = 36.0
dist_err = 2.5
distmod_mag = 32.78
distmod_err = 0.15
redshift = 0.00806
JD_offset = 2457600.
date_explosion = 2457644.60
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Other SN Details
# ------------------------------------------------------------------------------------------------------------------- #
# Leonard 1999em, Bose 2012aw, Bose 2013ej, Leonard 1999gi, 
# Sahu 2004et, Pastorello 2005cs, Bose 2013ab, Inserra 2009aw
dict_explosion = {'1999em': 2451475.6, '2012aw': 2456002.59, '2013ej': 2456497.30, '1999gi': 2451518.22, 
                  '2004et': 2453270.25, '2005cs': 2453549.0, '2013ab': 2456340.0, '2009bw': 2454916.5}
dict_maximum = {'1999em': 2451475.6, '2012aw': 2456002.59, '2013ej': 2456497.30, '1999gi': 2451518.22, 
                '2004et': 2453270.25, '2005cs': 2453549.0, '2013ab': 2456340.0, '2009bw': 2454925.3}
dict_EBV = {'1999em': [0.10, 0.05], '2012aw': [0.07, 0.01], '2013ej': [0.06, 0.001], '1999gi': [0.21, 0.09], 
            '2004et': [0.41, 0.0], '2005cs': [0.05, 0.00], '2013ab': [0.044, 0.066], '2009bw': [0.31, 0.03]}
dict_Mv = {'1999em': [-15.9, 0.2], '2012aw': [-16.67, 0.04], '2013ej': [-16.6, 0.1], '1999gi': [-16.4, 0.6], 
           '2004et': [-17.14, 0.0], '2005cs': [-15.2, 0.0], '2013ab': [-16.7, 0.0], '2009bw': [-16.87, 0.16]}
dict_dist = {'1999em': [8.2, 0.6], '2012aw': [9.9, 0.1], '2013ej': [9.57, 0.70], '1999gi': [13.3, 0.6], 
             '2004et': [5.6, 0.1], '2005cs': [8.9, 0.5], '2013ab': [24.3, 1.0], '2009bw': [20.2, 1.5]}
dict_snmark = {'1999em': ['o', 'c'], '2012aw': ['D', 'b'], '2013ej': ['s', 'g'], '1999gi': ['^', 'r'], 
               '2004et': ['p', 'm'], '2005cs': ['*', 'y'], '2013ab': ['P', 'coral'], '2009bw': ['X', 'violet']}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling & Displaying Task Completion
# ------------------------------------------------------------------------------------------------------------------- #

def group_similar_files(text_list, common_text, exceptions=''):
    """s
    Groups similar files based on the string "common_text". Writes the similar files
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


def display_text(text_to_display):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text_to_display : Text to be displayed
    Returns:
        None
    """
    print ("\n" + "# " + "-" * (12 + len(text_to_display)) + " #")
    print ("# " + "-" * 5 + " " + str(text_to_display) + " " + "-" * 5 + " #")
    print ("# " + "-" * (12 + len(text_to_display)) + " #" + "\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Convert Julian Date Into Calendar Date In String Format
# ------------------------------------------------------------------------------------------------------------------- #

def jd_to_cald(julian_day):
    """
    Converts julian day into calendar day in string format.
    Args:
        julian_day  : Julian day value to be converted to calendar day
    Returns:
        cal_date    : Calendar date corresponding to input julian day
    """
    time_tuple = jd2gcal(epoch, julian_day - epoch)
    cal_date = date(*time_tuple[0:3]).strftime("%Y-%m-%d")

    return cal_date


def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split("-")
    jd_tuple = gcal2jd(date_comp[0], date_comp[1], date_comp[2])
    julian_day = jd_tuple[0] + jd_tuple[1]

    return julian_day

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Manipulating Pandas DataFrames
# ------------------------------------------------------------------------------------------------------------------- #

def add_series(list_series, sub=False, err=False):
    """
    Adds multiple Pandas Series column wise and obtains a resultant Pandas Series.
    Args:
        list_series     : List of all Pandas Series to be added to obtain a single Pandas Series
        sub             : True, if the series needs to be subtracted
        err             : True, if the series contains error data
    Returns:
        output_series   : Output Pandas Series obtained after adding all the series
    """
    output_series = list_series[0]
    list_indices = output_series.index.values

    if err:
        sub = False

    for series in list_series[1:]:
        if not err:
            if not sub:
                append_data = [val_1 + val_2 if val_1 != 'INDEF' and val_2 != 'INDEF' else 'INDEF'
                               for val_1, val_2 in zip(output_series, series)]
            else:
                append_data = [val_1 - val_2 if val_1 != 'INDEF' and val_2 != 'INDEF' else 'INDEF'
                               for val_1, val_2 in zip(output_series, series)]
        else:
            append_data = [round((val_1 ** 2 + val_2 ** 2) ** 0.5,
                                 int(precision)) if val_1 != 'INDEF' and val_2 != 'INDEF' else 'INDEF'
                           for val_1, val_2 in zip(output_series, series)]

        output_series = pd.Series(data=append_data, index=list_indices)

    return output_series


def add_errseries(list_series):
    """
    Adds multiple Pandas Series containing error data in quadrature and obtains a resultant Pandas Series.
    Args:
        list_series     : List of all Pandas Series to be added to obtain a single Pandas Series
    Returns:
        output_series   : Output Pandas Series obtained after adding all the series
    """
    output_series = list_series[0]
    list_indices = output_series.index.values

    for series in list_series[1:]:
        append_data = [round((val_1 ** 2 + val_2 ** 2) ** 0.5, int(precision))
                       for val_1, val_2 in zip(output_series, series)]

        output_series = pd.Series(data=append_data, index=list_indices)

    return output_series


def unorgframe_to_orgframe(input_df, column):
    """
    Creates a Pandas DataFrame with magnitudes arranged epoch-wise from a DataFrame with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes
        column      : Column to be extracted from the Pandas DataFrame
    Returns:
        output_df   : Pandas DataFrame containing organised broadband magnitudes
    """
    dict_val = {}
    for index, row in input_df.iterrows():
        if index not in dict_val.keys():
            dict_val[index] = {}
        dict_val[index][row['FILTER']] = row[column]

    output_df = pd.DataFrame(dict_val).T
    output_df.index.name = 'Date'
    output_df = output_df.reset_index()
    output_df['Phase'] = output_df['Date'].apply(lambda x: input_df.loc[input_df.index == x, 'Phase'].iloc[0])

    return output_df


def ubvriframe_to_colormagframe(input_df, err=False):
    """
    Creates a Pandas DataFrame with color terms from an input DataFrame with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes and color terms
        err         : True, if the dataframe contains error data
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes and color terms
    """
    output_df = input_df[['Date', 'Phase']].copy()

    for band in [x for x in filters if x in input_df.columns.values]:
        if not err:
            input_df[band] = input_df[band].astype('float64') - dict_R[band] * EBV_mag
        else:
            input_df[band] = input_df[band].astype('float64').apply(lambda x: (x ** 2 +
                                                                               (dict_R[band] * EBV_err) ** 2) ** 0.5)
    output_df['B-V'] = add_series([input_df['B'], input_df['V']], sub=True, err=err)
    output_df['U-B'] = add_series([input_df['U'], input_df['B']], sub=True, err=err)
    output_df['V-R'] = add_series([input_df['V'], input_df['R']], sub=True, err=err)
    output_df['R-I'] = add_series([input_df['R'], input_df['I']], sub=True, err=err)
    output_df['V-I'] = add_series([input_df['V'], input_df['I']], sub=True, err=err)

    output_df[colors].apply(pd.to_numeric, errors='coerce').round(int(precision))

    return output_df


def calc_colormagframe(input_df, err=False):
    """
    Creates a Pandas DataFrame with color terms from an input DataFrame with broadband magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes and color terms
        err         : True, if the dataframe contains error data
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes and color terms
    """
    output_df = input_df[['Date', 'Phase']].copy()
    input_df = input_df.drop(['Date', 'Phase'], axis=1).apply(pd.to_numeric)

    if not err:
        output_df['B-V'] = input_df['B'] - input_df['V']
        output_df['U-B'] = input_df['U'] - input_df['B']
        output_df['V-R'] = input_df['V'] - input_df['R']
        output_df['R-I'] = input_df['R'] - input_df['I']
        output_df['V-I'] = input_df['V'] - input_df['I']
    else:
        output_df['B-V'] = add_errseries([input_df['B'], input_df['V']])
        output_df['U-B'] = add_errseries([input_df['U'], input_df['B']])
        output_df['V-R'] = add_errseries([input_df['V'], input_df['R']])
        output_df['R-I'] = add_errseries([input_df['R'], input_df['I']])
        output_df['V-I'] = add_errseries([input_df['V'], input_df['I']])

    output_df[colors] = output_df[colors].apply(pd.to_numeric, errors='coerce').round(int(precision))

    return output_df


def get_colorframe(name, data_df):
    """
    Converts a column-wise magnitude Pandas DataFrame to a row-wise Pandas DataFrame.
    Args:
        name        : Name of the SNe whose data is read
        data_df     : Input Pandas DataFrame
    Returns:
        output_df   : Output Pandas DataFrame
    """
    data_df = data_df.set_index('JD')
    mag_df = data_df[[x for x in data_df.columns.values if 'Err' not in x]].copy()
    err_df = data_df[['Date', 'Phase'] + [x for x in data_df.columns.values if 'Err' in x]].copy()
    err_df = err_df.rename(columns=lambda x: x.strip('Err'))

    for band in [x for x in filters if x not in data_df.columns.values]:
        mag_df[band] = np.nan
        err_df[band] = np.nan

    for band in [x for x in filters if x in mag_df.columns.values]:
        mag_df[band] = mag_df[band].astype('float64') - dict_R[band] * dict_EBV[name][0]

    mag_df = calc_colormagframe(mag_df)
    err_df = calc_colormagframe(err_df, err=True)

    return mag_df, err_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Subplots
# ------------------------------------------------------------------------------------------------------------------- #

def plot_color(ax_obj, input_df, name, color):
    """
    Plots the color evolution from a Pandas DataFrame containing color magnitudes for different phases of the SN.
    Args:
        ax_obj   : Axes object to be used for plotting and setting plot parameters
        input_df : Input Pandas DataFrame containing broadband magnitudes
        name     : Name of the SNe whose color terms are being plotted
        color    : Color term to be plotted from the Pandas DataFrame
    Returns:
        None
    """
    mag_df, err_df = get_colorframe(name, input_df)
    color_mag = mag_df.set_index('Phase')[color].dropna()
    ax_obj.plot(color_mag.index.values, color_mag.values, marker=dict_snmark[name][0], c=dict_snmark[name][1], 
                markersize=4, linestyle='--', label=name)

    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(25))
    ax_obj.set_ylabel(r'$(' + color + ')_0$', fontsize=12, labelpad=6)
    ax_obj.tick_params(which='both', direction='in', width=0.8, labelsize=12)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Other SNe Data From Archive Folder
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', DIR_SNe + 'LC_Data/*.asc')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Tabulate The Photometric Data Epoch-Wise Onto An Output File
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv('OUTPUT_FinalSNMagTemp', sep='\s+', engine='python')

data_df['Date'] = data_df['JD'].apply(lambda x: jd_to_cald(x))
data_df['Phase'] = data_df['JD'] - date_explosion
data_df['Phase'] = data_df['Phase'].round(int(precision))
data_df['Mag'] = data_df['FMAG'].apply(lambda x: '{:.2f}'.format(x)) + r'$\pm$' + \
                 data_df['FERR'].apply(lambda x: '{:.2f}'.format(x))

tabular_df = data_df[['FILTER', 'Date', 'JD', 'Phase', 'Mag']]
tabular_df.to_csv('OUTPUT_FinalTabularSNMag', sep=' ', index=False)
display_text("Magnitudes For Supernova Have Been Tabulated Epoch-Wise")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Output The Photometric Magnitudes Onto A Latex Table
# ------------------------------------------------------------------------------------------------------------------- #
# data_df = data_df[~((data_df['Date'] == '2017-01-13') & (data_df['FILTER'] == 'U'))]
# data_df = data_df[~((data_df['Date'] == '2016-12-29') & (data_df['FILTER'] == 'U'))]
# data_df = data_df[~((data_df['Date'] == '2017-02-12') & (data_df['FILTER'] == 'B'))]
# data_df = data_df[~((data_df['Date'] == '2017-02-18') & (data_df['FILTER'] == 'B'))]
# data_df = data_df[~((data_df['Date'] == '2017-02-09') & (data_df['FILTER'] == 'I'))]
# data_df = data_df[~((data_df['Date'] == '2017-08-11') & (data_df['FILTER'] == 'I'))]
# data_df = data_df[~((data_df['Date'] == '2017-01-28') & (data_df['FILTER'] == 'V'))]

# data_df = data_df[data_df['Phase'] < 250]
max_epoch = data_df['Phase'].max()

data_df = data_df.set_index('Date')
table_df = unorgframe_to_orgframe(data_df, column='Mag')

table_df['Phase'] = table_df['Date'].apply(lambda x: data_df.loc[data_df.index == x, 'Phase'].iloc[0])
table_df['JD'] = table_df['Phase'] + date_explosion - JD_offset
table_df['Phase'] = table_df['Phase'].apply(lambda x: '{:.2f}'.format(x) if x < 0 else '+{:.2f}'.format(x))
table_df = table_df.sort_values(by='Phase')

table_df = table_df.sort_values(by='JD')
table_df = table_df.replace(np.nan, '---', regex=True)
table_df = table_df[['Date', 'JD', 'Phase', 'U', 'B', 'V', 'R', 'I']]
table_df = table_df.rename(columns={'Phase': 'Phase$^*$', 'U': '$U$', 'B': '$B$', 'V': '$V$', 'R': '$R$', 'I': '$I$'})
table_df.to_latex('photsn.tex', escape=False, index=False)

display_text("Photometric Magnitudes (Epoch-Wise) For Supernova Have Been Written Onto A Latex File")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Get Color Terms From UBVRI Photometric Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
snmag_df = ubvriframe_to_colormagframe(unorgframe_to_orgframe(data_df, column='FMAG'))
snerr_df = ubvriframe_to_colormagframe(unorgframe_to_orgframe(data_df, column='FERR'), err=True)
snmag_df = snmag_df.set_index('Phase')
snerr_df = snerr_df.set_index('Phase')

fig = plt.figure(figsize=(6, 11))

ax1 = fig.add_subplot(511)
color_mag = snmag_df['B-V'].dropna()
color_err = snerr_df['B-V'].dropna()
ax1.set_ylim(-0.4, 2.1)
ax1.yaxis.set_major_locator(MultipleLocator(0.50))
ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.scatter(color_mag.index.values, color_mag.values, marker='o', c='k', s=12, alpha=0.5, label=None)
ax1.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='o', c='k', markersize=3, 
             alpha=0.5, linestyle='-', capsize=2, capthick=1, label=name_SNe)

ax2 = fig.add_subplot(512, sharex=ax1)
color_mag = snmag_df['U-B'].dropna()
color_err = snerr_df['U-B'].dropna()
ax2.yaxis.set_major_locator(MultipleLocator(1.00))
ax2.yaxis.set_minor_locator(MultipleLocator(0.50))
ax2.scatter(color_mag.index.values, color_mag.values, marker='o', c='k', s=12, alpha=0.5, label=None)
ax2.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='o', c='k', markersize=3, 
             alpha=0.5, linestyle='-', capsize=2, capthick=1, label=name_SNe)

ax3 = fig.add_subplot(513, sharex=ax1)
color_mag = snmag_df['V-R'].dropna()
color_err = snerr_df['V-R'].dropna()
ax3.set_ylim(-0.1, 1.2)
ax3.yaxis.set_major_locator(MultipleLocator(0.50))
ax3.yaxis.set_minor_locator(MultipleLocator(0.25))

ax3.scatter(color_mag.index.values, color_mag.values, marker='o', c='k', s=12, alpha=0.5, label=None)
ax3.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='o', c='k', markersize=3, 
             alpha=0.5, linestyle='-', capsize=2, capthick=1, label=name_SNe)

ax4 = fig.add_subplot(514, sharex=ax1)
color_mag = snmag_df['R-I'].dropna()
color_err = snerr_df['R-I'].dropna()
ax4.set_ylim(-0.25, 0.5)
ax4.yaxis.set_major_locator(MultipleLocator(0.20))
ax4.yaxis.set_minor_locator(MultipleLocator(0.10))

ax4.scatter(color_mag.index.values, color_mag.values, marker='o', c='k', s=12, alpha=0.5, label=None)
ax4.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='o', c='k', markersize=3, 
             alpha=0.5, linestyle='-', capsize=2, capthick=1, label=name_SNe)

ax5 = fig.add_subplot(515, sharex=ax1)
ax5.yaxis.set_major_locator(MultipleLocator(0.50))
ax5.yaxis.set_minor_locator(MultipleLocator(0.25))
color_mag = snmag_df['V-I'].dropna()
color_err = snerr_df['V-I'].dropna()
ax5.scatter(color_mag.index.values, color_mag.values, marker='o', c='k', s=12, alpha=0.5, label=None)
ax5.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='o', c='k', markersize=3, 
             alpha=0.5, linestyle='-', capsize=2, capthick=1, label=name_SNe)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name not in ['2005cs']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        datacomp_df = datacomp_df.replace('INDEF', np.nan)
        datacomp_df = datacomp_df.sort_values(by='Phase')[datacomp_df['Phase'] < 250]
        plot_color(ax1, datacomp_df, name, 'B-V')
        plot_color(ax2, datacomp_df, name, 'U-B')
        plot_color(ax3, datacomp_df, name, 'V-R')
        plot_color(ax4, datacomp_df, name, 'R-I')
        plot_color(ax5, datacomp_df, name, 'V-I')

ax2.legend(fontsize=9, markerscale=2, loc='upper right')
ax5.set_xlim(-10, max_epoch + 20)
ax5.set_xlabel('Time Since Explosion [Days]', fontsize=12)
fig.subplots_adjust(hspace=0.001, top=0.9, right=0.95)
fig.savefig('OUTPUT_PlotColor.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #

