#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxx---------------TABULTATE THE SN MAGNITUDE AND PLOT COLOR CURVES----------------xxxxxxxxxxxxxxxxxx #
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
JD_offset = 2457600
clip_epoch = 200

filters = ['U', 'B', 'V', 'R', 'I']
colors = ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
DIR_PHOT = "/home/avinash/Supernovae_Data/Photometry/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS and Other Type II SNe
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')

sndata_df = pd.read_csv(DIR_SNe + 'LC_Data/TypeIISNe.dat', sep='\s+', comment='#')
sndata_df = sndata_df.replace('INDEF', np.nan).set_index(['Name', 'Marker', 'Color']).astype('float64')
sndata_df = sndata_df.reset_index().set_index('Name')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of The SNe In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
EBV_mag = 0.0865
EBV_err = 0.0018
dist_val = 36.0
dist_err = 2.5
distmod_mag = 32.78
distmod_err = 0.15
redshift = 0.00806
date_explosion = 2457644.60


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
    cal_date = date(*time_tuple[0:3]).strftime('%Y-%m-%d')

    return cal_date


def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split('-')
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
            append_data = [round((val_1 ** 2 + val_2 ** 2) ** 0.5, int(precision))
                           if val_1 != 'INDEF' and val_2 != 'INDEF' else 'INDEF'
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
            input_df[band] = input_df[band].astype('float64') - filter_df.loc[band, 'RLambda'] * EBV_mag
        else:
            input_df[band] = input_df[band].astype('float64'). \
                apply(lambda x: (x ** 2 + (filter_df.loc[band, 'RLambda'] * EBV_err) ** 2) ** 0.5)

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
        output_df['B-V'] = add_series([input_df['B'], input_df['V']], err=True)
        output_df['U-B'] = add_series([input_df['U'], input_df['B']], err=True)
        output_df['V-R'] = add_series([input_df['V'], input_df['R']], err=True)
        output_df['R-I'] = add_series([input_df['R'], input_df['I']], err=True)
        output_df['V-I'] = add_series([input_df['V'], input_df['I']], err=True)

    output_df[colors] = output_df[colors].apply(pd.to_numeric, errors='coerce').round(int(precision))

    return output_df


def get_colorframe(name, data_df):
    """
    Converts a column-wise magnitude Pandas DataFrame to a row-wise Pandas DataFrame.
    Args:
        name        : Name of the SNe whose data is read
        data_df     : Input Pandas DataFrame with column-wise magnitudes
    Returns:
        mag_df      : Output Pandas DataFrame with color terms
        err_df      : Output Pandas DataFrame with errors on color terms
    """
    data_df = data_df.set_index('JD')
    mag_df = data_df[[x for x in data_df.columns.values if 'Err' not in x]].copy()
    err_df = data_df[['Date', 'Phase'] + [x for x in data_df.columns.values if 'Err' in x]].copy()
    err_df = err_df.rename(columns=lambda x: x.strip('Err'))

    for band in filters:
        if band not in data_df.columns.values:
            mag_df[band] = np.nan
            err_df[band] = np.nan
        else:
            mag_df[band] = mag_df[band].astype('float64') - filter_df.loc[band, 'RLambda'] * sndata_df.loc[name, 'EBV']

    mag_df = calc_colormagframe(mag_df)
    err_df = calc_colormagframe(err_df, err=True)

    return mag_df, err_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Subplots
# ------------------------------------------------------------------------------------------------------------------- #

def plot_color(ax_obj, input_df, name, color):
    """
    Plots the color evolution from a Pandas DataFrame containing color magnitudes for SNe to be compared.
    Args:
        ax_obj   : Axes object to be used for plotting
        input_df : Input Pandas DataFrame containing broadband magnitudes
        name     : Name of the SNe whose color terms are being plotted
        color    : Color term to be plotted from the Pandas DataFrame
    Returns:
        None
    """
    mag_df, err_df = get_colorframe(name, input_df)
    #     if color == 'U-B':
    #         mag_df = mag_df[mag_df['Phase'] < 100]
    color_mag = mag_df.set_index('Phase')[color].dropna()

    ax_obj.plot(color_mag.index.values, color_mag.values, marker=sndata_df.loc[name, 'Marker'], markersize=5,
                c=sndata_df.loc[name, 'Color'], linestyle='', alpha=0.9, label=name, markerfacecolor='None')


def plot_sncolor(ax_obj, inpmag_df, inperr_df, color):
    """
    Plots the color evolution from a Pandas DataFrame containing color magnitudes of the SN in study.
    Args:
        ax_obj   : Axes object to be used for plotting and setting plot parameters
        inpmag_df   : Input Pandas DataFrame containing broadband magnitudes
        inperr_df   : Input Pandas DataFrame containing errors in broadband magnitudes
        color       : Color term to be plotted from the Pandas DataFrame
    Returns:
        None
    """
    color_mag = inpmag_df[color].dropna()
    color_err = inperr_df[color].dropna()

    ax_obj.scatter(color_mag.index.values, color_mag.values, marker='*', c='k', s=40, label=name_SN)
    ax_obj.errorbar(color_mag.index.values, color_mag.values, yerr=color_err.values, marker='*', c='k', markersize=8,
                    alpha=0.8, linestyle='-', capsize=2, capthick=1, label=None)

    #     ax_obj.set_ylabel(r'$\rm (' + color + ')_0$', fontsize=14, labelpad=4)
    ax_obj.set_title(r'$\rm (' + color + ')_0$', fontsize=14, loc='center')


def set_plotparams(ax_obj):
    """
    Sets the plot parameters for the color evolution plots.
    Args:
        ax_obj   : Axes object to be used for plotting and setting plot parameters
    Returns:
        None
    """
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_major_locator(MultipleLocator(0.50))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.125))
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.tick_params(which='major', direction='in', length=6, width=1, labelsize=12)
    ax_obj.tick_params(which='minor', direction='in', length=3, width=1, labelsize=12)


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
data_df['Phase'] = (data_df['JD'] - date_explosion).round(int(precision))
data_df['Mag'] = data_df['FMAG'].apply(lambda x: "{:.2f}".format(x)) + r"$\pm$" + \
                 data_df['FERR'].apply(lambda x: "{:.2f}".format(x))

tabular_df = data_df[['FILTER', 'Date', 'JD', 'Phase', 'Mag']]
tabular_df.to_csv('OUTPUT_FinalTabularSNMag', sep=' ', index=False)
display_text("Magnitudes For Supernova Have Been Tabulated Epoch-Wise")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Output The Photometric Magnitudes Onto A Latex Table
# ------------------------------------------------------------------------------------------------------------------- #
data_df = data_df.set_index('Date')
table_df = unorgframe_to_orgframe(data_df, column='Mag')

table_df['Phase'] = table_df['Date'].apply(lambda x: data_df.loc[data_df.index == x, 'Phase'].iloc[0])
table_df['JD'] = table_df['Phase'] + date_explosion - JD_offset
table_df['Phase'] = table_df['Phase'].apply(lambda x: "{:.2f}".format(x) if x < 0 else "+{:.2f}".format(x))

table_df = table_df[['Date', 'JD', 'Phase', 'U', 'B', 'V', 'R', 'I']].sort_values(by='JD')
table_df = table_df.rename(columns={'Phase': 'Phase$^*$', 'U': '$U$', 'B': '$B$', 'V': '$V$', 'R': '$R$', 'I': '$I$'})
table_df = table_df.replace(np.nan, '---', regex=True)
table_df.to_latex('_PhotTable.tex', escape=False, index=False)
display_text("Photometric Magnitudes (Epoch-Wise) For Supernova Have Been Written Onto A Latex File")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Get Color Terms From UBVRI Photometric Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
snmag_df = ubvriframe_to_colormagframe(unorgframe_to_orgframe(data_df, column='FMAG'))
snerr_df = ubvriframe_to_colormagframe(unorgframe_to_orgframe(data_df, column='FERR'), err=True)

snmag_df = snmag_df[snmag_df['Phase'] < clip_epoch].set_index('Phase')
snerr_df = snerr_df[snerr_df['Phase'] < clip_epoch].set_index('Phase')

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10), sharex=True)

set_plotparams(ax1)
plot_sncolor(ax1, snmag_df, snerr_df, color='B-V')

set_plotparams(ax2)
plot_sncolor(ax2, snmag_df, snerr_df, color='U-B')
ax2.set_yticklabels([])

set_plotparams(ax3)
plot_sncolor(ax3, snmag_df, snerr_df, color='V-R')

set_plotparams(ax4)
plot_sncolor(ax4, snmag_df, snerr_df, color='V-I')
ax4.set_yticklabels([])

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['1999em', '2012ec', '2012aw', '2013ab', '2013ej']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
        datacomp_df = datacomp_df[datacomp_df['Phase'] < clip_epoch]
        plot_color(ax1, datacomp_df, name, 'B-V')
        plot_color(ax2, datacomp_df, name, 'U-B')
        plot_color(ax3, datacomp_df, name, 'V-R')
        plot_color(ax4, datacomp_df, name, 'V-I')

ax1.set_ylim(-1.1, 2.1)
ax2.set_ylim(-1.1, 2.1)

ax3.set_ylim(-0.1, 1.55)
ax4.set_ylim(-0.1, 1.55)

ax4.set_xlim(-10, clip_epoch + 5)
ax1.legend(fontsize=10, markerscale=1.8, loc='lower right', frameon=False)
ax4.set_xlabel('Time Since Explosion [Days]', fontsize=14)

fig.subplots_adjust(hspace=0.09, wspace=0.01)
fig.savefig('PLOT_ColorEvolution.pdf', format='pdf', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
