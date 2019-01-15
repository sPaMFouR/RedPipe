#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxx--------------------FIT SPLINE TO SUPERNOVA LIGHT CURVE-----------------------xxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
from scipy import interpolate
from datetime import date
import matplotlib.pyplot as plt
from jdcal import jd2gcal, gcal2jd
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
name_SNe = '2016gfy'
precision = 3
first_epoch = 0
fit_epoch = 385
last_epoch = 385
epoch = 2400000.5
date_explosion = 2457641.40
filters = ['U', 'B', 'V', 'R', 'I']
dict_markers = {'U': ['o', 'k'], 'B': ['D', 'b'], 'V': ['s', 'g'], 'R': ['^', 'r'],
                'I': ['p', 'm'], 'g': ['*', 'y'], 'r': ['8', 'c'], 'i': ['.', 'indigo']}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Convert Julian Date Into Calendar Date In String Format & Vice Versa
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
# Function For Writing Interpolated Supernova Magnitudes Onto OUTPUT Files For Different Bands
# ------------------------------------------------------------------------------------------------------------------- #

def write_snmag(func, epoch_series, band, output_prefix="OUTPUT_InterpSNMag_"):
    """
    Writes interpolated SN Magnitudes onto a text file for the specified function and the array of dates.
    Args:
        func            : Functional fit to the SN light curve
        epoch_series    : Pandas Series of epochs at which the object was observed
        band            : Band of observation of the SN
        output_prefix   : Prefix of the output file onto which the magnitudes have to be written
    Returns:
        data_df         : Pandas DataFrame containing interpolated supernova magnitudes for a single band
    """
    data_df = pd.DataFrame()
    data_df['Phase'] = np.arange(first_epoch, epoch_series.max(), 0.5)
    data_df['JD'] = data_df['Phase'] + date_explosion
    data_df['Mag'] = data_df['Phase'].apply(lambda x: "{0:.3f}".format(float(func(x))))
    data_df = data_df.round(int(precision))
    data_df.to_csv(str(output_prefix) + str(band), sep=' ', index=False)

    return data_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions To Fit Spline Over The Light Curve
# ------------------------------------------------------------------------------------------------------------------- #

def fit_spline1d(ax_obj, list_jd, list_mag, label, smooth=0.01):
    """
    Fits 1-Dimensional spline onto the SN light curve.
    Args:
        ax_obj      : Axes object on which the plot parameters are to be set
        list_jd     : List of JD at which SN was observed
        list_mag    : List of magnitude corresponding to the epoch (JD) the SN was observed
        label       : Label to be given to the light curve
        smooth      : Smoothing parameter to be used for the spline fit
    Returns:
        spline      : Functional form of the spline function fitted to the light curve
    """
    xaxis = np.linspace(list_jd.min(), list_jd.max(), 1000)
    spline = interpolate.UnivariateSpline(list_jd, list_mag, k=1, s=smooth)

    ax_obj.plot(xaxis, spline(xaxis), 'g-')
    ax_obj.scatter(list_jd, list_mag, label=label, marker=dict_markers[label][0], c=dict_markers[label][1])

    return spline


def fit_spline2d(ax_obj, list_jd, list_mag, label, smooth=0.5):
    """
    Fits 2-Dimensional spline onto the SN light curve.
    Args:
        ax_obj      : Axes object on which the plot parameters are to be set
        list_jd     : List of JD at which SN was observed
        list_mag    : List of magnitude corresponding to the epoch (JD) the SN was observed
        label       : Label to be given to the light curve
        smooth      : Smoothing parameter to be used for the spline fit
    Returns:
        spline      : Functional form of the spline function fitted to the light curve
    """
    xaxis = np.linspace(list_jd.min(), list_jd.max(), 1000)
    spline = interpolate.UnivariateSpline(list_jd, list_mag, k=2, s=smooth)

    ax_obj.plot(xaxis, spline(xaxis), 'k-')
    ax_obj.scatter(list_jd, list_mag, label=label, marker=dict_markers[label][0], c=dict_markers[label][1])

    return spline


def fit_spline3d(ax_obj, list_jd, list_mag, label, smooth=0.04):
    """
    Fits 3-Dimensional spline onto the SN light curve.
    Args:
        ax_obj      : Axes object on which the plot parameters are to be set
        list_jd     : List of JD at which SN was observed
        list_mag    : List of magnitude corresponding to the epoch (JD) the SN was observed
        label       : Label to be given to the light curve
        smooth      : Smoothing parameter to be used for the spline fit
    Returns:
        spline      : Functional form of the spline function fitted to the light curve
    """
    xaxis = np.linspace(list_jd.min(), list_jd.max(), 1000)
    spline = interpolate.UnivariateSpline(list_jd, list_mag, k=3, s=smooth)

    ax_obj.plot(xaxis, spline(xaxis), 'r-')
    ax_obj.scatter(list_jd, list_mag, label=label, marker=dict_markers[label][0], c=dict_markers[label][1])

    return spline


def fit_cubicspline(ax_obj, list_jd, list_mag, label):
    """
    Fits 4-Dimensional spline onto the SN light curve.
    Args:
        ax_obj      : Axes object on which the plot parameters are to be set
        list_jd     : List of JD at which SN was observed
        list_mag    : List of magnitude corresponding to the epoch (JD) the SN was observed
        label       : Label to be given to the light curve
    Returns:
        spline      : Functional form of the spline function fitted to the light curve
    """
    xaxis = np.linspace(list_jd.min(), list_jd.max(), 1000)
    spline = interpolate.CubicSpline(list_jd, list_mag, bc_type='natural')

    ax_obj.plot(xaxis, spline(xaxis), 'b-')
    ax_obj.scatter(list_jd, list_mag, label=label, marker=dict_markers[label][0], c=dict_markers[label][1])

    return spline


def fit_curve(input_df, band, func, smooth):
    """
    Fits spline function specified by input 'func' onto the light curve band specified the input 'band'.
    Args:
        input_df : Pandas DataFrame comprising of the light curve data of a single band
        band     : Band of observation
        func     : Spline function to be used for fitting the light curve
        smooth   : Smoothening parameter to be used for the fit
    Returns:
        None
    """
    global column_df, row_df

    band_df = input_df[input_df['FILTER'] == band].copy()
    funcfit = func(ax, band_df['Phase'], band_df['FMAG'], label=band, smooth=smooth)
    data_df = write_snmag(funcfit, epoch_series=band_df['Phase'], band=band)

    data2_df = data_df.copy()
    data2_df['FILTER'] = band
    row_df = pd.concat([row_df, data2_df], axis=0)

    data_df = data_df.rename(columns={'Mag': band})
    data_df = data_df.set_index('Phase')
    column_df = pd.concat([column_df, data_df[band]], axis=1)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Setting Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plot_params(ax_obj):
    """
    Apply plot parameters to the axes object specified in the input.
    Args:
        ax_obj  : Axes object on which the plot parameters are to be set
    Returns:
        None
    """
    ax.grid()
    ax_obj.legend(fontsize=14, markerscale=2, frameon=False)
    ax_obj.set_ylim(22.5, 15)
    ax_obj.set_xlim(-10, fit_epoch + 20)

    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_major_locator(MultipleLocator(1))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.25))
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))

    ax_obj.set_ylabel('Apparent Magnitude [mag]', fontsize=12)
    ax_obj.set_xlabel('Time Since Explosion [Days]', fontsize=12)
    ax_obj.set_title("Spline Fits To The Light Curves", fontsize=12)
    ax_obj.tick_params(which='both', direction='in', width=0.5, labelsize=12)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit A Spline To Apparent Magnitude Light Curve And Record The Spline Fit Data Onto Log Files
# ------------------------------------------------------------------------------------------------------------------- #
output_df = pd.read_csv("OUTPUT_FinalSNMagTemp", sep='\s+')
output_df['Phase'] = output_df['JD'] - date_explosion
output_df['Date'] = output_df['JD'].apply(lambda x: jd_to_cald(x))

output_df = output_df[output_df['Phase'] < fit_epoch]

output_df = output_df[~((output_df['Date'] == '2017-04-10') & (output_df['FILTER'] == 'V'))]
output_df = output_df[~((output_df['Date'] == '2017-05-05') & (output_df['FILTER'] == 'V'))]
output_df = output_df[~((output_df['Date'] == '2017-10-01') & (output_df['FILTER'] == 'U'))]
output_df = output_df[~((output_df['Date'] == '2016-12-05') & (output_df['FILTER'] == 'U'))]

fig = plt.figure(figsize=(15, 12))
ax = fig.add_subplot(111)

column_df = pd.DataFrame(index=(np.arange(first_epoch, last_epoch, 0.5)))
column_df['JD'] = column_df.index.values + date_explosion
row_df = pd.DataFrame(columns=['Phase', 'JD', 'Mag', 'FILTER'])

fit_curve(output_df, 'U', fit_spline1d, smooth=0.002)
fit_curve(output_df, 'B', fit_spline1d, smooth=0.01)
fit_curve(output_df, 'V', fit_spline3d, smooth=0.03)
fit_curve(output_df, 'R', fit_spline3d, smooth=0.022)
fit_curve(output_df, 'I', fit_spline2d, smooth=0.03)

column_df.to_csv('OUTPUT_InterpSNMagCol', sep=' ', na_rep='INDEF')
row_df.to_csv('OUTPUT_InterpSNMagRow', sep=' ', na_rep='INDEF', index=False)

set_plot_params(ax)
fig.savefig('PLOT_InterpLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
