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
import pwlf
import numpy as np
import pandas as pd
from datetime import date
import uncertainties as unc
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from jdcal import jd2gcal, gcal2jd
from scipy.optimize import minimize
from numpy.polynomial import legendre as ls
from matplotlib.ticker import MultipleLocator

import seaborn as sns
sns.set_style('ticks')
# sns.set_palette(sns.color_palette('Paired', 10))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
epoch_mjd = 2400000.5
JD_offset = 2457600
clip_epoch = 150
fmt_flt = '{0:>7.3f}'
fmt_exp = '{0:>7.4e}'
filters = ['U', 'B', 'V', 'R', 'I']
colors = ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of The SNe In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
EBVMW_mag = 0.07
EBVMW_err = 0.01
EBV_mag = 0.21
EBV_err = 0.05
dist_val = 29.64
dist_err = 2.65
redshift = 0.008059
date_explosion = 2457641.40
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
DIR_PHOT = "/home/avinash/Supernovae_Data/2016gfy/Photometry/"
DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
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
# Function For Calculating Reduced ChiSquare
# ------------------------------------------------------------------------------------------------------------------- #

def redchisq(ydata, ymod, sd=np.empty(0), n=8):
    """
    Args:
        ydata    : Observed data
        ymod     : Model data
          sd     : Uncertainties in ydata
          n      : Number of free parameters in the model
    Returns:
        RedChiSq : Reduced Chi-Square
    """
    ydata = np.array(ydata)
    ymod = np.array(ymod)
    sd = np.array(sd)

    if not sd.any():
        chisq = np.sum((ydata - ymod) ** 2)
    else:
        chisq = np.sum(((ydata - ymod) / sd) ** 2)

    nu = ydata.size - 1 - n
    redchisq = chisq / nu

    return redchisq


def calc_mag(mag, err, band, name=name_SN):
    mag = float(mag)
    err = float(err)
    zp = filter_df.loc[band, 'ZeroPoint']
    rlambda = filter_df.loc[band, 'RLambda']

    if name != name_SN:
        dval = sndata_df.loc[name, 'D']
        derr = sndata_df.loc[name, 'DErr']
        ebvmag = sndata_df.loc[name, 'EBV']
        ebverr = sndata_df.loc[name, 'EBVErr']
    else:
        dval = dist_val
        derr = dist_err
        ebvmag = EBV_mag
        ebverr = EBV_err

    distmod_mag = 5 * np.log10(dval * 1e6) - 5
    distmod_err = 5 * np.log10((dval + derr) * 1e6) - 5 - distmod_mag

    absmag = fmt_flt.format(mag - rlambda * ebvmag - distmod_mag)
    abserr = fmt_flt.format((err ** 2 + (rlambda * ebverr) ** 2 + distmod_err ** 2) ** 0.5)

    extcormag = fmt_flt.format(mag - rlambda * ebvmag)
    extcorerr = fmt_flt.format(err)

    return float(absmag), float(abserr), float(extcormag), float(extcorerr)

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
                append_data = [round(val_1 + val_2, int(precision)) if val_1 != 'INDEF' and val_2 != 'INDEF'
                               else 'INDEF' for val_1, val_2 in zip(output_series, series)]
            else:
                append_data = [round(val_1 - val_2, int(precision)) if val_1 != 'INDEF' and val_2 != 'INDEF'
                               else 'INDEF' for val_1, val_2 in zip(output_series, series)]
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


def calc_extcorcolorframe(input_df, err=False, ebv=EBV_mag, ebverr=EBV_err):
    """
    Creates a Pandas DataFrame with color terms from an input DataFrame with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes and color terms
        err         : True, if the dataframe contains error data
        ebv         : The value of reddening adopted
        ebverr      : The error in the reddening adopted
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes and color terms
    """
    output_df = input_df[['Date', 'Phase']].copy()

    for band in [x for x in filters if x in input_df.columns.values]:
        if not err:
            input_df[band] = input_df[band].astype('float64') - filter_df.loc[band, 'RLambda'] * ebv
        else:
            input_df[band] = input_df[band].astype('float64'). \
                apply(lambda x: (x ** 2 + (filter_df.loc[band, 'RLambda'] * ebverr) ** 2) ** 0.5)

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


def organise_sndf(input_df, column):
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


def calc_swiftcolordf(file_name, name):
    """
    Creates a Pandas DataFrame with SWIFT UVOT magnitudes arranged epoch-wise from a file 'file_name'.
    Args:
        file_name   : Text file containing SWIFT UVOT magnitudes
        name        : Name of the supernova for which SWIFT UVOT magntiudes are to be plotted
    Returns:
        mag_df      : Pandas DataFrame containing organised SWIFT UVOT magnitudes
        err_df      : Pandas DataFrame containing organised SWIFT UVOT errors
    """
    input_df = pd.read_csv(file_name, sep='\s+', comment='#', usecols=[0, 1, 2, 3], header=None)
    input_df = input_df.rename(columns={0: 'FILTER', 1: 'MJD', 2: 'FMAG', 3: 'FERR'}).replace('INDEF', np.nan)
    input_df['FILTER'] = input_df['FILTER'].apply(lambda x: 'uv' + x.lower() if 'uv' not in x.lower() else x.lower())

    output_df = input_df.set_index('FILTER').astype('float64').reset_index().copy()
    output_df['JD'] = (output_df['MJD'] + epoch_mjd).round(1)
    output_df['Date'] = output_df['JD'].apply(jd_to_cald)
    output_df['Phase'] = (output_df['JD'] - sndata_df.loc[name, 'DateExp']).round(1)

    for index, band in output_df['FILTER'].iteritems():
        data_magflux = calc_mag(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'],
                                band=band, name=name)
        output_df.loc[index, 'ExtCorMag'] = data_magflux[2]
        output_df.loc[index, 'ExtCorErr'] = data_magflux[3]

    output_df = output_df.set_index('JD').dropna(axis=0, how='any')
    mag_df = organise_sndf(output_df, column='ExtCorMag')
    mag_df['uvw2-uvw1'] = mag_df['uvw2'] - mag_df['uvw1']
    mag_df['uvw2-v'] = mag_df['uvw2'] - mag_df['uvv']
    err_df = organise_sndf(output_df, column='ExtCorErr')
    err_df['uvw2-uvw1'] = (err_df['uvw2'] ** 2 + err_df['uvw1'] ** 2) ** 0.5
    err_df['uvw2-v'] = (err_df['uvw2'] ** 2 + err_df['uvv'] ** 2) ** 0.5

    return mag_df, err_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Subplots
# Plot Confidence Intervals And Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xarr, fcolor='grey'):
    """
    Plots 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xarr    : Array of X-Values over which confidence intervals are to be plotted
        fcolor  : Fill color for the confidence intervals
    Returns:
        None
    """
    coeff = unc.correlated_values(optpar, covpar)
    func = ls.legval(xarr, coeff)
    fit = unp.nominal_values(func)
    sigma = unp.std_devs(func)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr, fitlow, ls='-.', c='k', lw=0.7, alpha=0.3, label='_nolegend_')
    ax_obj.plot(xarr, fithigh, ls='-.', c='k', lw=0.7, alpha=0.3, label='_nolegend_')
    ax_obj.fill_between(xarr, fitlow, fithigh, facecolor=fcolor, alpha=0.3)


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
    color_mag = mag_df.set_index('Phase')[color].dropna()

    ax_obj.plot(color_mag.index.values, color_mag.values, marker=sndata_df.loc[name, 'Marker'], c='k', ms=8,
                lw=0.9, ls='--', markerfacecolor='None', markeredgewidth=1, alpha=0.6, label='_nolegend_')
    ax_obj.plot(color_mag.index.values, color_mag.values, marker=sndata_df.loc[name, 'Marker'], ms=8, ls='',
                alpha=0.7, label=name)


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
    snGalmag_df = calc_extcorcolorframe(unorgframe_to_orgframe(data_df, column='FMAG'), ebv=EBVMW_mag)
    snGalerr_df = calc_extcorcolorframe(unorgframe_to_orgframe(data_df, column='FERR'), err=True, ebverr=EBVMW_err)
    snGalmag_df = snGalmag_df[snGalmag_df['Phase'] < clip_epoch].set_index('Phase')
    snGalerr_df = snGalerr_df[snGalerr_df['Phase'] < clip_epoch].set_index('Phase')

    x = inpmag_df[color].dropna().index.values
    y = inpmag_df[color].dropna().values
    yerr = inperr_df[color].dropna().values
    yGal = snGalmag_df[color].dropna().values
    yGalerr = snGalerr_df[color].dropna().values

    number_segments = 3
    myPWLF = pwlf.PiecewiseLinFit(x, y, sorted_data=True)
    myPWLF.fit(number_segments)

    xarr = np.arange(min(x), max(x), 0.5)
    xbreaks = [np.min(x)] + dict_colors[color][0] + [np.max(x)]
    fit = myPWLF.fit_with_breaks(xbreaks)
    yfit = myPWLF.predict(xarr)
    fitsigma = np.sqrt(myPWLF.prediction_variance(xarr))
    yticks = ax_obj.get_yticks(minor=True)

    print "{0} Color".format(color)
    print "Number of Parameters: {0}".format(myPWLF.n_parameters)
    print "Manual Fit: ", myPWLF.slopes * 100
    print redchisq(y, myPWLF.predict(x), n=myPWLF.n_parameters)

    for xbreak in dict_colors[color][0]:
        ax_obj.axvline(xbreak, ls='--', lw=0.8, c='k')
        if color in ['B-V', 'U-B']:
            ax_obj.text(xbreak + 1, yticks[6], r'$\sim${:.0f} d'.format(xbreak), rotation=90, color='r', fontsize=13)
        else:
            ax_obj.text(xbreak + 1, yticks[-4], r'$\sim${:.0f} d'.format(xbreak), rotation=90, color='r', fontsize=13)

    s1 = (myPWLF.slopes)[0]
    s2 = (myPWLF.slopes)[1]

    ax_obj.text((xbreaks[0] + xbreaks[1]) * 0.17, yticks[2], '$s_1$={:.2f}'.format(s1 * 100), color='b', fontsize=13)
    ax_obj.text((xbreaks[1] + xbreaks[2]) * 0.45, yticks[2], '$s_2$={:.2f}'.format(s2 * 100), color='b', fontsize=13)
    ax_obj.annotate(r'$\rm (' + color + ')_0$', xy=(1, 0), xycoords='axes fraction', fontsize=16, xytext=(-10, 10),
                    textcoords='offset points', ha='right', va='bottom')

    coeff = ls.legfit(x, yGal, 4)
    ax_obj.plot(x, ls.legval(x, coeff), lw=2, ls='--', c='r', label='_nolegend_')
    ax_obj.plot(x, y, ls='', lw=1, marker='o', c='k', ms=9, markerfacecolor='dimgrey', markeredgewidth=2,
                alpha=0.7, label=name_SN)
    ax_obj.plot(x, y, ls='', marker='o', c='k', ms=4, markerfacecolor='None', markeredgewidth=1,
                alpha=0.8, label='_nolegend_')
    ax_obj.fill_between(xarr, yfit - (fitsigma * 3), yfit + (fitsigma * 3), alpha=0.3, color='red')
    ax_obj.plot(xarr, yfit, ls='-', lw=3, c='blue', label='_nolegend_')

    data_legfit = pd.DataFrame(index=xarr)
    data_legfit.index.name = 'Phase'
    data_legfit['Color'] = np.round(ls.legval(xarr, coeff), 3)
    data_legfit['ColorErr'] = np.round(fitsigma, 3)
    data_legfit.to_csv('OUTPUT_InterpColor{0}'.format(color), sep=' ', index=True, header=True)


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
    ax_obj.yaxis.set_major_locator(MultipleLocator(0.5))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(5))
    ax_obj.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=16)
    ax_obj.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=16)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Other Type II SNe Data From The SN Archive Folder
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', DIR_SNe + 'LC_Data/*.asc', exceptions='SWIFT')
list_uvfiles = group_similar_files('', DIR_SNe + 'LC_Data/*SWIFT.asc')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Tabulate The Photometric Data Epoch-Wise
# Get The Data Ready For The Plots
# ------------------------------------------------------------------------------------------------------------------- #
list_colors = ['U-B', 'B-V', 'V-R', 'V-I']
dict_colors = {'U-B': [[39, 89], 1.0, 0.1, 15.2, 19.2], 'B-V': [[38, 92], 0.1, 0.01, 15.92, 16.38],
               'V-R': [[39, 92], 0.1, 0.01, 15.6, 16.14], 'V-I': [[38, 92], 0.1, 0.01, 15.33, 15.97]}

data_df = pd.read_csv('OUTPUT_FinalSNMagTemp', sep='\s+', engine='python')
data_df['Date'] = data_df['JD'].apply(lambda x: jd_to_cald(x))
data_df['Phase'] = (data_df['JD'] - date_explosion).round(int(precision))
data_df = data_df.set_index('Date')
data_df = data_df[data_df.index != '2017-01-10']

snmag_df = calc_extcorcolorframe(unorgframe_to_orgframe(data_df, column='FMAG'))
snerr_df = calc_extcorcolorframe(unorgframe_to_orgframe(data_df, column='FERR'), err=True)

snmag_df = snmag_df[snmag_df['Phase'] < clip_epoch].set_index('Phase')
snerr_df = snerr_df[snerr_df['Phase'] < clip_epoch].set_index('Phase')

mag_df, err_df = calc_swiftcolordf(DIR_PHOT + '2016gfy_SWIFT.dat', name=name_SN)
uvw2w1 = mag_df[['Phase', 'uvw2-uvw1']].dropna(how='any').sort_values('Phase')
uvw2v = mag_df[['Phase', 'uvw2-v']].dropna(how='any').sort_values('Phase')
uvw2w1Err = err_df[['Phase', 'uvw2-uvw1']].dropna(how='any').sort_values('Phase')
uvw2vErr = err_df[['Phase', 'uvw2-v']].dropna(how='any').sort_values('Phase')
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plot Color Terms From UBVRI Photometric Magnitudes
# # ------------------------------------------------------------------------------------------------------------------- #

# sns.set_palette(sns.color_palette('Dark2', 10)[5:9])
# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12), sharex=True)

# for file_name in list_files:
#     name = file_name.split('/')[-1].split('.')[0]
#     if name in ['2005cs', '1999em', '2004et', '2013ab']:
#         datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
#         datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
#         datacomp_df = datacomp_df[datacomp_df['Phase'] < clip_epoch]
#         plot_color(ax1, datacomp_df, name, 'U-B')
#         plot_color(ax2, datacomp_df, name, 'B-V')
#         plot_color(ax3, datacomp_df, name, 'V-R')
#         plot_color(ax4, datacomp_df, name, 'V-I')

# ax1.set_ylim(-1.3, 1.9)
# ax2.set_ylim(-1.3, 1.9)

# ax3.set_ylim(-0.15, 1.4)
# ax4.set_ylim(-0.15, 1.4)

# set_plotparams(ax1)
# ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
# plot_sncolor(ax1, snmag_df, snerr_df, color='U-B')

# set_plotparams(ax2)
# ax2.set_yticklabels([])
# ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
# plot_sncolor(ax2, snmag_df, snerr_df, color='B-V')

# set_plotparams(ax3)
# plot_sncolor(ax3, snmag_df, snerr_df, color='V-R')

# set_plotparams(ax4)
# ax4.set_yticklabels([])
# plot_sncolor(ax4, snmag_df, snerr_df, color='V-I')

# ax4.set_xlim(-10, clip_epoch + 5)
# ax3.legend(fontsize=14, markerscale=1.7, loc=2, frameon=False)
# ax3.set_xlabel('Time Since Explosion [Days]', fontsize=16)
# ax4.set_xlabel('Time Since Explosion [Days]', fontsize=16)

# fig.subplots_adjust(hspace=0.01, wspace=0.01)
# fig.savefig('PLOT_ColorEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Color Terms From UBVRI and Swift-UVOT Photometric Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('Dark2', 10)[5:9])
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 13))

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['2005cs', '1999em', '2004et', '2013ab']:
        datacomp_df = pd.read_csv(file_name, sep='\s+', comment='#')
        datacomp_df = datacomp_df.replace('INDEF', np.nan).sort_values(by='Phase')
        datacomp_df = datacomp_df[datacomp_df['Phase'] < clip_epoch]
        plot_color(ax1, datacomp_df, name, 'U-B')
        plot_color(ax3, datacomp_df, name, 'V-R')

for file_name in list_uvfiles:
    name = file_name.split('/')[-1].split('_')[0]
    if name in ['2016X', '2014cx', '2013ab', '2005cs']:
        magcomp_df, _ = calc_swiftcolordf(file_name, name=name)
        coluvw2w1 = magcomp_df[['Phase', 'uvw2-uvw1']].dropna(how='any').sort_values('Phase')
        coluvw2v = magcomp_df[['Phase', 'uvw2-v']].dropna(how='any').sort_values('Phase')
        ax4.plot(coluvw2w1['Phase'], coluvw2w1['uvw2-uvw1'], ls=':', lw=1,
                 marker=sndata_df.loc[name, 'Marker'], ms=12, alpha=0.8, label=name)
        ax4.plot(coluvw2w1['Phase'], coluvw2w1['uvw2-uvw1'], marker=sndata_df.loc[name, 'Marker'], ms=12, c='k',
                 markerfacecolor='None', markeredgewidth=0.8, ls='', alpha=0.7, label='_nolegend_')
        ax2.plot(coluvw2v['Phase'], coluvw2v['uvw2-v'], ls=':', lw=1,
                 marker=sndata_df.loc[name, 'Marker'], ms=12, alpha=0.8, label=name)
        ax2.plot(coluvw2v['Phase'], coluvw2v['uvw2-v'], marker=sndata_df.loc[name, 'Marker'], ms=12, c='k',
                 markerfacecolor='None', markeredgewidth=0.8, ls='', alpha=0.7, label='_nolegend_')

ax1.set_ylim(-1.3, 1.9)
ax1.set_xlim(-5, clip_epoch + 10)
ax3.set_ylim(-0.1, 1.0)
ax3.set_xlim(-5, clip_epoch + 10)

ax4.set_ylim(-0.45, 1.3)
ax4.set_xlim(2, 21)
ax2.set_ylim(-2.7, 2.5)
ax2.set_xlim(2, 21)

set_plotparams(ax1)
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
plot_sncolor(ax1, snmag_df, snerr_df, color='U-B')
ax1.set_title('HCT-HFOSC', fontsize=14)

set_plotparams(ax3)
ax3.yaxis.set_major_locator(MultipleLocator(0.3))
ax3.yaxis.set_minor_locator(MultipleLocator(0.03))
plot_sncolor(ax3, snmag_df, snerr_df, color='V-R')

set_plotparams(ax2)
ax2.set_title('SWIFT-UVOT', fontsize=14)
ax2.yaxis.set_major_locator(MultipleLocator(1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
ax2.xaxis.set_major_locator(MultipleLocator(4))
ax2.xaxis.set_minor_locator(MultipleLocator(0.4))

ax2.plot(uvw2v['Phase'], uvw2v['uvw2-v'], ls='', marker='o', c='k', ms=10,
         markerfacecolor='dimgrey', markeredgewidth=2, alpha=0.7, label=name_SN)
ax2.plot(uvw2v['Phase'], uvw2v['uvw2-v'], ls='-', lw=2, marker='o', c='k', ms=4,
         markerfacecolor='None', markeredgewidth=1, alpha=0.8, label='_nolegend_')
ax2.fill_between(uvw2v['Phase'], uvw2v['uvw2-v'] - uvw2vErr['uvw2-v'],
                 uvw2v['uvw2-v'] + uvw2vErr['uvw2-v'], color='dodgerblue', alpha=0.4)

set_plotparams(ax4)
ax4.yaxis.set_major_locator(MultipleLocator(0.4))
ax4.yaxis.set_minor_locator(MultipleLocator(0.04))
ax4.xaxis.set_major_locator(MultipleLocator(4))
ax4.xaxis.set_minor_locator(MultipleLocator(0.4))

# ax4.axvline(3.2, ls='--', lw=1, color='orangered')
# ax4.text(2.8, 0.0, 'Discovery Epoch', rotation=90, fontsize=12)
ax4.axvline(8, ls='--', lw=1, color='dimgrey')
ax4.axvline(25, ls='--', lw=1, color='dimgrey')
ax4.axvspan(8, 25, color='grey', alpha=0.2)
ax4.plot(uvw2w1['Phase'], uvw2w1['uvw2-uvw1'], ls='', marker='o', c='k', ms=10,
         markerfacecolor='dimgrey', markeredgewidth=2, alpha=0.7, label=name_SN)
ax4.plot(uvw2w1['Phase'], uvw2w1['uvw2-uvw1'], ls='-', lw=2, marker='o', c='k', ms=4,
         markerfacecolor='None', markeredgewidth=1, alpha=0.8, label='_nolegend_')
ax4.fill_between(uvw2w1['Phase'], uvw2w1['uvw2-uvw1'] - uvw2w1Err['uvw2-uvw1'],
                 uvw2w1['uvw2-uvw1'] + uvw2w1Err['uvw2-uvw1'], color='dodgerblue', alpha=0.4)

ax4.annotate(r'$\rm (uvw2-uvw1)_0$', xy=(1, 0), xycoords='axes fraction', fontsize=16, xytext=(-10, 10),
             textcoords='offset points', ha='right', va='bottom')
ax2.annotate(r'$\rm (uvw2-v)_0$', xy=(1, 0), xycoords='axes fraction', fontsize=16, xytext=(-10, 10),
             textcoords='offset points', ha='right', va='bottom')
ax1.legend(fontsize=14, markerscale=1.7, loc=1, frameon=False)
ax2.legend(fontsize=14, markerscale=1.3, loc=2, frameon=False)
ax3.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax4.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax4.text(11, 0.0, s='Ejecta-CSM Interaction', color='navy', alpha=0.7, fontsize=16)

fig.subplots_adjust(hspace=0.01, wspace=0.1)
fig.savefig('PLOT_OpticalColor.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
