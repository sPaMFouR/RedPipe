#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx--------------PLOT THE SUPERNOVA LIGHT CURVES FROM INPUT MAGNITUDE FILES--------------xxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import math
import pwlf
import numpy as np
import pandas as pd
import seaborn as sns
from datetime import date
import matplotlib.pyplot as plt
from jdcal import jd2gcal, gcal2jd
from matplotlib.ticker import MultipleLocator
from scipy.optimize import curve_fit, minimize
from scipy.ndimage.filters import median_filter
from astropy.modeling.blackbody import blackbody_lambda
from scipy.interpolate import CubicSpline, UnivariateSpline
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
fmt_flt = '{0:>7.3f}'
fmt_exp = '{0:>7.4e}'
epoch_mjd = 2400000.5
solar_rad = 6.957e10
Rv = 3.1
wave_data = np.linspace(3100, 9200, 1000)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
file_photvel = 'OUTPUT_PhotVel'
core_file = 'core.out'
shell_file = 'shell.out'
comb_file = 'comb.out'
obsbol_file = 'OUTPUT_DateWiseSNBolFlux'
comb_cols = ['Phase', 'LumCore', 'LumShell']

DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
DIR_PHOT = "/home/avinash/Supernovae_Data/2016gfy/Photometry/"
DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
DIR_Model = "/home/avinash/Dropbox/ModelSNe/Nagy_LC2/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of The SNe In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
EBV_mag = 0.21
EBV_err = 0.05
dist_val = 29.64
dist_err = 2.65
distmod_mag = 32.36
distmod_err = 0.18
redshift = 0.008059
JD_offset = 2457600
phase_nebstart = 115
date_explosion = 2457641.40
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS and Other Type II SNe
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
list_filters = filter_df.index.tolist()

for index, row in filter_df.iterrows():
    if len(index) == 3 and index[0:2] == 'uv':
        name = index[-1].upper()
    else:
        name = index
    if row['Offset'] > 0:
        filter_df.loc[index, 'Label'] = name + ' + ' + str(row['Offset'])
    elif row['Offset'] == 0:
        filter_df.loc[index, 'Label'] = name
    else:
        filter_df.loc[index, 'Label'] = name + ' - ' + str(abs(row['Offset']))

sndata_df = pd.read_csv(DIR_SNe + 'LC_Data/TypeIISNe.dat', sep='\s+', comment='#')
sndata_df = sndata_df.replace('INDEF', np.nan).set_index(['Name', 'Marker', 'Color']).astype('float64')
sndata_df = sndata_df.reset_index().set_index('Name')
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
# Useful Functions For Calculating Luminosity, Absolute Magnitudes & Fluxes
# ------------------------------------------------------------------------------------------------------------------- #

def blackbody_flux(wave, amp, temp):
    """
    Calculates blackbody flux as a function of wavelength (um) and temperature (K).
    Args:
        wave    : Wavelength (In Angstroms)
        amp     : Amplitude of the blackbody flux
        temp    : Temperature (In Kelvin)
    Returns:
         units of erg/s/cm^2/Angstrom
    """
    return amp * blackbody_lambda(in_x=np.asarray(wave), temperature=temp).value


def calc_radius(lum, temp):
        return (lum / (4 * np.pi * 5.67e-5 * (temp ** 4))) ** 0.5


def calc_objlum(flux, name=name_SN):
    if name != name_SN:
        dval = sndata_df.loc[name, 'D']
        derr = sndata_df.loc[name, 'DErr']
    else:
        dval = dist_val
        derr = dist_err

    val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
    lum = fmt_exp.format(val * dval ** 2)
    lumerr = fmt_exp.format(val * ((dval + derr) ** 2 - dval ** 2))

    return float(lum), float(lumerr)


def get_nickelmass(lum, phase):
    hamuy_ni = 7.866e-44 * lum * np.exp((phase * (1 + redshift) - 6.1) / 111.26)
    jerk_ni = (lum * 0.07 / 9.92e41) / (math.exp(-phase / 111.4) - math.exp(-phase / 8.8))
    return hamuy_ni, jerk_ni


def calc_magflux(mag, err, band, name=name_SN):
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

    flux = float(fmt_exp.format(10 ** (-0.4 * (mag - rlambda * ebvmag + zp + 21.100))))
    fluxerr = fmt_exp.format(abs(flux - 10 ** (-0.4 * (mag + err - rlambda * ebvmag + zp + 21.100))))

    return float(absmag), float(abserr), float(flux), float(fluxerr)

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
    time_tuple = jd2gcal(epoch_mjd, julian_day - epoch_mjd)
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
# Functions For Manipulating Pandas DataFrame Containing Data Of Other Well Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def multicol_to_fluxdf(name, input_df):
    """
    Converts a column-wise magnitude Pandas DataFrame to a row-wise Pandas DataFrame with
    the flux values and the absolute magnitudes.
    Args:
        name        : Name of the SN whose data is read
        input_df    : Input Pandas DataFrame
    Returns:
        output_df   : Output Pandas DataFrame
    """
    input_df = input_df.set_index('JD')
    for column in ['Date', 'Phase']:
        if column in input_df.columns.values:
            input_df = input_df.drop(column, axis=1)
    data_arr = input_df.values

    size = data_arr.shape
    list_jd = np.repeat(input_df.index.values, (size[1] / 2))
    list_filters = [x for x in input_df.columns.values if 'Err' not in x]
    data_arr = np.reshape(data_arr, [size[0] * size[1] / 2, 2])

    output_df = pd.DataFrame(data_arr, index=list_jd, columns=['FMAG', 'FERR'])
    output_df.index.name = 'JD'
    output_df = output_df.reset_index(drop=False)

    output_df['FILTER'] = list_filters * size[0]
    output_df['Date'] = output_df['JD'].apply(jd_to_cald)
    output_df['Phase'] = output_df['JD'] - sndata_df.loc[name, 'DateExp']

    output_df = output_df.replace('INDEF', np.nan).dropna(axis=0, how='any')
    output_df = output_df[['Date', 'JD', 'Phase', 'FILTER', 'FMAG', 'FERR']].reset_index(drop=True)
    output_df['ALambda'] = output_df['FILTER'].apply(lambda x: float(fmt_flt.format(
        filter_df.loc[x, 'RLambda'] * sndata_df.loc[name, 'EBV'])))

    for index, band in output_df['FILTER'].items():
        magflux = calc_magflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'],
                               band=band, name=name)
        output_df.loc[index, 'AbsMag'] = magflux[0]
        output_df.loc[index, 'AbsErr'] = magflux[1]
        output_df.loc[index, 'Flux'] = magflux[2]
        output_df.loc[index, 'FluxErr'] = magflux[3]

    return output_df


def calc_boldf(name, input_df, flux='Flux', fluxerr='FluxErr', plot=False):
    """
    Creates a Pandas DataFrame with Bolometric Fluxes.
    Args:
        name        : Name of the SNe whose bolometric magnitudes are to be computed
        input_df    : Input Pandas DataFrame containing individual band fluxes
        flux        : Name of the Flux column in the Pandas DataFrame
        fluxerr     : Name of the Flux Error column in the Pandas DataFrame
        plot        : Whether the spline fits to the fluxes should be plotted
    Returns:
        output_df   : Output Pandas DataFrame containing bolometric fluxes
    """
    input_df = input_df.set_index('JD')

    dict_val = {}
    for index, row in input_df.iterrows():
        if index not in dict_val.keys():
            dict_val[index] = {}
        if row['FILTER'] not in dict_val[index]:
            dict_val[index][row['FILTER']] = []
            dict_val[index][row['FILTER'] + 'Err'] = []

        dict_val[index][row['FILTER']].append(row[flux])
        dict_val[index][row['FILTER'] + 'Err'].append(row[fluxerr])

    for (day, dict_date) in dict_val.items():
        for (band, list_flux) in dict_date.items():
            if 'Err' not in band:
                dict_val[day][band] = np.mean(list_flux)
            else:
                dict_val[day][band] = np.sqrt(np.sum(np.square(list_flux)))

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 2:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if 'Err' not in band:
                    dict_flux[day][filter_df.loc[band, 'CentreWave']] = flux
                else:
                    dict_flux[day][str(filter_df.loc[band[:-3], 'CentreWave']) + 'Err'] = flux

    flux_df = pd.DataFrame(dict_flux).T
    flux_df.index.name = 'JD'
    flux_df = flux_df.interpolate(method='polynomial', order=1, limit=3).T

    dict_bolflux = {}
    for jd in flux_df.columns.values:
        series = flux_df[jd].dropna()
        if len(series) < 3:
            continue

        mag = series.loc[[x for x in list(series.index) if type(x) != str]]
        err = series.loc[[x for x in list(series.index) if type(x) == str]]

        spline = CubicSpline(mag.index.values.tolist(), mag.values.tolist(), bc_type='natural', extrapolate=True)

        if jd - sndata_df.loc[name, 'DateExp'] > 99:
            wave_data = np.linspace(4000, 9200, 1000)
        else:
            wave_data = np.linspace(3100, 9200, 1000)

        flux_data = spline(wave_data)
        flux_data[flux_data < 0] = 0
        netflux = np.trapz(flux_data, wave_data)

        dict_bolflux[jd] = {}
        dict_bolflux[jd]['Date'] = jd_to_cald(jd)
        dict_bolflux[jd]['Phase'] = jd - sndata_df.loc[name, 'DateExp']
        dict_bolflux[jd]['Flux'] = netflux
        dict_bolflux[jd]['Lum'] = calc_objlum(netflux, name=name)[0]
        dict_bolflux[jd]['LumErr'] = calc_objlum(netflux, name=name)[1]

        if plot:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wave_data, spline(wave_data), 'r-', label='CubicSpline Fit')
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig)

    if len(dict_bolflux) != 0:
        output_df = pd.DataFrame(dict_bolflux).T
        output_df.index.name = 'JD'
        output_df = output_df.reset_index().set_index(keys='Date', drop=True)
    else:
        output_df = pd.DataFrame()

    return output_df


def calc_swiftuvdf(file_name, name):
    """
    Creates a Pandas DataFrame with SWIFT UVOT magnitudes arranged epoch-wise from a file 'file_name'.
    Args:
        file_name   : Text file containing SWIFT UVOT magnitudes
        name        : Name of the supernova for which SWIFT UVOT magntiudes are to be plotted
    Returns:
        output_df   : Pandas DataFrame containing organised SWIFT UVOT magnitudes
    """
    input_df = pd.read_csv(file_name, sep='\s+', comment='#', usecols=[0, 1, 2, 3], header=None)
    input_df = input_df.rename(columns={0: 'FILTER', 1: 'MJD', 2: 'FMAG', 3: 'FERR'}).replace('INDEF', np.nan)
    input_df['FILTER'] = input_df['FILTER'].apply(lambda x: 'uv' + x.lower() if 'uv' not in x.lower() else x.lower())

    output_df = input_df.set_index('FILTER').astype('float64').reset_index().copy()
    output_df['JD'] = (output_df['MJD'] + epoch_mjd).round(1)
    output_df['Date'] = output_df['JD'].apply(jd_to_cald)
    output_df['Phase'] = (output_df['JD'] - sndata_df.loc[name, 'DateExp']).round(1)

    for index, band in output_df['FILTER'].iteritems():
        data_magflux = calc_magflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'],
                                    band=band, name=name)
        output_df.loc[index, 'AbsMag'] = data_magflux[0]
        output_df.loc[index, 'AbsErr'] = data_magflux[1]
        output_df.loc[index, 'Flux'] = data_magflux[2]
        output_df.loc[index, 'FluxErr'] = data_magflux[3]
        output_df.loc[index, 'ALambda'] = float(fmt_flt.format(filter_df.loc[band, 'RLambda'] * EBV_mag))

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Manipulating Pandas DataFrames Containing Data From SN In Study
# ------------------------------------------------------------------------------------------------------------------- #

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


def calc_snboldf(input_df, mag='FMAG', magerr='FERR', flux='Flux', fluxerr='FluxErr', tempsub=False, plot=False):
    """
    Creates a Pandas DataFrame with bolometric luminosity arranged epoch-wise from a DataFrame
    with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes
        mag         : Name of the Magnitude column in the Pandas DataFrame
        magerr      : Name of tge Magnitude Error column in the Pandas DataFrame
        flux        : Name of the Flux column in the Pandas DataFrame
        fluxerr     : Name of the Flux Error column in the Pandas DataFrame
        tempsub     : Whether the magnitudes are template subtracted?
        plot        : Whether the spline fits to the fluxes should be plotted
    Returns:
        bolm_df     : Pandas DataFrame containing organised
        flux_df     : Pandas DataFrame containing organised broadband fluxes
    """
    if not tempsub:
        add_text = ''
    else:
        add_text = 'Temp'

    input_df = input_df.set_index('JD')

    dict_mag = {}
    for index, row in input_df.iterrows():
        if index not in dict_mag.keys():
            dict_mag[index] = {}
        if row['FILTER'] not in dict_mag[index]:
            dict_mag[index][row['FILTER']] = []
            dict_mag[index][row['FILTER'] + 'Err'] = []

        dict_mag[index][row['FILTER']].append(row['FMAG'])
        dict_mag[index][row['FILTER'] + 'Err'].append(row['FERR'])

    for (day, dict_date) in dict_mag.items():
        for (band, list_mag) in dict_date.items():
            dict_mag[day][band] = float(np.mean(list_mag))

    mag_df = pd.DataFrame(dict_mag).T
    mag_df.index.name = 'JD'
    mag_df = mag_df.reset_index()
    mag_df['Date'] = mag_df['JD'].apply(jd_to_cald)
    mag_df['Phase'] = mag_df['JD'].apply(lambda x: x - date_explosion).round(int(precision))
    mag_df.to_csv('OUTPUT_DateWiseSNAppMag' + add_text, sep=' ', index=True, na_rep='INDEF')

    dict_val = {}
    for index, row in input_df.iterrows():
        if index not in dict_val.keys():
            dict_val[index] = {}
        if row['FILTER'] not in dict_val[index]:
            dict_val[index][row['FILTER']] = []
            dict_val[index][row['FILTER'] + 'Err'] = []

        dict_val[index][row['FILTER']].append(row[flux])
        dict_val[index][row['FILTER'] + 'Err'].append(row[fluxerr])

    for (day, dict_date) in dict_val.items():
        for (band, list_flux) in dict_date.items():
            if 'Err' not in band:
                dict_val[day][band] = np.mean(list_flux)
            else:
                dict_val[day][band] = np.sqrt(np.sum(np.square(list_flux)))

    flux_df = pd.DataFrame(dict_val).T
    flux_df.index.name = 'JD'
    flux_df = flux_df.reset_index()
    flux_df['Date'] = flux_df['JD'].apply(jd_to_cald)
    flux_df['Phase'] = flux_df['JD'].apply(lambda x: x - date_explosion).round(int(precision))
    flux_df.to_csv('OUTPUT_DateWiseSNAppFlux' + add_text, sep=' ', index=False, na_rep='INDEF')

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 3:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if 'Err' not in band:
                    dict_flux[day][filter_df.loc[band, 'CentreWave']] = flux
                else:
                    dict_flux[day][str(filter_df.loc[band.rstrip('Err'), 'CentreWave']) + 'Err'] = flux

    bolm_df = pd.DataFrame(dict_flux).T
    bolm_df.index.name = 'JD'
    bolm_df = bolm_df.interpolate(method='polynomial', order=1, limit=3).T

    dict_bolflux = {}
    for caljd in bolm_df.columns.values:
        series = bolm_df[caljd].dropna()
        mag = series.loc[[x for x in list(series.index) if type(x) != str]]
        err = series.loc[[x for x in list(series.index) if type(x) == str]]
        spline = CubicSpline(mag.index.values.tolist(), mag.values.tolist(), bc_type='natural', extrapolate=True)

        if caljd - date_explosion > 99:
            wave_data = np.linspace(4000, 9200, 1000)
        else:
            wave_data = np.linspace(3100, 9200, 1000)

        flux_data = spline(wave_data)
        flux_data[flux_data < 0] = 0
        netflux = np.trapz(flux_data, wave_data)

        popt, pcov = curve_fit(blackbody_flux, mag.index.values.tolist(), mag.values.tolist(),
                               sigma=err.values.tolist(), p0=[1e-15, 10000])
#         print ("Best-Fit Temperature[{0:.1f}]= {1:.2f}+/-{2:.2f}".format(caljd - date_explosion, popt[1],
#                                                                          np.sqrt(np.diag(pcov)[1])))

        lum, lumerr = calc_objlum(netflux, name=name_SN)
        radius = calc_radius(lum, popt[1])
        dict_bolflux[caljd] = {}
        dict_bolflux[caljd]['Date'] = jd_to_cald(caljd)
        dict_bolflux[caljd]['Phase'] = caljd - date_explosion
        dict_bolflux[caljd]['Flux'] = netflux
        dict_bolflux[caljd]['Lum'] = lum
        dict_bolflux[caljd]['LumErr'] = lumerr
        dict_bolflux[caljd]['Temp'] = popt[1]
        dict_bolflux[caljd]['TempErr'] = np.sqrt(np.diag(pcov)[1])
        dict_bolflux[caljd]['Rad'] = radius
        dict_bolflux[caljd]['RadErr'] = (radius / 2) * ((lumerr / lum) - (4 * np.sqrt(np.diag(pcov)[1]) / popt[1]))

        if plot:
            fig_temp = plt.figure(figsize=(10, 8))
            ax = fig_temp.add_subplot(111)
            wave2_data = np.linspace(1600, 9200, 1000)
            ax.plot(mag.index.values, mag.values, 'o', label='Data Points')
            ax.plot(wave_data, flux_data, 'r-', label='CubicSpline Fit')

            if len(mag.index.values) > 5:
                ax.plot(wave2_data, blackbody_flux(wave2_data, *popt), 'k--', label='Blackbody Fit')
            else:
                ax.plot(wave_data, blackbody_flux(wave_data, *popt), 'k--', label='Blackbody Fit')

            ax.set_ylabel(r'Apparent Flux [$\rm erg\ s^{-1}\ cm^{-2}\ \AA^{-1}$]')
            ax.set_xlabel(r'Wavelength [$\rm \AA$]')
            ax.set_title("Temperature = {0:.2f}+/- {1:.2f}".format(popt[1], np.sqrt(np.diag(pcov)[1])))
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig_temp)

    fbolm_df = pd.DataFrame(dict_bolflux).T
    fbolm_df.index.name = 'JD'
    fbolm_df.to_csv('OUTPUT_DateWiseSNBolFlux' + add_text, sep=' ', index=True)

    for index, row in fbolm_df.iterrows():
        if row['Phase'] > 150:
            fbolm_df.loc[index, 'MNi'] = get_nickelmass(row['Lum'], phase=row['Phase'])[0]
            fbolm_df.loc[index, 'MNi2'] = get_nickelmass(row['Lum'], phase=row['Phase'])[1]
            fbolm_df.loc[index, 'MNiErr'] = get_nickelmass(row['LumErr'], phase=row['Phase'])[0]

    return bolm_df, fbolm_df


def calc_snfluxdf(file_name, concat_df=pd.DataFrame(), tempsub=False):
    """
    Creates a Pandas DataFrame with Individual Band Fluxes.
    Args:
        file_name   : Name of the file which has the photometric magnitudes for the SNe
        concat_df   : Input Pandas DataFrame to be appended to the original data
        tempsub     : Boolean describing whether these magnitudes are template subtracted?
    Returns:
        output_df   : Output Pandas DataFrame containing individual band fluxes
    """
    if not tempsub:
        add_text = ''
    else:
        add_text = 'Temp'

    output_df = pd.read_csv(file_name, sep="\s+", engine='python')
    if not concat_df.empty:
        output_df = pd.concat([concat_df, output_df], axis=0)

    output_df = output_df.sort_values(by=['FILTER', 'JD'], kind='mergesort')
    output_df['Date'] = output_df['JD'].apply(lambda x: jd_to_cald(x))
    output_df['Phase'] = output_df['JD'].apply(lambda x: x - date_explosion).round(int(precision))

    output_df = output_df[['Date', 'JD', 'Phase', 'FILTER', 'FMAG', 'FERR']].reset_index(drop=True)
    output_df.to_csv('OUTPUT_NetSNMag' + add_text, sep=' ', index=False)

    for index, band in output_df['FILTER'].items():
        data_magflux = calc_magflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'], band=band)
        output_df.loc[index, 'AbsMag'] = data_magflux[0]
        output_df.loc[index, 'AbsErr'] = data_magflux[1]
        output_df.loc[index, 'Flux'] = data_magflux[2]
        output_df.loc[index, 'FluxErr'] = data_magflux[3]
        output_df.loc[index, 'ALambda'] = float(fmt_flt.format(filter_df.loc[band, 'RLambda'] * EBV_mag))

    output_df.to_csv('OUTPUT_NetSNFlux' + add_text, sep=' ', index=False)

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Setting Plot Parameters And Plotting
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj):
    """
    Sets plot parameters to the axes object 'ax_obj'.
    Args:
        ax_obj  : Axes object to be used for plotting and setting plot parameters
    Returns:
        None
    """
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(100))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
    ax_obj.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)
    ax_obj.set_xlabel('Time Since Explosion [Days]', fontsize=16)


def plot_uvband(ax_obj, name, data_df, band, label=True):
    """
    Plot the SWIFT UVOT band to the axes object 'ax_obj'.
    Args:
        ax_obj  : Axes object to be used for plotting and setting plot parameters
        name    : Name of Supernova whose magnitudes are to be plotted
        data_df : Pandas DataFrame containing SWIFT UVOT magnitudes
        band    : SWIFT UVOT band to be plotted
        label   : Boolean stating whether the plot is to be labelled
    Returns:
        None
    """
    plot_df = data_df[data_df['FILTER'] == band].copy()

    if label:
        ax_obj.plot(plot_df['Phase'], plot_df['AbsMag'], ls=':', markersize=8, c=sndata_df.loc[name, 'Color'],
                    marker=sndata_df.loc[name, 'Marker'], label=name)
    else:
        ax_obj.plot(plot_df['Phase'], plot_df['AbsMag'], ls=':', markersize=8, c=sndata_df.loc[name, 'Color'],
                    marker=sndata_df.loc[name, 'Marker'], label='_nolegend_')

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Other SNe Data From Archive Folder
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', DIR_SNe + 'LC_Data/*.asc', exceptions='SWIFT')
list_uvfiles = group_similar_files('', DIR_SNe + 'LC_Data/*SWIFT*.asc')
list_modnifiles = group_similar_files('', DIR_SNe + 'NiMix_Data/*.dat')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
rawopt_df = calc_snfluxdf('OUTPUT_FinalSNMag', tempsub=False)
outputopt_df = calc_snfluxdf('OUTPUT_FinalSNMagTemp', tempsub=True)

max_epoch = outputopt_df['Phase'].max()
outputopt_df = outputopt_df[~((outputopt_df['Date'] == '2016-12-29') & (outputopt_df['FILTER'] == 'U'))]
outputopt_df = outputopt_df[~((outputopt_df['Date'] == '2017-10-01') & (outputopt_df['FILTER'] == 'U'))]
outputopt_df[['JD', 'Phase']] = outputopt_df[['JD', 'Phase']].round(1)

vabs_df = outputopt_df[outputopt_df['FILTER'] == 'V'].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read SWIFT Data And Combine It With Optical Data To Calculate Bolometric Light Curve
# Read Photospheric Velocity For Determining Radius
# Read Nagy Fit To The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
outputuv_df = calc_swiftuvdf(DIR_PHOT + '2016gfy_SWIFT.dat', name=name_SN)
outputnet_df = pd.concat([outputopt_df, outputuv_df], sort=True)
outputuv_df = outputuv_df.set_index('JD').dropna(axis=0, how='any')
bolm_df, fbolm_df = calc_snboldf(outputopt_df)
# _, fbolm2_df = calc_snboldf(outputnet_df)

temprad_df = fbolm_df[fbolm_df['Phase'] <= phase_nebstart].copy()
temprad_df[['Temp', 'TempErr']] = temprad_df[['Temp', 'TempErr']] / 1000.
temprad_df[['Rad', 'RadErr']] = temprad_df[['Rad', 'RadErr']] / (1000 * solar_rad)

vel_df = pd.read_csv(DIR_SPEC + file_photvel, sep='\s+').replace('INDEF', np.nan).drop('Date', axis=1)
vel_df = vel_df.astype('float64')[['Phase', '5169', '5169Err']].dropna()
vel_df['Rad'] = vel_df['Phase'] * vel_df['5169'] * 86400e5 / (1000 * solar_rad)
vel_df['RadErr'] = vel_df['Phase'] * vel_df['5169Err'] * 86400e5 / (1000 * solar_rad)

obsbol_df = pd.read_csv(DIR_PHOT + obsbol_file, sep='\s+')
obsbol_df = obsbol_df[obsbol_df['Phase'] < 205]
obsbol_df['LogLum'] = obsbol_df['Lum'].apply(lambda x: np.log10(x))

comb_df = pd.read_csv(DIR_Model + comb_file, sep='\s+', header=None, names=comb_cols, comment='#')
comb_df['Lum'] = comb_df['LumCore'] + comb_df['LumShell']
comb_df['LogLumCore'] = comb_df['LumCore'].apply(lambda x: np.log10(x))
comb_df['LogLumShell'] = comb_df['LumShell'].apply(lambda x: np.log10(x))
comb_df['LogLum'] = comb_df['Lum'].apply(lambda x: np.log10(x))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Peacewise Polynomial To Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
x = np.array(fbolm_df['Phase'])
y = np.array(fbolm_df['Lum'].apply(lambda lum: np.log10(lum)))

number_segments = 5
myPWLF = pwlf.PiecewiseLinFit(x, y, sorted_data=True)
myPWLF.fit(number_segments)

xarr = np.linspace(np.min(x), np.max(x), num=1000)
xguess = np.zeros(number_segments - 1)
xguess[0:] = [8, 39, 91, 124]

myPWLF.fit_with_breaks([np.min(x)] + [8, 39, 91, 124] + [np.max(x)])
fit = myPWLF.predict(xarr)
fitsigma = np.sqrt(myPWLF.prediction_variance(xarr))

s1 = -(myPWLF.slopes)[1] * 100
s2 = -(myPWLF.slopes)[2] * 100
s3 = -(myPWLF.slopes)[4] * 100

print "Number of Parameters: {0}".format(myPWLF.n_parameters)
print "Best Fit: ", myPWLF.slopes * 100
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Tabulate The Photometric Data Epoch-Wise Onto An Output File
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv('OUTPUT_FinalSNMagTemp', sep='\s+')

data_df['Date'] = data_df['JD'].apply(lambda x: jd_to_cald(x))
data_df['Phase'] = (data_df['JD'] - date_explosion).round(int(precision))
data_df['Mag'] = data_df['FMAG'].apply(lambda x: '{:.2f}'.format(x)) + r'$\pm$' + data_df['FERR'].apply(lambda x: '{:.2f}'.format(x))

tabular_df = data_df[['FILTER', 'Date', 'JD', 'Phase', 'Mag']]
tabular_df.to_csv('OUTPUT_FinalTabularSNMag', sep=' ', index=False)
display_text("HCT Magnitudes For SN {0} Have Been Tabulated Epoch-Wise".format(name_SN))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Output The HCT Photometric Magnitudes Onto A Latex Table
# ------------------------------------------------------------------------------------------------------------------- #
data_df = data_df.set_index('Date')
tabopt_df = organise_sndf(data_df, column='Mag')

tabopt_df['Phase'] = tabopt_df['Date'].apply(lambda x: data_df.loc[data_df.index == x, 'Phase'].iloc[0])
tabopt_df['JD'] = np.round(tabopt_df['Phase'] + date_explosion - JD_offset, 2)
tabopt_df['Phase'] = tabopt_df['Phase'].apply(lambda x: "{:.2f}".format(x) if x < 0 else "+{:.2f}".format(x))

tabopt_df = tabopt_df[['Date', 'JD', 'Phase', 'U', 'B', 'V', 'R', 'I']].sort_values(by='JD')
tabopt_df = tabopt_df.rename(columns={'Phase': 'Phase$^*$', 'U': '$U$', 'B': '$B$', 'V': '$V$',
                                      'R': '$R$', 'I': '$I$'}).dropna(how='all', axis=1)
tabopt_df = tabopt_df.replace(np.nan, '---', regex=True)

tabopt_df.to_latex('_PhotHCT.tex', escape=False, index=False)
display_text("HCT Magnitudes For SN {0} Have Been Logged Onto A Latex File".format(name_SN))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Output The SWIFT Photometric Magnitudes Onto A Latex Table
# ------------------------------------------------------------------------------------------------------------------- #
outputuv_df['Mag'] = outputuv_df['FMAG'].apply(lambda x: "{:.2f}".format(x)) + r'$\pm$' + outputuv_df['FERR'].apply(lambda x: "{:.2f}".format(x))
outputuv_df = outputuv_df.set_index('Date', drop=False)
tabuv_df = organise_sndf(outputuv_df, column='Mag')

tabuv_df['Phase'] = tabuv_df['Date'].apply(lambda x: outputuv_df.loc[outputuv_df.index == x, 'Phase'].iloc[0])
tabuv_df['JD'] = tabuv_df['Phase'] + date_explosion - JD_offset
tabuv_df['Phase'] = tabuv_df['Phase'].apply(lambda x: "{:.2f}".format(x) if x < 0 else "+{:.2f}".format(x))

tabuv_df = tabuv_df[['Date', 'JD', 'Phase', 'uvw2', 'uvm2', 'uvw1', 'uvu', 'uvv']].sort_values(by='JD')
tabuv_df = tabuv_df.rename(columns={'Phase': 'Phase$^*$', 'uvw2': '$uvw2$', 'uvm2': '$uvm2$', 'uvw1': '$uvw1$',
                                    'uvu': '$uvu$', 'uvv': '$uvv$'})
tabuv_df = tabuv_df.replace(np.nan, '---', regex=True).dropna(how='all', axis=1)

tabuv_df.to_latex('_PhotSWIFT.tex', escape=False, index=False)
display_text("SWIFT Magnitudes For SN {0} Have Been Logged Onto A Latex File".format(name_SN))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Compute The Eeta_Ni For The SN In Study
# ------------------------------------------------------------------------------------------------------------------- #
mNi = 0.033
phase_nebbol = 124
platarr = np.arange(0, phase_nebbol, 0.1)
bollum = UnivariateSpline(fbolm_df['Phase'], fbolm_df['Lum'], k=1)


def calc_nilum(t):
    return mNi * (6.45 * np.exp(-t / 8.8) + 1.45 * np.exp(-t / 111.3)) * 1e43


def calc_eeta():
    return np.trapz(platarr * calc_nilum(platarr), platarr) / np.trapz(platarr * (bollum(platarr) - calc_nilum(platarr)), platarr)


def calc_deltae():
    return bollum(25) * (80 ** 2) / np.trapz(platarr * (bollum(platarr) - calc_nilum(platarr)), platarr)


def calc_delta():
    return bollum(25) * (80 ** 2) / np.trapz(platarr * (bollum(platarr)), platarr)


def calc_deltaM():
    return -2.5 * np.log10(bollum(75) / bollum(25))


print calc_eeta()
print calc_deltae(), 2.5 * np.log10(calc_deltae()), 0.44 + 0.82 * 2.5 * np.log10(calc_delta())
print calc_delta(), 2.5 * np.log10(calc_delta())
print calc_deltaM()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The HCT Apparent Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_app = plt.figure(figsize=(8, 8))
ax_app = fig_app.add_subplot(111)

for band, band_df in outputopt_df.groupby('FILTER'):
    ax_app.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                   marker=filter_df.loc[band, 'Marker'], s=30, alpha=0.6, label=filter_df.loc[band, 'Label'])
    ax_app.plot(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c='k', markerfacecolor='None',
                ls='', markeredgewidth=1, marker=filter_df.loc[band, 'Marker'], ms=6, alpha=0.6, label='_nolegend_')
    ax_app.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                    c=filter_df.loc[band, 'Color'], ls='--', lw=0.5, capsize=2, capthick=1, label='_nolegend_')

handles, labels = ax_app.get_legend_handles_labels()
handles = [handles[3], handles[0], handles[4], handles[2], handles[1]]
labels = [labels[3], labels[0], labels[4], labels[2], labels[1]]
ax_app.legend(handles, labels, fontsize=12, markerscale=2, loc=1, frameon=False)

set_plotparams(ax_app)
ax_app.set_ylim(22.0, 14.5)
ax_app.set_xlim(-15, max_epoch + 20)
ax_app.yaxis.set_major_locator(MultipleLocator(1))
ax_app.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_app.set_ylabel('Apparent Magnitude [mag]', fontsize=16)

fig_app.savefig('PLOT_ApparentLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_app)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The SWIFT Apparent Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_swift = plt.figure(figsize=(8, 8))
ax_swift = fig_swift.add_subplot(111)

for band, band_df in outputuv_df.groupby('FILTER'):
    ax_swift.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                     marker=filter_df.loc[band, 'Marker'], s=35, label=filter_df.loc[band, 'Label'])
    ax_swift.plot(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c='k', markerfacecolor='None',
                  ls='', markeredgewidth=1, marker=filter_df.loc[band, 'Marker'], ms=7, alpha=0.6, label='_nolegend_')
    ax_swift.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                      c=filter_df.loc[band, 'Color'], ls='--', lw=1, capsize=2, capthick=1, label='_nolegend_')

handles, labels = ax_swift.get_legend_handles_labels()
handlesuv = [handles[5], handles[1], handles[4], handles[2], handles[0], handles[3]]
labelsuv = [labels[5], labels[1], labels[4], labels[2], labels[0], labels[3]]
ax_swift.legend(handlesuv, labelsuv, fontsize=11, markerscale=2, loc=4, frameon=False)

ax_swift.set_ylim(22.0, 14.5)
ax_swift.set_xlim(-1, 31)
ax_swift.yaxis.set_ticks_position('both')
ax_swift.xaxis.set_ticks_position('both')
ax_swift.yaxis.set_major_locator(MultipleLocator(1))
ax_swift.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_swift.xaxis.set_major_locator(MultipleLocator(5))
ax_swift.xaxis.set_minor_locator(MultipleLocator(0.5))
ax_swift.set_ylabel('Apparent Magnitude [mag]', fontsize=16)
ax_swift.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_swift.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_swift.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

fig_swift.savefig('PLOT_ApparentUVLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_swift)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Apparent Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_comb, (ax_opt, ax_uv) = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

for band, band_df in outputopt_df.groupby('FILTER'):
    ax_opt.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                   marker=filter_df.loc[band, 'Marker'], s=30, alpha=0.6, label=filter_df.loc[band, 'Label'])
    ax_opt.plot(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c='k', markerfacecolor='None',
                ls='', markeredgewidth=0.7, marker=filter_df.loc[band, 'Marker'], ms=6, alpha=0.6, label='_nolegend_')
    ax_opt.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                    c=filter_df.loc[band, 'Color'], ls='--', lw=0.5, capsize=2, capthick=1, label='_nolegend_')

handles, labels = ax_opt.get_legend_handles_labels()
handles = [handles[3], handles[0], handles[4], handles[2], handles[1]]
labels = [labels[3], labels[0], labels[4], labels[2], labels[1]]
ax_opt.legend(handles, labels, fontsize=16, markerscale=2.5, frameon=False, loc=1)

set_plotparams(ax_opt)
ax_opt.set_ylim(22.0, 14.5)
ax_opt.set_xlim(-15, max_epoch + 20)
ax_opt.yaxis.set_major_locator(MultipleLocator(1))
ax_opt.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_opt.set_ylabel('Apparent Magnitude [mag]', fontsize=16)
ax_opt.set_title('HCT-HFOSC', fontsize=16)

for band, band_df in outputuv_df.groupby('FILTER'):
    ax_uv.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                  marker=filter_df.loc[band, 'Marker'], s=50, label=filter_df.loc[band, 'Label'])
    ax_uv.plot(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c='k', markerfacecolor='None',
               ls='', markeredgewidth=0.7, marker=filter_df.loc[band, 'Marker'], ms=7, alpha=0.6, label='_nolegend_')
    ax_uv.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                   c=filter_df.loc[band, 'Color'], ls='--', lw=0.5, capsize=2, capthick=1, label='_nolegend_')

handles, labels = ax_uv.get_legend_handles_labels()
handlesuv = [handles[5], handles[1], handles[4], handles[2], handles[0], handles[3]]
labelsuv = [labels[5], labels[1], labels[4], labels[2], labels[0], labels[3]]
ax_uv.legend(handlesuv, labelsuv, fontsize=14, markerscale=2, loc=4, frameon=False)

ax_uv.set_title('SWIFT-UVOT', fontsize=16)
ax_uv.set_xlim(3, 32)
# ax_uv.set_yticklabels([])
ax_uv.yaxis.set_ticks_position('both')
ax_uv.xaxis.set_ticks_position('both')
ax_uv.yaxis.set_major_locator(MultipleLocator(1))
ax_uv.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_uv.xaxis.set_major_locator(MultipleLocator(5))
ax_uv.xaxis.set_minor_locator(MultipleLocator(0.5))
ax_uv.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_uv.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_uv.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

fig_comb.subplots_adjust(wspace=0.01)
fig_comb.savefig('PLOT_ApparentCombLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_comb)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Comparison Of The Template Subtracted Magnitudes With Original Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
fig_comp = plt.figure(figsize=(9, 6))
ax_comp = fig_comp.add_subplot(111)

for band, band_df in outputopt_df.groupby('FILTER'):
    ax_comp.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'],
                    marker=filter_df.loc[band, 'Marker'],
                    c=filter_df.loc[band, 'Color'], s=10, label=filter_df.loc[band, 'Label'], alpha=0.5)
    ax_comp.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'],
                     c=filter_df.loc[band, 'Color'], fmt='', ls='', lw=0.5, capsize=2, capthick=1,
                     label='_nolegend_')

for band, band_df in rawopt_df.groupby('FILTER'):
    ax_comp.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], marker='+',
                    c='lightsalmon', s=10, label=filter_df.loc[band, 'Label'], alpha=0.5)
    ax_comp.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'],
                     c='lightsalmon', fmt='', ls='--', lw=0.5, capsize=2, capthick=1, label='_nolegend_')

set_plotparams(ax_comp)
ax_comp.set_ylim(22.0, 14.5)
ax_comp.set_xlim(-10, max_epoch + 20)
ax_comp.legend(markerscale=2, frameon=False)
ax_comp.yaxis.set_major_locator(MultipleLocator(2))
ax_comp.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_comp.set_ylabel('Apparent Magnitude [mag]', fontsize=16)

fig_comp.savefig('PLOT_TempCompLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_comp)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The SWIFT UV Absolute Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_uvabs = plt.figure(figsize=(8, 12))
ax_uvw1 = fig_uvabs.add_subplot(311)
ax_uvm2 = fig_uvabs.add_subplot(312, sharex=ax_uvw1)
ax_uvw2 = fig_uvabs.add_subplot(313, sharex=ax_uvw1)

for file_name in list_uvfiles:
    name = file_name.split('/')[-1].split('_')[0]
    if name in ['2016X', '2013ej', '2012aw', '2014cx', '2013ab']:
        data_df = calc_swiftuvdf(file_name, name)
        data_df = data_df[data_df['Phase'] < 50].sort_values(by='Phase')
        plot_uvband(ax_uvw1, name, data_df, 'uvw1')
        plot_uvband(ax_uvm2, name, data_df, 'uvm2', label=False)
        plot_uvband(ax_uvw2, name, data_df, 'uvw2', label=False)

set_plotparams(ax_uvw1)
ax_uvw1.invert_yaxis()
ax_uvw1.xaxis.set_major_locator(MultipleLocator(10))
ax_uvw1.xaxis.set_minor_locator(MultipleLocator(1))
ax_uvw1.yaxis.set_major_locator(MultipleLocator(2))
ax_uvw1.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_uvw1.set_ylabel(r'$\rm M_{uvw1}$ [mag]', fontsize=16)

set_plotparams(ax_uvm2)
ax_uvm2.invert_yaxis()
ax_uvm2.xaxis.set_major_locator(MultipleLocator(10))
ax_uvm2.xaxis.set_minor_locator(MultipleLocator(1))
ax_uvm2.yaxis.set_major_locator(MultipleLocator(2))
ax_uvm2.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_uvm2.set_ylabel(r'$\rm M_{uvm2}$ [mag]', fontsize=16)

set_plotparams(ax_uvw2)
ax_uvw2.invert_yaxis()
ax_uvw2.xaxis.set_major_locator(MultipleLocator(10))
ax_uvw2.xaxis.set_minor_locator(MultipleLocator(1))
ax_uvw2.yaxis.set_major_locator(MultipleLocator(2))
ax_uvw2.yaxis.set_minor_locator(MultipleLocator(0.2))

for band, band_df in outputuv_df.groupby('FILTER'):
    if band == 'uvw1':
        ax_uvw1.scatter(band_df['Phase'], band_df['AbsMag'], c='k', marker='*', s=180, label=name_SN)
        ax_uvw1.errorbar(band_df['Phase'], band_df['AbsMag'], yerr=band_df['AbsErr'], fmt='',
                         marker='*', c='k', ls='-', lw=1.2, capsize=3, capthick=1, label='_nolegend_')
    elif band == 'uvm2':
        ax_uvm2.scatter(band_df['Phase'], band_df['AbsMag'], c='k', marker='*', s=180, label=None)
        ax_uvm2.errorbar(band_df['Phase'], band_df['AbsMag'], yerr=band_df['AbsErr'], fmt='',
                         marker='*', c='k', ls='-', lw=1.2, capsize=3, capthick=1, label='_nolegend_')
    elif band == 'uvw2':
        ax_uvw2.scatter(band_df['Phase'], band_df['AbsMag'], c='k', marker='*', s=180, label=None)
        ax_uvw2.errorbar(band_df['Phase'], band_df['AbsMag'], yerr=band_df['AbsErr'], fmt='',
                         marker='*', c='k', ls='-', lw=1.2, capsize=3, capthick=1, label='_nolegend_')
    else:
        continue

ax_uvw1.legend(fontsize=14, markerscale=1.5, loc=1, frameon=False)
ax_uvw2.set_ylabel(r'$\rm M_{uvw2}$ [mag]', fontsize=16)

fig_uvabs.subplots_adjust(hspace=0.01)
fig_uvabs.savefig('PLOT_SWIFTUVAbsLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_uvabs)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The V-Band Absolute Magnitude Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket'))

fig_vabs = plt.figure(figsize=(9, 9))
ax_vabs = fig_vabs.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['1999em', '2004et', '1987A', '2013ab', '2009ib']:
        data_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        data_df = multicol_to_fluxdf(name, data_df)
        temp_df = data_df[data_df['FILTER'] == 'V'].copy()
        temp_df = temp_df[temp_df['Phase'] < max_epoch].sort_values(by='Phase')
        if not temp_df.empty:
            ax_vabs.plot(temp_df['Phase'], temp_df['AbsMag'], ls=':', markersize=7,
                         marker=sndata_df.loc[name, 'Marker'], label=name)

ax_vabs.plot(vabs_df['Phase'], vabs_df['AbsMag'], color='k', ms=8, markerfacecolor='dimgrey', markeredgewidth=1,
             marker='o', ls='', alpha=0.8, label=name_SN)
ax_vabs.plot(vabs_df['Phase'], vabs_df['AbsMag'], color='k', ms=3, markerfacecolor='None', markeredgewidth=1,
             marker='o', ls='-', lw=1, alpha=0.8, label='_nolegend_')

nebarr = np.arange(phase_nebstart + 50, phase_nebstart + 250, 1)
ax_vabs.plot(nebarr, 0.0098 * nebarr - 16.5, ls='--', lw=1.2, c='indianred')

ax_vabs.axvline(0, ls='--', c='dimgrey', lw=1.2)
ax_vabs.axvline(90, ls='--', c='dimgrey', lw=1.2)
ax_vabs.axvspan(0, 90, color='orange', alpha=0.1)

ax_vabs.text(20, -12.6, s='OPTd~{0:.0f} d'.format(90), fontsize=12)
ax_vabs.annotate(s='', xy=(0, -12.5), xytext=(90, -12.5), arrowprops=dict(arrowstyle='<->'))
ax_vabs.text(phase_nebstart + 60, -14.9, s=r'$\rm ^{56}Co \rightarrow\ ^{56}Fe$', fontsize=14, rotation=-32)

set_plotparams(ax_vabs)
ax_vabs.invert_yaxis()
ax_vabs.set_ylim(-11.4, -17.6)
ax_vabs.set_xlim(-10, max_epoch + 20)
ax_vabs.legend(fontsize=14, markerscale=2, loc=1, frameon=False)
ax_vabs.yaxis.set_major_locator(MultipleLocator(1))
ax_vabs.yaxis.set_minor_locator(MultipleLocator(0.1))
ax_vabs.set_ylabel(r'Absolute V-Band Magnitude, $\rm M_V$ [mag]', fontsize=16)

fig_vabs.savefig('PLOT_VAbsLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_vabs)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig_bol = plt.figure(figsize=(9, 9))
ax_bol = fig_bol.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['1999em', '2004et', '2013ab', '2005cs', '2009N']:
        data_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        temp_df = calc_boldf(name, multicol_to_fluxdf(name, data_df))
        if not temp_df.empty:
            temp_df = temp_df[temp_df['Phase'] < max_epoch].sort_values(by='Phase')
            temp_df['Lum'] = temp_df['Lum'].apply(lambda x: np.log10(x))
            ax_bol.plot(temp_df['Phase'], temp_df['Lum'], ls=':', markersize=6, c=sndata_df.loc[name, 'Color'],
                        alpha=0.9, marker=sndata_df.loc[name, 'Marker'], label=name)

ax_bol.plot(x, y, ls='', color='k', marker='o', markerfacecolor='grey',
            markeredgewidth=1, markersize=8, alpha=0.8, label=name_SN)
ax_bol.plot(x, y, ls='-.', lw=0.8, color='k', marker='o', markerfacecolor='None',
            markeredgewidth=1, markersize=3, alpha=0.8, label='_nolegend_')
ax_bol.plot(xarr, fit, ls='-', lw=2, c='blue', label='_nolegend_')

set_plotparams(ax_bol)
ax_bol.set_ylim(39.6, 42.5)
ax_bol.yaxis.set_major_locator(MultipleLocator(1))
ax_bol.yaxis.set_minor_locator(MultipleLocator(0.1))

yticks = ax_bol.get_yticks(minor=True)
ax_bol.text(xguess[0] - 13, 40.1, 'Rising Phase', color='dimgrey', alpha=0.9, fontsize=11, rotation=90)
ax_bol.text(xguess[1] - 21, 41.04, 'Early Plateau Phase [$s_1$={:.2f} dex / 100 d]'.format(s1), alpha=0.8, fontsize=11,
            rotation=90)
ax_bol.text(xguess[2] - 33, 41.02, 'Late Plateau Phase [$s_2$={:.2f} dex / 100 d]'.format(s2), alpha=0.8, fontsize=11,
            rotation=90)
ax_bol.text(xguess[3] - 20, 40.22, 'Transition Phase', color='dimgrey', alpha=0.9, fontsize=11, rotation=90)
ax_bol.text(xguess[3] + 42, 41.17, 'Nebular Phase [$s_3$={:.2f} dex / 100 d]'.format(s3), alpha=0.8, fontsize=11,
            rotation=-32)

# ax_bol.axvspan(-10, xguess[0], color='blue', alpha=0.2)
ax_bol.axvspan(xguess[0], xguess[1], color='grey', alpha=0.2)
ax_bol.axvspan(xguess[1], xguess[2], color='grey', alpha=0.2)
# ax_bol.axvspan(xguess[2], xguess[3], color='blue', alpha=0.3)
ax_bol.axvspan(xguess[3], max_epoch + 20, color='grey', alpha=0.2)

for xbreak in xguess:
    ax_bol.axvline(xbreak, ls='--', lw=0.8, c='k')
    ax_bol.axvline(xbreak, ls='--', lw=0.8, c='k')

ax_bol.set_xlim(-12, max_epoch + 20)
ax_bol.legend(fontsize=13, markerscale=2, loc=1, frameon=False)
ax_bol.set_ylabel(r'Log [$\rm L_{UBVRI}\ (erg\ s^{-1})$]', fontsize=16)

fig_bol.savefig('PLOT_BolometricLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_bol)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette(["#95a5a6", "#9b59b6", "#3498db", "#e74c3c", "#34495e", "#2ecc71"])[1:])

# sns.set_palette(sns.xkcd_palette(["windows blue", "amber", "greyish", "faded green", "dusty purple"]))
fig_bol = plt.figure(figsize=(16, 12))
ax_bol = fig_bol.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['1999em', '2004et', '2013ab', '2005cs', '2009N']:
        data_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        temp_df = calc_boldf(name, multicol_to_fluxdf(name, data_df))
        if not temp_df.empty:
            temp_df = temp_df[temp_df['Phase'] < max_epoch].set_index('Phase').sort_index()
            temp_df['Lum'] = median_filter(temp_df['Lum'].apply(lambda x: np.log10(x)), size=3)
            ax_bol.plot(temp_df.index.values, temp_df['Lum'], ls=':', ms=9, alpha=0.7,
                        marker=sndata_df.loc[name, 'Marker'], label=name)
            ax_bol.plot(temp_df.index.values, temp_df['Lum'], ls=':', c='dimgrey', markerfacecolor='None', ms=9,
                        alpha=0.7, markeredgewidth=0.7, marker=sndata_df.loc[name, 'Marker'], label='_nolegend_')

ax_bol.plot(x, y, ls='', color='k', marker='o', markerfacecolor='grey',
            markeredgewidth=1, markersize=12, alpha=0.9, label=name_SN)
ax_bol.plot(x, y, ls='-.', lw=0.8, color='k', marker='o', markerfacecolor='None',
            markeredgewidth=1, markersize=6, alpha=0.9, label='_nolegend_')
ax_bol.plot(xarr, fit, ls='-', lw=2.5, c='blue', label='_nolegend_')

ax_bol.set_ylim(39.6, 42.6)
ax_bol.set_xlim(-12, max_epoch + 100)
ax_bol.legend(fontsize=15, markerscale=2, loc='upper center', bbox_to_anchor=(0.5, 1.06),
              fancybox=True, shadow=True, ncol=6)
ax_bol.yaxis.set_ticks_position('both')
ax_bol.xaxis.set_ticks_position('both')
ax_bol.yaxis.set_major_locator(MultipleLocator(1))
ax_bol.yaxis.set_minor_locator(MultipleLocator(0.1))
ax_bol.xaxis.set_major_locator(MultipleLocator(100))
ax_bol.xaxis.set_minor_locator(MultipleLocator(10))
ax_bol.tick_params(axis='both', which='major', direction='in', length=10, width=1.8, labelsize=18)
ax_bol.tick_params(axis='both', which='minor', direction='in', length=5, width=0.9, labelsize=18)
ax_bol.set_xlabel('Time Since Explosion [Days]', fontsize=18)
ax_bol.set_ylabel(r'Log [$\rm L_{UBVRI}\ (erg\ s^{-1})$]', fontsize=18)

yticks = ax_bol.get_yticks(minor=True)
ax_bol.text(xguess[0] - 12, 40.08, 'Rising Phase', color='navy', alpha=0.9, fontsize=14, rotation=90)
ax_bol.text(xguess[1] - 21, 41.03, 'Early Plateau Phase [$s_1$={:.2f} dex / 100 d]'.format(s1), color='brown',
            alpha=0.8, fontsize=14, rotation=90)
ax_bol.text(xguess[2] - 33, 41.00, 'Late Plateau Phase [$s_2$={:.2f} dex / 100 d]'.format(s2), color='brown',
            alpha=0.8, fontsize=14, rotation=90)
ax_bol.text(xguess[3] - 20, 40.21, 'Transition Phase', color='navy', alpha=0.9, fontsize=14, rotation=90)
ax_bol.text(xguess[3] + 132, 40.77, 'Nebular Phase [$s_3$={:.2f} dex / 100 d]'.format(s3), color='brown',
            alpha=0.8, fontsize=14, rotation=-30)

# ax_bol.axvspan(-10, xguess[0], color='blue', alpha=0.2)
# ax_bol.axvspan(xguess[2], xguess[3], color='blue', alpha=0.3)
ax_bol.axvspan(xguess[0], xguess[1], color='grey', alpha=0.2)
ax_bol.axvspan(xguess[1], xguess[2], color='grey', alpha=0.2)
ax_bol.axvspan(xguess[3], max_epoch + 100, color='grey', alpha=0.2)

for xbreak in xguess:
    ax_bol.axvline(xbreak, ls='--', lw=0.8, c='k')
    ax_bol.axvline(xbreak, ls='--', lw=0.8, c='k')

axins = inset_axes(ax_bol, height=4.5, width=7.7, loc=1)

axins.plot(obsbol_df['Phase'], obsbol_df['LogLum'], c='k', markerfacecolor='dimgrey',
           marker='*', ls='', ms=15, label='Observed LC')
axins.plot(comb_df['Phase'], comb_df['LogLumCore'], c='r', ls='-', lw=0.9, label='Model Core LC (a=0)')
axins.plot(comb_df['Phase'], comb_df['LogLumShell'], c='g', ls=':', lw=0.9, label='Model Shell LC (n=2)')
axins.plot(comb_df['Phase'], comb_df['LogLum'], c='b', ls='--', lw=1.5, label='Model Combined LC')

axins.set_ylim(40.75, 42.35)
axins.set_xlim(-5, 210)
axins.legend(markerscale=1.5, fontsize=13, frameon=False)
axins.yaxis.set_ticks_position('both')
axins.xaxis.set_ticks_position('both')
axins.xaxis.set_major_locator(MultipleLocator(50))
axins.xaxis.set_minor_locator(MultipleLocator(5))
axins.yaxis.set_major_locator(MultipleLocator(0.5))
axins.yaxis.set_minor_locator(MultipleLocator(0.05))
axins.tick_params(axis='both', which='major', direction='in', width=1.4, length=8, labelsize=16)
axins.tick_params(axis='both', which='minor', direction='in', width=0.8, length=4, labelsize=16)

fig_bol.savefig('PLOT_BolometricNagyLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_bol)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Temperature And Radius Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig_temp = plt.figure(figsize=(9, 8))
ax_temp = fig_temp.add_subplot(111)
ax_rad = ax_temp.twinx()

ax_temp.plot(temprad_df['Phase'], temprad_df['Temp'], ls='', lw=0.8, c='orange', marker='D', ms=7, label='Temperature')
ax_temp.errorbar(temprad_df['Phase'], temprad_df['Temp'], yerr=temprad_df['TempErr'], c='grey', marker='o', capsize=3,
                 capthick=1, ms=2, alpha=0.8, elinewidth=1, label='_nolegend_')

ax_rad.plot(temprad_df['Phase'], temprad_df['Rad'], ls='', lw=0.8, c='dodgerblue', marker='o', ms=7,
            label='Radius [BB Fit]')
ax_rad.errorbar(temprad_df['Phase'], temprad_df['Rad'], yerr=temprad_df['RadErr'], c='grey', marker='o', capsize=3,
                capthick=1, ms=2, alpha=0.8, elinewidth=1, label='_nolegend_')

ax_rad.errorbar(vel_df['Phase'], vel_df['Rad'], yerr=vel_df['RadErr'], c='r', marker='s', capsize=3, lw=1, ls=':',
                capthick=1, ms=6, alpha=0.8, elinewidth=1, label=r'Radius [$\rm v_{Fe\,II}\ t$]')

ax_temp.axvline(8, ls='--', lw=0.8, c='k')
ax_temp.axvline(25, ls='--', lw=0.8, c='k')
ax_temp.axvspan(8, 25, color='grey', alpha=0.2)

ax_temp.xaxis.set_ticks_position('both')
ax_temp.xaxis.set_major_locator(MultipleLocator(25))
ax_temp.xaxis.set_minor_locator(MultipleLocator(5))
ax_temp.yaxis.set_major_locator(MultipleLocator(2))
ax_temp.yaxis.set_minor_locator(MultipleLocator(0.4))
ax_rad.yaxis.set_major_locator(MultipleLocator(4))
ax_rad.yaxis.set_minor_locator(MultipleLocator(0.8))
ax_temp.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_temp.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)
ax_rad.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_rad.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)

ax_temp.set_xlim(-5, 150)
ax_temp.set_ylim(3.5, 19)
ax_temp.legend(fontsize=12, markerscale=1.5, loc=2)
ax_rad.legend(fontsize=12, markerscale=1.5, loc=1)

ax_temp.text(15, 17, s='Ejecta-CSM Interaction', color='k', alpha=0.7, rotation=90, fontsize=12)
ax_temp.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_temp.set_ylabel(r'Temperature [$\rm \times 10^3\ K$]', fontsize=16)
ax_rad.set_ylabel(r'Radius [$\rm \times 10^3\ R_{\odot}$]', fontsize=16)

fig_temp.savefig('PLOT_TempRadEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_temp)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Temperature And Radius Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig_temp = plt.figure(figsize=(10, 9))
ax_temp = fig_temp.add_subplot(111)
ax_rad = ax_temp.twinx()

phase = temprad_df['Phase'].astype('float64').values
temp = temprad_df['Temp'].astype('float64').values
temperr = temprad_df['TempErr'].astype('float64').values
rad = temprad_df['Rad'].astype('float64').values
raderr = temprad_df['RadErr'].astype('float64').values

ax_temp.plot(phase, temp, ls='-', lw=0.8, c='forestgreen', markerfacecolor='dimgrey', markeredgewidth=1.5, marker='D',
             ms=9, label='Temperature')
ax_temp.fill_between(phase, temp - temperr, temp + temperr, color='dimgrey', alpha=0.4)

ax_rad.plot(phase, rad, ls='--', lw=0.8, c='blue', markerfacecolor='dimgrey', marker='o', ms=9, markeredgewidth=2,
            label='Radius [BB Fit]')
ax_rad.fill_between(phase, rad - raderr, rad + raderr, color='grey', alpha=0.4)

ax_rad.plot(vel_df['Phase'], vel_df['Rad'], c='k', markerfacecolor='dimgrey', markeredgewidth=2, marker='*', lw=1,
            ls=':', ms=12, label=r'Radius [$\rm v_{Fe\,II}\ t$]')
ax_rad.fill_between(vel_df['Phase'], vel_df['Rad'] - vel_df['RadErr'], vel_df['Rad'] + vel_df['RadErr'],
                    color='tomato', alpha=0.4)

ax_temp.axvline(8, ls='--', lw=0.8, c='k')
ax_temp.axvline(25, ls='--', lw=0.8, c='k')
ax_temp.axvspan(8, 25, color='grey', alpha=0.2)

ax_temp.xaxis.set_ticks_position('both')
ax_temp.xaxis.set_major_locator(MultipleLocator(25))
ax_temp.xaxis.set_minor_locator(MultipleLocator(5))
ax_temp.yaxis.set_major_locator(MultipleLocator(2))
ax_temp.yaxis.set_minor_locator(MultipleLocator(0.4))
ax_rad.yaxis.set_major_locator(MultipleLocator(4))
ax_rad.yaxis.set_minor_locator(MultipleLocator(0.8))
ax_temp.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_temp.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)
ax_rad.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_rad.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)

ax_temp.set_xlim(-5, 115)
ax_temp.set_ylim(3.5, 21)
ax_temp.legend(fontsize=12, markerscale=1.5, loc=2)
ax_rad.legend(fontsize=12, markerscale=1.5, loc='upper center')

ax_temp.text(15, 18.5, s='Ejecta-CSM Interaction', color='k', alpha=0.7, rotation=90, fontsize=12)
ax_temp.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_temp.set_ylabel(r'Temperature [$\rm \times 10^3\ K$]', fontsize=16)
ax_rad.set_ylabel(r'Radius [$\rm \times 10^3\ R_{\odot}$]', fontsize=16)

fig_temp.savefig('PLOT_TempRadEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_temp)
# ------------------------------------------------------------------------------------------------------------------- #
