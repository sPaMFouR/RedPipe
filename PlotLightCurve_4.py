#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx-------------------PLOT THE SUPERNOVA LIGHT CURVES----------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import math
import numpy as np
import pandas as pd
from datetime import date
import matplotlib.pyplot as plt
from jdcal import jd2gcal, gcal2jd
from scipy.interpolate import CubicSpline
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
fmt_flt = '{0:>7.3f}'
fmt_exp = '{0:>7.4e}'
epoch = 2400000.5
wave_data = np.linspace(3100, 9200, 1000)
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
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/DataExt/IIP_Data/"
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

data = pd.read_csv(DIR_SNe + 'LC_Data/TypeIISNe.dat', sep='\s+', comment='#')
data = data.replace('INDEF', np.nan).set_index(['Name', 'Marker', 'Color']).astype('float64')
data = data.reset_index().set_index('Name')
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
# Functions For Manipulating Pandas DataFrame Containing Data From Well Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def coltorow_df(name, data_df):
    """
    Converts a column-wise magnitude Pandas DataFrame to a row-wise Pandas DataFrame.
    Args:
        name        : Name of the SNe whose data is read
        data_df     : Input Pandas DataFrame
    Returns:
        output_df   : Output Pandas DataFrame
    """

    def calc_magflux(name, mag, err, band):
        mag = float(mag)
        err = float(err)
        zp = filter_df.loc[band, 'ZeroPoint']
        rlambda = filter_df.loc[band, 'RLambda']

        distmod_mag = 5 * np.log10(data.loc[name, 'D'] * 10 ** 6) - 5
        distmod_err = 5 * np.log10((data.loc[name, 'D'] + data.loc[name, 'DErr']) * 10 ** 6) - 5 - distmod_mag

        absmag = fmt_flt.format(mag - rlambda * data.loc[name, 'EBV'] - distmod_mag)
        abserr = fmt_flt.format((err ** 2 + (rlambda * data.loc[name, 'EBVErr']) ** 2 + distmod_err ** 2) ** 0.5)

        flux = float(fmt_exp.format(10 ** (-0.4 * (mag - rlambda * data.loc[name, 'EBV'] + zp + 21.10))))
        fluxerr = fmt_exp.format(abs(flux - 10 ** (-0.4 * (mag + err - rlambda * data.loc[name, 'EBV'] + zp + 21.10))))

        return float(absmag), float(abserr), float(flux), float(fluxerr)

    data_df = data_df.set_index('JD')
    data_df = data_df.drop(['Date', 'Phase'], axis=1)
    data_arr = data_df.as_matrix()

    size = data_arr.shape
    list_jd = np.repeat(data_df.index.values, (size[1] / 2))
    list_filters = [x for x in data_df.columns.values if 'Err' not in x]
    data_arr = np.reshape(data_arr, [size[0] * size[1] / 2, 2])

    input_df = pd.DataFrame(data_arr, index=list_jd, columns=['FMAG', 'FERR'])
    input_df.index.name = 'JD'
    input_df = input_df.reset_index(drop=False)

    input_df['FILTER'] = list_filters * size[0]
    input_df['Date'] = input_df['JD'].apply(jd_to_cald)
    input_df['Phase'] = input_df['JD'] - data.loc[name, 'DateExp']

    input_df = input_df.replace('INDEF', np.nan).dropna(axis=0, how='any')
    input_df = input_df[['Date', 'JD', 'Phase', 'FILTER', 'FMAG', 'FERR']].reset_index(drop=True)
    input_df['ALambda'] = input_df['FILTER'].apply(lambda x: float(fmt_flt.format(
        filter_df.loc[x, 'RLambda'] * data.loc[name, 'EBV'])))

    for index, band in input_df['FILTER'].items():
        magflux = calc_magflux(name, mag=input_df.loc[index, 'FMAG'], err=input_df.loc[index, 'FERR'], band=band)
        input_df.loc[index, 'AbsMag'] = magflux[0]
        input_df.loc[index, 'AbsErr'] = magflux[1]
        input_df.loc[index, 'Flux'] = magflux[2]
        input_df.loc[index, 'FluxErr'] = magflux[3]

    return input_df


def obtain_epochwisedf(input_df, flux='Flux', fluxerr='FluxErr'):
    """
    Creates a Pandas DataFrame with magnitudes arranged epoch-wise from a DataFrame with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing fluxes
        flux        : Name of the Flux column in the Pandas DataFrame
        fluxerr     : Name of the Flux Error column in the Pandas DataFrame
    Returns:
        bolm_df     : Pandas DataFrame containing organised magnitudes
        flux_df     : Pandas DataFrame containing organised broadband fluxes
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
            dict_val[day][band] = float(np.mean(list_flux))

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 2:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if len(band) == 1:
                    dict_flux[day][filter_df.loc[band, 'CentreWave']] = flux

    bolm_df = pd.DataFrame(dict_flux).T
    bolm_df.index.name = 'JD'
    bolm_df = bolm_df.interpolate(method='linear', limit=3)
    bolm_df = bolm_df.T

    return bolm_df


def calc_boldf(name, input_df, plot=False):
    """
    Creates a Pandas DataFrame with Bolometric Fluxes.
    Args:
        name        : Name of the SNe whose bolometric magnitudes are to be computed
        input_df    : Input Pandas DataFrame containing individual band fluxes
        plot        : Whether the spline fits to the fluxes should be plotted
    Returns:
        output_df   : Output Pandas DataFrame containing bolometric fluxes
    """

    def calc_luminosity(name, flux):
        val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
        lum = fmt_exp.format(val *
                             [name, 'D'] ** 2)
        lumerr = fmt_exp.format(val * ((data.loc[name, 'D'] + data.loc[name, 'DErr']) ** 2 - data.loc[name, 'D'] ** 2))
        return float(lum), float(lumerr)

    dict_flux = {}

    for jd in input_df.columns.values:
        series = input_df[jd].dropna().apply(lambda x: float(x))
        if len(series) <= 2:
            continue

        spline = CubicSpline(series.index.values.tolist(), series.values.tolist(), bc_type='natural', extrapolate=True)
        #         spline = Rbf(series.index.values.tolist(), series.values.tolist())

        if jd - data.loc[name, 'DateExp'] > 99:
            wave_data = np.linspace(4000, 9200, 1000)
        else:
            wave_data = np.linspace(3100, 9200, 1000)

        flux_data = spline(wave_data)
        flux_data[flux_data < 0] = 0
        netflux = np.trapz(flux_data, wave_data)

        dict_flux[jd] = {}
        dict_flux[jd]['Date'] = jd_to_cald(jd)
        dict_flux[jd]['Phase'] = jd - data.loc[name, 'DateExp']
        dict_flux[jd]['Flux'] = netflux
        dict_flux[jd]['Lum'] = calc_luminosity(name, netflux)[0]
        dict_flux[jd]['LumErr'] = calc_luminosity(name, netflux)[1]

        if plot:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wave_data, spline(wave_data), 'r-', label='CubicSpline Fit')
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig)

    if len(dict_flux) != 0:
        output_df = pd.DataFrame(dict_flux).T
        output_df.index.name = 'JD'
        output_df = output_df.reset_index().set_index(keys='Date', drop=True)
    else:
        output_df = pd.DataFrame()

    return output_df


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Manipulating Pandas DataFrames Containing Data From SN In Study
# ------------------------------------------------------------------------------------------------------------------- #

def organise_dataframe(input_df, flux='Flux', fluxerr='FluxErr', tempsub=False, plot=False):
    """
    Creates a Pandas DataFrame with magnitudes arranged epoch-wise from a DataFrame with unorganised magnitudes.
    Args:
        input_df    : Pandas DataFrame containing magnitudes
        flux        : Name of Flux column in the Pandas DataFrame
        fluxerr     : Name of Flux Error column in the Pandas DataFrame
        tempsub     : Whether the magnitudes are template subtracted?
        plot        : Whether the spline fits to the fluxes should be plotted
    Returns:
        bolm_df     : Pandas DataFrame containing organised
        flux_df     : Pandas DataFrame containing organised broadband fluxes
    """

    def calc_lum(flux):
        val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
        lum = fmt_exp.format(val * dist_val ** 2)
        lumerr = fmt_exp.format(val * ((dist_val + dist_err) ** 2 - dist_val ** 2))
        return float(lum), float(lumerr)

    def get_nickelmass(lum, phase):
        hamuy_ni = 7.866e-44 * lum * np.exp((phase * (1 + redshift) - 6.1) / 111.26)
        jerk_ni = (lum * 0.07 / 9.92e41) / (math.exp(-phase / 111.4) - math.exp(-phase / 8.8))
        return hamuy_ni, jerk_ni

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
            dict_val[day][band] = float(np.mean(list_flux))

    flux_df = pd.DataFrame(dict_val).T
    flux_df.index.name = 'JD'
    flux_df = flux_df.reset_index()
    flux_df['Date'] = flux_df['JD'].apply(jd_to_cald)
    flux_df['Phase'] = flux_df['JD'].apply(lambda x: x - date_explosion).round(int(precision))
    flux_df.to_csv('OUTPUT_DateWiseSNAppFlux' + add_text, sep=' ', index=True, na_rep='INDEF')

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 2:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if len(band) == 1:
                    dict_flux[day][filter_df.loc[band, 'CentreWave']] = flux

    bolm_df = pd.DataFrame(dict_flux).T
    bolm_df.index.name = 'JD'
    bolm_df = bolm_df.interpolate(method='linear', limit=2)
    bolm_df = bolm_df.T

    dict_flux = {}
    for caljd in bolm_df.columns.values:
        series = bolm_df[caljd].dropna().apply(lambda x: float(x))
        if len(series) <= 2:
            continue
        spline = CubicSpline(series.index.values.tolist(), series.values.tolist(), bc_type='natural', extrapolate=True)
        #         spline3 = Rbf(series.index.values.tolist(), series.values.tolist())

        if caljd - date_explosion > 99:
            wave_data = np.linspace(4000, 9200, 1000)
        else:
            wave_data = np.linspace(3100, 9200, 1000)

        flux_data = spline(wave_data)
        flux_data[flux_data < 0] = 0
        netflux = np.trapz(flux_data, wave_data)

        dict_flux[caljd] = {}
        dict_flux[caljd]['Date'] = jd_to_cald(caljd)
        dict_flux[caljd]['Phase'] = caljd - date_explosion
        dict_flux[caljd]['Flux'] = netflux
        dict_flux[caljd]['Lum'] = calc_lum(netflux)[0]
        dict_flux[caljd]['LumErr'] = calc_lum(netflux)[1]

        if plot:
            fig_temp = plt.figure(figsize=(8, 6))
            ax = fig_temp.add_subplot(111)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wave_data, spline(wave_data), 'r-', label='CubicSpline Fit')
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig_temp)

    fbolm_df = pd.DataFrame(dict_flux).T
    fbolm_df.index.name = 'Date'
    fbolm_df.to_csv('OUTPUT_DateWiseSNBolFlux' + add_text, sep=' ', index=True)

    for index, row in fbolm_df.iterrows():
        if row['Phase'] > 160:
            fbolm_df.loc[index, 'MNi'] = get_nickelmass(row['Lum'], phase=row['Phase'])[0]
            fbolm_df.loc[index, 'MNi2'] = get_nickelmass(row['Lum'], phase=row['Phase'])[1]
            fbolm_df.loc[index, 'MNiErr'] = get_nickelmass(row['LumErr'], phase=row['Phase'])[0]

    return bolm_df, fbolm_df


def calc_fluxdf(file_name, concat_df=pd.DataFrame(), tempsub=False):
    """
    Creates a Pandas DataFrame with Individual Band Fluxes.
    Args:
        file_name   : Name of the file which has the photometric magnitudes for the SNe
        concat_df   : Input Pandas DataFrame to be appended to the original data
        tempsub     : Boolean describing whether these magnitudes are template subtracted?
    Returns:
        output_df   : Output Pandas DataFrame containing individual band fluxes
    """

    def calc_absmagflux(mag, err, band):
        zp = filter_df.loc[band, 'ZeroPoint']
        rlambda = filter_df.loc[band, 'RLambda']
        absmag = fmt_flt.format(mag - rlambda * EBV_mag - distmod_mag)
        abserr = fmt_flt.format((err ** 2 + (rlambda * EBV_err) ** 2 + distmod_err ** 2) ** 0.5)
        flux = fmt_exp.format(10 ** (-0.4 * (mag - rlambda * EBV_mag + zp + 21.100)))
        fluxerr = fmt_exp.format(abs(float(flux) - 10 ** (-0.4 * (mag + err - rlambda * EBV_mag + zp + 21.100))))

        return float(absmag), float(abserr), float(flux), float(fluxerr)

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
        data_magflux = calc_absmagflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'], band=band)
        output_df.loc[index, 'AbsMag'] = data_magflux[0]
        output_df.loc[index, 'AbsErr'] = data_magflux[1]
        output_df.loc[index, 'Flux'] = data_magflux[2]
        output_df.loc[index, 'FluxErr'] = data_magflux[3]

    output_df['ALambda'] = output_df['FILTER'].apply(
        lambda x: float(fmt_flt.format(filter_df.loc[x, 'RLambda'] * EBV_mag)))
    output_df.to_csv('OUTPUT_NetSNFlux' + add_text, sep=' ', index=False)

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Other SNe Data From Archive Folder
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', DIR_SNe + 'LC_Data/*.asc')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
rawopt_df = calc_fluxdf('OUTPUT_FinalSNMag', tempsub=False)
outputopt_df = calc_fluxdf('OUTPUT_FinalSNMagTemp', tempsub=True)

max_epoch = outputopt_df['Phase'].max()
outputopt_df = outputopt_df[~((outputopt_df['Date'] == '2016-12-29') & (outputopt_df['FILTER'] == 'U'))]
outputopt_df = outputopt_df[~((outputopt_df['Date'] == '2017-10-01') & (outputopt_df['FILTER'] == 'U'))]
outputopt_df[['JD', 'Phase']] = outputopt_df[['JD', 'Phase']].round(1)

vabs_df = outputopt_df[outputopt_df['FILTER'] == 'V'].copy()
nebular_df = outputopt_df[(outputopt_df['Phase'] > 110) & (outputopt_df['Phase'] < 400)].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read SWIFT Data And Combine It With Optical Data To Calculate Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
outputuv_df = pd.read_csv(DIR_PHOT + '2016gfy_Swift.dat', sep='\s+', comment='#')
outputuv_df = outputuv_df[['FILTER', 'JD', 'FMAG', 'FERR']].replace('INDEF', np.nan)
outputuv_df = outputuv_df.set_index('FILTER').astype('float64').reset_index()
outputuv_df['JD'] = (outputuv_df['JD'] + 2400000).round(1)
outputuv_df['Date'] = outputuv_df['JD'].apply(jd_to_cald)
outputuv_df['Phase'] = (outputuv_df['JD'] - date_explosion).round(1)
outputuv_df = outputuv_df.set_index('JD')

outputnet_df = pd.concat([outputopt_df.set_index('JD'), outputuv_df])

dict_mag = {}
for index, row in outputnet_df.iterrows():
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

outputnet_df = pd.DataFrame(dict_mag).T
outputnet_df.index.name = 'JD'
outputnet_df = outputnet_df.reset_index()
outputnet_df['Phase'] = (outputnet_df['JD'] - date_explosion).round(1)
interpnet_df = outputnet_df.interpolate(method='linear', limit=2)

# bolm_df, fbolm_df = organise_dataframe(interpnet_df)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Optical Apparent Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_app = plt.figure(figsize=(8, 8))
ax_app = fig_app.add_subplot(111)

for band, band_df in outputopt_df.groupby('FILTER'):
    ax_app.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                   marker=filter_df.loc[band, 'Marker'], s=20, label=filter_df.loc[band, 'Label'], alpha=0.5)
    ax_app.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                    c=filter_df.loc[band, 'Color'], linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

handles, labels = ax_app.get_legend_handles_labels()
handles = [handles[3], handles[0], handles[4], handles[2], handles[1]]
labels = [labels[3], labels[0], labels[4], labels[2], labels[1]]
ax_app.legend(handles, labels, fontsize=12, markerscale=2, loc=1, frameon=False)

ax_app.set_ylim(22.6, 14.5)
ax_app.set_xlim(-15, max_epoch + 20)
ax_app.yaxis.set_ticks_position('both')
ax_app.xaxis.set_ticks_position('both')
ax_app.yaxis.set_major_locator(MultipleLocator(2))
ax_app.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_app.xaxis.set_major_locator(MultipleLocator(100))
ax_app.xaxis.set_minor_locator(MultipleLocator(10))
ax_app.set_ylabel('Apparent Magnitude [mag]', fontsize=16)
ax_app.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_app.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_app.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

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
                     marker=filter_df.loc[band, 'Marker'], s=30, label=filter_df.loc[band, 'Label'])
    ax_swift.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                      c=filter_df.loc[band, 'Color'], linestyle='--', linewidth=1, capsize=2, capthick=1, label=None)

handles, labels = ax_swift.get_legend_handles_labels()
handlesuv = [handles[5], handles[1], handles[4], handles[2], handles[0], handles[3]]
labelsuv = [labels[5], labels[1], labels[4], labels[2], labels[0], labels[3]]
ax_swift.legend(handlesuv, labelsuv, fontsize=11, markerscale=2, loc=4, frameon=False)

ax_swift.set_ylim(22.6, 14.5)
ax_swift.set_xlim(-1, 25)
ax_swift.yaxis.set_ticks_position('both')
ax_swift.xaxis.set_ticks_position('both')
ax_swift.yaxis.set_major_locator(MultipleLocator(2))
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
fig_comb, (ax_opt, ax_swift) = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

for band, band_df in outputopt_df.groupby('FILTER'):
    ax_opt.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                   marker=filter_df.loc[band, 'Marker'], s=30, label=filter_df.loc[band, 'Label'], alpha=0.5)
    ax_opt.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                    c=filter_df.loc[band, 'Color'], linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

handles, labels = ax_opt.get_legend_handles_labels()
handles = [handles[3], handles[0], handles[4], handles[2], handles[1]]
labels = [labels[3], labels[0], labels[4], labels[2], labels[1]]
ax_opt.legend(handles, labels, fontsize=12, markerscale=2, loc=1)

ax_opt.set_title('HCT-HFOSC', fontsize=16)
ax_opt.set_ylim(22.6, 14.5)
ax_opt.set_xlim(-15, max_epoch + 20)
ax_opt.yaxis.set_ticks_position('both')
ax_opt.xaxis.set_ticks_position('both')
ax_opt.yaxis.set_major_locator(MultipleLocator(2))
ax_opt.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_opt.xaxis.set_major_locator(MultipleLocator(100))
ax_opt.xaxis.set_minor_locator(MultipleLocator(10))
ax_opt.set_ylabel('Apparent Magnitude [mag]', fontsize=16)
ax_opt.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_opt.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_opt.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

for band, band_df in outputuv_df.groupby('FILTER'):
    ax_swift.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], c=filter_df.loc[band, 'Color'],
                     marker=filter_df.loc[band, 'Marker'], s=30, label=filter_df.loc[band, 'Label'])
    ax_swift.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'], fmt='',
                      c=filter_df.loc[band, 'Color'], linestyle='', linewidth=1, capsize=2, capthick=1, label=None)

handles, labels = ax_swift.get_legend_handles_labels()
handlesuv = [handles[5], handles[1], handles[4], handles[2], handles[0], handles[3]]
labelsuv = [labels[5], labels[1], labels[4], labels[2], labels[0], labels[3]]
ax_swift.legend(handlesuv, labelsuv, fontsize=12, markerscale=2, loc=4)

ax_swift.set_title('SWIFT-UVOT', fontsize=16)
ax_swift.set_xlim(-1, 25)
ax_swift.set_yticklabels([])
ax_swift.yaxis.set_ticks_position('both')
ax_swift.xaxis.set_ticks_position('both')
ax_swift.yaxis.set_major_locator(MultipleLocator(2))
ax_swift.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_swift.xaxis.set_major_locator(MultipleLocator(5))
ax_swift.xaxis.set_minor_locator(MultipleLocator(0.5))
ax_swift.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_swift.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_swift.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

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
                     c=filter_df.loc[band, 'Color'], fmt='', linestyle='', linewidth=0.5, capsize=2, capthick=1,
                     label=None)

for band, band_df in rawopt_df.groupby('FILTER'):
    ax_comp.scatter(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], marker='+',
                    c='lightsalmon', s=10, label=filter_df.loc[band, 'Label'], alpha=0.5)
    ax_comp.errorbar(band_df['Phase'], band_df['FMAG'] + filter_df.loc[band, 'Offset'], yerr=band_df['FERR'],
                     c='lightsalmon', fmt='', linestyle='--', linewidth=0.5, capsize=2, capthick=1, label=None)

ax_comp.legend(markerscale=2, frameon=False)
ax_comp.set_ylim(22.5, 14.5)
ax_comp.set_xlim(-10, max_epoch + 20)
ax_comp.yaxis.set_ticks_position('both')
ax_comp.xaxis.set_ticks_position('both')
ax_comp.yaxis.set_major_locator(MultipleLocator(2))
ax_comp.yaxis.set_minor_locator(MultipleLocator(0.2))
ax_comp.xaxis.set_major_locator(MultipleLocator(50))
ax_comp.xaxis.set_minor_locator(MultipleLocator(5))
ax_comp.set_ylabel('Apparent Magnitude [mag]', fontsize=16)
ax_comp.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_comp.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_comp.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

fig_comp.savefig('PLOT_TempCompLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_comp)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The V-Band Absolute Magnitude Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig_vabs = plt.figure(figsize=(8, 8))
ax_vabs = fig_vabs.add_subplot(111)

for file_name in list_files:
    name = file_name.split('/')[-1].split('.')[0]
    if name in ['1999em', '2004et', '2009N', '2012ec', '2012aw', '2013ab']:
        data_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')
        data_df = coltorow_df(name, data_df)
        temp_df = data_df[data_df['FILTER'] == 'V'].copy()
        temp_df = temp_df[temp_df['Phase'] < max_epoch].sort_values(by='Phase')

        ax_vabs.plot(temp_df['Phase'], temp_df['AbsMag'], linestyle=':', markersize=6, marker=data.loc[name, 'Marker'],
                     c=data.loc[name, 'Color'], label=name)

ax_vabs.plot(vabs_df['Phase'], vabs_df['AbsMag'], color='k', markersize=5, marker='o', linestyle='-', label=name_SN)

ax_vabs.invert_yaxis()
ax_vabs.set_ylim(-10.8, -17.7)
ax_vabs.set_xlim(-5, max_epoch + 20)
ax_vabs.legend(fontsize=12, markerscale=2, loc=1, frameon=False)

ax_vabs.yaxis.set_ticks_position('both')
ax_vabs.xaxis.set_ticks_position('both')
ax_vabs.yaxis.set_major_locator(MultipleLocator(1))
ax_vabs.yaxis.set_minor_locator(MultipleLocator(0.1))
ax_vabs.xaxis.set_major_locator(MultipleLocator(100))
ax_vabs.xaxis.set_minor_locator(MultipleLocator(10))
ax_vabs.set_ylabel(r'Absolute V-Band Magnitude [$\rm M_V$]', fontsize=16)
ax_vabs.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax_vabs.tick_params(axis='both', which='major', direction='in', length=8, width=1, labelsize=14)
ax_vabs.tick_params(axis='both', which='minor', direction='in', length=4, width=1, labelsize=14)

fig_vabs.savefig('PLOT_VAbsLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_vabs)
# ------------------------------------------------------------------------------------------------------------------- #
