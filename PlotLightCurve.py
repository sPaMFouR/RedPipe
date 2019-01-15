#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx--------------Plot The Supernova Light Curve From Input Magnitude Files---------------xxxxxxxxxxxxxxx #
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
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import CubicSpline, Rbf
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
fmt_flt = "{0:>7.3f}"
fmt_exp = "{0:>7.4e}"
epoch = 2400000.5
wave_data = np.linspace(3100, 9200, 1000)
list_filters = ['U', 'B', 'V', 'R', 'I', 'g', 'r', 'i']
dict_markers = {'U': [+2.0, 'o', 'k', 'U + 2.0'], 'B': [+0.5, 'D', 'b', 'B + 0.5'], 'V': [+0.5, 's', 'g', 'V + 0.5'],
                'R': [+0.0, '^', 'r', 'R      '], 'I': [-0.5, 'p', 'm', 'I - 0.5'], 'g': [+1.0, '*', 'y', 'g + 1.0'],
                'r': [-1.0, 'P', 'c', 'r - 1.0'], 'i': [-3.0, 'X', 'indigo', 'i - 3.0']}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/"
DIR_PHOT = "/home/avinash/Supernovae_Data/2016gfy/PhotTemp/"
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
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def group_similar_files(text_list, common_text, exceptions=''):
    """
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
                test = re.search(str(text), file_name)
                if test:
                    try:
                        list_files.remove(file_name)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(str(text_list), "w") as f:
            for index in range(0, len(list_files)):
                f.write(str(list_files[index]) + "\n")

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
        zp = dict_zp[band]
        rlambda = dict_R[band]
        
        distmod_mag = 5 * np.log10(dict_dist[name][0] * 10 ** 6) - 5
        distmod_err = 5 * np.log10((dict_dist[name][0] +  dict_dist[name][1]) * 10 ** 6) - 5 - distmod_mag

        absmag = fmt_flt.format(mag - rlambda * dict_EBV[name][0] - distmod_mag)
        abserr = fmt_flt.format((err ** 2 + (rlambda * dict_EBV[name][1]) ** 2 + distmod_err ** 2) ** 0.5)
        
        flux = fmt_exp.format(10 ** (-0.4 * (mag - rlambda * dict_EBV[name][0] + zp + 21.100)))
        fluxerr = fmt_exp.format(abs(float(flux) - 10 ** (-0.4 * (mag + err - rlambda * dict_EBV[name][0] + 
                                                                  zp + 21.100))))

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
    input_df['Phase'] = input_df['JD'] - dict_explosion[name]
    input_df = input_df.replace('INDEF', np.nan).dropna(axis=0, how='any')
    input_df = input_df[['Date', 'JD', 'Phase', 'FILTER', 'FMAG', 'FERR']].reset_index(drop=True)
    input_df['ALambda'] = input_df['FILTER'].apply(lambda x: float(fmt_flt.format(dict_R[x] * dict_EBV[name][0])))

    for index, band in input_df['FILTER'].iteritems():
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
                    dict_flux[day][dict_centralwav[band]] = flux

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
    def calc_lumerr(name, flux):
        val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
        lum = fmt_exp.format(val * dict_dist[name][0] ** 2)
        lumerr = fmt_exp.format(val * ((dict_dist[name][0] + dict_dist[name][1]) ** 2 - dict_dist[name][0] ** 2))
        return float(lum), float(lumerr)

    dict_flux = {}

    for jd in input_df.columns.values:
        series = input_df[jd].dropna().apply(lambda x: float(x))
#         spline = CubicSpline(series.index.values.tolist(), series.values.tolist(), bc_type='natural', extrapolate=True)
        spline = Rbf(series.index.values.tolist(), series.values.tolist())

        netflux = np.trapz(spline(wave_data), wave_data)
        lumerr = calc_lumerr(name, netflux)
        
        dict_flux[jd] = {}
        dict_flux[jd]['Date'] = jd_to_cald(jd)
        dict_flux[jd]['Phase'] = jd - dict_explosion[name]
        dict_flux[jd]['Flux'] = netflux
        dict_flux[jd]['Lum'] = lumerr[0]
        dict_flux[jd]['Err'] = lumerr[1]

        if plot:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wave_data, spline(wave_data), 'r-', label='Rbf Fit')
#             ax.plot(wave_data, spline(wave_data), 'r-', label='CubicSpline Fit')
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig)

    output_df = pd.DataFrame(dict_flux).T
    output_df.index.name = 'JD'
    output_df = output_df.reset_index().set_index(keys='Date', drop=True)

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Manipulating Pandas DataFrames Containing Data From SNe In Study
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
        jerk_ni = (lum * 0.07 / 9.92e41) / (math.exp(-phase/111.4) - math.exp(-phase/8.8))
        return hamuy_ni, jerk_ni

    if not tempsub:
        add_text = ''
    else:
        add_text = 'Temp'
        
    input_df = input_df.set_index('Date')
    
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
    flux_df.index.name = 'Date'
    flux_df = flux_df.reset_index()
    flux_df['Phase'] = flux_df['Date'].apply(lambda x: input_df.loc[input_df.index == x, 'Phase'].iloc[0])
    flux_df.to_csv("OUTPUT_DateWiseSNAppFlux" + add_text, sep=" ", index=True, na_rep="INDEF")

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 2:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if len(band) == 1:
                    dict_flux[day][dict_centralwav[band]] = flux

    bolm_df = pd.DataFrame(dict_flux).T
    bolm_df.index.name = 'Date'

    bolm_df['JD'] = bolm_df.apply(lambda x: cald_to_jd(x.name), axis=1)
    bolm_df = bolm_df.reset_index().set_index(keys='JD')
    bolm_df = bolm_df.interpolate(method='linear', limit=2)
    bolm_df = bolm_df.set_index(keys='Date', drop=True).T

    dict_flux = {}

    for caldate in bolm_df.columns.values[1:]:
        series = bolm_df[caldate].dropna().apply(lambda x: float(x))
#         spline = CubicSpline(series.index.values.tolist(), series.values.tolist(), bc_type='natural', extrapolate=True)
        spline = Rbf(series.index.values.tolist(), series.values.tolist())
        
        if cald_to_jd(caldate) - date_explosion > 100:
            wave_data = np.linspace(4000, 9200, 1000)
        else:
            wave_data = np.linspace(3100, 9200, 1000)

        netflux = np.trapz(spline(wave_data), wave_data)
        lumerr = calc_lum(netflux)
        
        dict_flux[caldate] = {}
        dict_flux[caldate]['JD'] = cald_to_jd(caldate)
        dict_flux[caldate]['Phase'] = dict_flux[caldate]['JD'] - date_explosion
        dict_flux[caldate]['Flux'] = netflux
        dict_flux[caldate]['Lum'] = lumerr[0]
        dict_flux[caldate]['Err'] = lumerr[1]

        if plot:
            fig_temp = plt.figure(figsize=(8, 6))
            ax = fig_temp.add_subplot(111)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wave_data, spline(wave_data), 'r-', label='Rbf Fit')
#             ax.plot(wave_data, spline(wave_data), 'r-', label='CubicSpline Fit')
            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig_temp)

    fbolm_df = pd.DataFrame(dict_flux).T
    fbolm_df.index.name = 'Date'
    fbolm_df.to_csv("OUTPUT_DateWiseSNBolFlux" + add_text, sep=" ", index=True)

    for index, row in fbolm_df.iterrows():
        if row['Phase'] > 150:
            fbolm_df.loc[index, 'MNi'] = get_nickelmass(row['Lum'], phase=row['Phase'])[0]
            fbolm_df.loc[index, 'MNi2'] = get_nickelmass(row['Lum'], phase=row['Phase'])[1]
            fbolm_df.loc[index, 'MNiErr'] = get_nickelmass(row['Err'], phase=row['Phase'])[0]

    return fbolm_df


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
        zp = dict_zp[band]
        rlambda = dict_R[band]
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
    output_df.to_csv("OUTPUT_NetSNMag" + add_text, sep=" ", index=False)
    
    for index, band in output_df['FILTER'].iteritems():
        magflux = calc_absmagflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'], band=band)
        output_df.loc[index, 'AbsMag'] = magflux[0]
        output_df.loc[index, 'AbsErr'] = magflux[1]
        output_df.loc[index, 'Flux'] = magflux[2]
        output_df.loc[index, 'FluxErr'] = magflux[3]

    output_df['ALambda'] = output_df['FILTER'].apply(lambda x: float(fmt_flt.format(dict_R[x] * EBV_mag)))
    output_df.to_csv("OUTPUT_NetSNFlux" + add_text, sep=" ", index=False)
    
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
raw_df = calc_fluxdf("OUTPUT_FinalSNMag", tempsub=False)
output_df = calc_fluxdf("OUTPUT_FinalSNMagTemp", tempsub=True)

output_df = output_df[~((output_df['Date'] == '2017-01-13') & (output_df['FILTER'] == 'U'))]
output_df = output_df[~((output_df['Date'] == '2017-02-12') & (output_df['FILTER'] == 'B'))]
output_df = output_df[~((output_df['Date'] == '2017-02-18') & (output_df['FILTER'] == 'B'))]
output_df = output_df[~((output_df['Date'] == '2017-02-09') & (output_df['FILTER'] == 'I'))]

max_epoch = output_df['Phase'].max()

main_df = output_df[output_df['Phase'] < 400].copy()
mainraw_df = raw_df[raw_df['Phase'] < 400].copy()

fbolm_df = organise_dataframe(output_df)
vabs_df = output_df[output_df['FILTER'] == 'V'].copy()
nebular_df = output_df[(output_df['Phase'] > 110) & (output_df['Phase'] < 400)].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Apparent Magnitude Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig_app = plt.figure(figsize=(8, 5))
ax_main = fig_app.add_subplot(111)

for band, band_df in main_df.groupby('FILTER'):
    ax_main.scatter(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], marker=dict_markers[band][1],
                    c=dict_markers[band][2], s=10, label=dict_markers[band][3])
    ax_main.errorbar(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], yerr=band_df['FERR'],
                     c=dict_markers[band][2], fmt='', linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

handles, labels = ax_main.get_legend_handles_labels()
handles = [handles[3], handles[0], handles[4], handles[2], handles[1]]
labels = [labels[3], labels[0], labels[4], labels[2], labels[1]]
ax_main.legend(handles, labels, fontsize=11, markerscale=2, loc=1, frameon=False)

ax_main.set_ylim(24.5, 13.9)
ax_main.set_xlim(-15, max_epoch + 20)
ax_main.yaxis.set_ticks_position('both')
ax_main.xaxis.set_ticks_position('both')
ax_main.yaxis.set_major_locator(MultipleLocator(2))
ax_main.yaxis.set_minor_locator(MultipleLocator(0.5))
ax_main.xaxis.set_major_locator(MultipleLocator(50))
ax_main.xaxis.set_minor_locator(MultipleLocator(10))
ax_main.set_ylabel('Apparent Magnitude [mag]', fontsize=10)
ax_main.set_xlabel('Time Since Explosion [Days]', fontsize=10)
ax_main.tick_params(which='both', direction='in', width=1, labelsize=10)

fig_app.savefig("OUTPUT_PlotApparentLC.eps", format="eps", dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig_app)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Comparison Of The Template Subtracted Magnitudes With Original Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
fig_comp = plt.figure(figsize=(9, 6))
ax_comp = fig_comp.add_subplot(111)

for band, band_df in main_df.groupby('FILTER'):
    ax_comp.scatter(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], marker=dict_markers[band][1],
                    c=dict_markers[band][2], s=10, label=dict_markers[band][3])
    ax_comp.errorbar(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], yerr=band_df['FERR'],
                     c=dict_markers[band][2], fmt='', linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

for band, band_df in mainraw_df.groupby('FILTER'):
    ax_comp.scatter(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], marker='+',
                    c='lightsalmon', s=10, label=None)
    ax_comp.errorbar(band_df['Phase'], band_df['FMAG'] + dict_markers[band][0], yerr=band_df['FERR'],
                     c='lightsalmon', fmt='', linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

ax_comp.legend(markerscale=2, frameon=False)
ax_comp.set_ylim(24, 13.9)
ax_comp.set_xlim(-10, max_epoch + 20)
ax_comp.yaxis.set_ticks_position('both')
ax_comp.xaxis.set_ticks_position('both')
ax_comp.yaxis.set_major_locator(MultipleLocator(2))
ax_comp.yaxis.set_minor_locator(MultipleLocator(0.5))
ax_comp.xaxis.set_major_locator(MultipleLocator(50))
ax_comp.xaxis.set_minor_locator(MultipleLocator(10))
ax_comp.set_ylabel('Apparent Magnitude [mag]', fontsize=12)
ax_comp.set_xlabel('Time Since Explosion [Days]', fontsize=12)
ax_comp.tick_params(which='both', direction='in', width=1, labelsize=12)

fig_comp.savefig("OUTPUT_PlotTempCompLC.eps", format="eps", dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig_comp)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The V-Band Absolute Magnitude Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig_vabs = plt.figure(figsize=(8, 7))
ax_vabs = fig_vabs.add_subplot(111)

for file_name in list_files:
    name  = file_name.split('/')[-1].split('.')[0]
    if name!= '2005cs':
        data_df = pd.read_csv(file_name, sep="\s+", comment='#', engine='python')
        data_df = coltorow_df(name, data_df)
        temp_df = data_df[data_df['FILTER'] == 'V'].copy()
        temp_df = temp_df[temp_df['Phase'] < 210].sort_values(by='Phase')

        ax_vabs.plot(temp_df['Phase'], temp_df['AbsMag'], linestyle=':', markersize=4, marker=dict_snmark[name][0], 
                        c=dict_snmark[name][1], label=name)
    #     ax_vabs.errorbar(temp_df['Phase'], temp_df['AbsMag'], yerr=temp_df['AbsErr'], linestyle='', linewidth=0.5,
    #                      marker=dict_snmark[name][0], c=dict_snmark[name][1], capsize=1, capthick=1, label=None)

ax_vabs.plot(vabs_df['Phase'], vabs_df['AbsMag'], color='k', markersize=5, marker='o', linestyle='-', label=name_SNe)
ax_vabs.errorbar(vabs_df['Phase'], vabs_df['AbsMag'], yerr=vabs_df['AbsErr'], color='k', linestyle='', linewidth=0.5,
                 capsize=1, capthick=1, label=None)

ax_vabs.legend(fontsize=12, markerscale=2, loc=1, frameon=False)
ax_vabs.set_ylim(-11.8, -18)
ax_vabs.set_xlim(-5, max_epoch + 20)
# ax_vabs.axvline(x=0, ymin=-17, ymax=-16, linestyle='--', color='k')
# ax_vabs.axvline(x=40, ymin=-17, ymax=-16, linestyle='--', color='k')
# ax_vabs.axvline(x=94, ymin=-17, ymax=-16, linestyle='--', color='k')
# ax_vabs.text(0, ax.get_ybound()[0] + 0.5, r'$\rm t_{exp}$', rotation=-90, color='k')
# ax_vabs.text(45, ax.get_ybound()[0] + 0.5, r'$\rm t_{tran}$', rotation=-90, color='k')
# ax_vabs.text(90, ax.get_ybound()[0] + 0.5, r'$\rm t_{end}$', rotation=-90, color='k')
ax_vabs.yaxis.set_ticks_position('both')
ax_vabs.xaxis.set_ticks_position('both')
ax_vabs.yaxis.set_major_locator(MultipleLocator(1))
ax_vabs.yaxis.set_minor_locator(MultipleLocator(0.25))
ax_vabs.xaxis.set_major_locator(MultipleLocator(50))
ax_vabs.xaxis.set_minor_locator(MultipleLocator(10))
ax_vabs.set_ylabel(r'Absolute V-Band Magnitude [$\rm M_V$]', fontsize=14)
ax_vabs.set_xlabel('Time Since Explosion [Days]', fontsize=14)
ax_vabs.tick_params(which='both', direction='in', width=1, labelsize=14)

fig_vabs.savefig("OUTPUT_PlotVAbsLC.eps", format="eps", dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig_vabs)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig_bol = plt.figure(figsize=(8, 8))
ax_bol = fig_bol.add_subplot(111)

for file_name in list_files:
    name  = file_name.split('/')[-1].split('.')[0]
    data_df = pd.read_csv(file_name, sep="\s+", comment='#', engine='python')
    data_df = obtain_epochwisedf(coltorow_df(name, data_df))
    temp_df = calc_boldf(name, data_df)
    temp_df = temp_df[temp_df['Phase'] < max_epoch].sort_values(by='Phase')
    
    ax_bol.semilogy(temp_df['Phase'], temp_df['Lum'], linestyle=':', markersize=4, marker=dict_snmark[name][0], 
                    c=dict_snmark[name][1], label=name)

ax_bol.semilogy(fbolm_df['Phase'], fbolm_df['Lum'], linestyle='-', color='k', marker='o', markersize=5, label=name_SNe)
# ax_bol.errorbar(fbolm_df['Phase'], fbolm_df['Lum'], yerr=fbolm_df['Err'], linestyle='', color='k', markersize=5, 
#                 capsize=2, capthick=1, label=None)

ax_bol.legend(fontsize=12, markerscale=3, loc=1, frameon=False)
ax_bol.set_ylim(1e39, 7e42)
ax_bol.set_xlim(-10, max_epoch + 20)
ax_bol.yaxis.set_ticks_position('both')
ax_bol.xaxis.set_ticks_position('both')
ax_bol.xaxis.set_major_locator(MultipleLocator(100))
ax_bol.xaxis.set_minor_locator(MultipleLocator(20))
ax_bol.set_ylabel(r'Quasi-Bolometric Luminosity [$\rm erg\ s^{-1}$]', fontsize=14)
ax_bol.set_xlabel('Time Since Explosion [Days]', fontsize=14)
ax_bol.tick_params(which='both', direction='in', width=1, labelsize=14)

fig_bol.savefig("OUTPUT_PlotBolometricLC.eps", format="eps", dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig_bol)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Absolute V-Band Magnitude Vs Nickel Mass At Mid-Plateau
# ------------------------------------------------------------------------------------------------------------------- #
data_sn = pd.read_csv(DIR_SNe + 'Param_Data/Log_MNi', sep='\s+', comment='#', engine='python')
data_sn = data_sn.replace('INDEF', np.nan)
data_sn = data_sn.set_index('Name', drop=True)
data_sn = data_sn.astype('float64')

data_sn['logMNi'] = data_sn['MNi'].apply(lambda x: np.log10(x))
data_sn['logMNiErr+'] = data_sn['MNiErr+'] / data_sn['MNi']
data_sn['logMNiErr-'] = data_sn['MNiErr-'] / data_sn['MNi']

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)

ax.set_xscale('log')
ax.set_ylim(-12.3, -19)
ax.plot([0, 1], [0, 1], transform=ax.transAxes, linestyle='--', color='g')
ax.errorbar(0.029, -16.9, xerr=0.005, yerr=0.2, color='r', fmt='*', markersize=12, capsize=3, label=name_SNe)
ax.errorbar(data_sn['MNi'], data_sn['Mv'], yerr=data_sn['MvErr'], color='k', fmt='o', capsize=3, markersize=5, 
            capthick=0.5, xerr=np.vstack((data_sn['MNiErr-'], data_sn['MNiErr+'])), elinewidth=1, label=None)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
ax.legend(handles, labels, fontsize=12, markerscale=2, frameon=False)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.5))
ax.set_ylabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}$', fontsize=16)
ax.set_xlabel(r'$\rm M_{Ni}\ [M_{\odot}]$', fontsize=16)
ax.tick_params(which='both', direction='in', width=1, labelsize=14)

fig.savefig('OUTPUT_PlotNi.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Mid-Plateau Velocity Vs Nickel Mass (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)

ax.set_xscale('log')
ax.set_ylim(900, 8400)
ax.errorbar(0.029, 4910, xerr=0.005, yerr=100, color='r', fmt='*', markersize=12, capsize=3, label=name_SNe)
ax.errorbar(data_sn['MNi'], data_sn['v50'], yerr=data_sn['v50Err'], color='k', fmt='o', capsize=3, markersize=5, 
            xerr=np.vstack((data_sn['MNiErr-'], data_sn['MNiErr+'])), elinewidth=1, capthick=0.5, label=None)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
ax.legend(handles, labels, fontsize=12, markerscale=2, frameon=False)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(2000))
ax.yaxis.set_minor_locator(MultipleLocator(500))
ax.set_ylabel(r'Mid-Plateau Velocity, $\rm V_{50}$', fontsize=16)
ax.set_xlabel(r'$\rm M_{Ni}\ [M_{\odot}]$', fontsize=16)
ax.tick_params(which='both', direction='in', width=1, labelsize=14)

fig.savefig('OUTPUT_PlotV50.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Absolute V-Band Magnitude Vs Mid-Plateau Velocity (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)

ax.errorbar(4910, -16.9, xerr=100, yerr=0.2, color='r', fmt='*', markersize=12, capsize=3, label=name_SNe)
ax.errorbar(data_sn['v50'], data_sn['Mv'], yerr=data_sn['MvErr'], xerr=data_sn['v50Err'], color='k', markersize=5,
            fmt='o', capsize=3, elinewidth=1, capthick=0.5, label=None)

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] for h in handles]
ax.legend(handles, labels, fontsize=12, markerscale=2, frameon=False)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.5))
ax.xaxis.set_major_locator(MultipleLocator(2000))
ax.xaxis.set_minor_locator(MultipleLocator(500))
ax.set_ylabel(r'Absolute Mid-Plateau Magnitude, $\rm M^{50}_{V}$', fontsize=16)
ax.set_xlabel(r'Mid-Plateau Velocity, $\rm V_{50}$', fontsize=16)
ax.tick_params(which='both', direction='in', width=1, labelsize=14)

fig.savefig('OUTPUT_PlotMv.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #

