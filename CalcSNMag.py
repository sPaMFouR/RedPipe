#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxx---------Calculate SN Magnitudes For All The Epochs Using Secondary Standard Magnitudes---------xxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
from pyraf import iraf
from astropy.io import fits

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
date_explosion = 2456841.50
filters = ['U', 'B', 'V', 'R', 'I']
list_psfcol = ["ID", "IMAGE", "IFILTER", "XCENTER", "YCENTER", "SKY_COUNTS", "AIRMASS", "PSFRAD", "PSFMAG", "PSFERR"]
list_magcol = ["ID", "IMAGE", "IFILTER", "XCENTER", "YCENTER", "SKY_COUNTS", "AIRMASS",
               "APER_1", "APER_2", "MAG_1", "MAG_2", "MERR_1", "MERR_2"]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
JD_keyword = 'JD'
date_keyword = 'DATE-OBS'
filter_keyword = 'FILTER'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.ptools(_doprint=0)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Handling Files & Lists
# ------------------------------------------------------------------------------------------------------------------- #

def remove_file(file_name):
    """
    Removes the file "file_name" in the constituent directory.
    Args:
         file_name  : Name of the file to be removed from the current directory
    Returns:
        None
    """
    try:
        os.remove(file_name)
    except OSError:
        pass


def remove_similar_files(common_text):
    """
    Removes similar files based on the string "common_text".
    Args:
        common_text : String containing partial name of the files to be deleted
    Returns:
        None
    """
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)


# noinspection PyShadowingNames
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


def list_statistics(list_values):
    """
    Returns the statistics of the list of elements in the input 'list_values'.
    Args:
        list_values : Input list of elements
    Returns:
        value_mean  : Mean of the list of elements
        value_median: Median of the list of elements
        value_std   : Standard Deviation of the list of elements
    """
    value_mean = np.mean(list_values)
    value_median = np.median(list_values)
    value_std = np.std(list_values)

    return value_mean, value_median, value_std


def reject(list_values, iterations=2):
    """
    Rejects outliers from the input 'list_values'.
    Args:
        list_values : Input list of elements
        iterations  : No. of iterations of rejection to be run on the input list
    Returns:
        list_reject : Output list after rejecting outliers from the input 'list_values'
    """
    list_reject = filter(lambda x: x != 'INDEF', list_values)
    list_reject = [round(float(val), int(precision)) for val in list_reject]
    list_reject.sort()

    for _ in range(0, iterations):
        if len(list_values) > 2:
            value_mean, value_median, value_std = list_statistics(list_reject)

            if abs(list_reject[0] - value_median) < abs(list_reject[-1] - value_median):
                remove_index = -1
            else:
                remove_index = 0

            if abs(list_reject[remove_index] - value_median) > value_std:
                list_reject.pop(remove_index)

    return list_reject


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
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def txdump(common_text, output_file):
    """
    Performs TXDUMP task on the MAG or ALS files generated by photometry tasks. This extracts
    useful data from magnitude files.
    Args:
        common_text : Partial name of the MAG or ALS files from which data is to be extracted
        output_file : Output file where data from the list of input files is to be written
    Returns:
        None
    """
    if re.search('mag', common_text):
        fields = "ID, IMAGE, IFILTER, XCENTER, YCENTER, MSKY, XAIRMASS, RAPERT, MAG, MERR"
    else:
        fields = "ID, IMAGE, IFILTER, XCENTER, YCENTER, MSKY, XAIRMASS, PSFRAD, MAG, MERR"

    task = iraf.noao.digiphot.ptools.txdump
    task.unlearn()

    file_temp = 'temp_dump'
    group_similar_files(str(file_temp), common_text=common_text)
    task(textfile='@' + str(file_temp), fields=fields, expr="yes", Stdout=str(output_file))
    remove_file(file_temp)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Calculate JD, Date, FILTER From The Image Header
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_jd(file_name):
    """
    Extracts Julian day, date of observation, and filter from the header of the FITS file.
    Args:
        file_name   : Name of the FITS file from which header details have to be extracted
    Returns:
        jd          : Julian Day of observation of the FITS file
        date        : Date of observation of the FITS file
        ifilter     : Band of observation of the FITS file
    """
    file_header = fits.getheader(file_name)

    jd = file_header[JD_keyword]
    date = file_header[date_keyword]
    ifilter = file_header[filter_keyword][-1]

    return jd, date, ifilter

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Manipulating Pandas DataFrames
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


def append_missing_data(input_df):
    """
    Appends missing data for a filter as a column of 'INDEF' to the DataFrame.
    Args:
        input_df    : Pandas DataFrame containing star magnitudes
    Returns:
        output_df   : Pandas DataFrame containing appended columns for missing data
    """
    star_id = set(input_df.index.values)

    for band in filters:
        if band not in set(input_df['FILTER'].values):
            data_ext = [[str(band)] + ['INDEF'] * (len(input_df.columns.values) - 1) for _ in range(0, len(star_id))]
            input_df = pd.concat([pd.DataFrame(data_ext, columns=input_df.columns.values, index=star_id), input_df])

    output_df = input_df.sort_values(by='FILTER').sort_index(kind='mergesort')
    output_df = output_df.replace('INDEF', np.nan, regex=True)

    return output_df


def colormag_to_ubvriframe(input_df, err=False):
    """
    Creates a pandas DataFrame with magnitudes and color terms from an input DataFrame with V-band magnitude and
    color terms.
    Args:
        input_df    : Pandas DataFrame containing standard star magnitudes and color terms with errors
        err         : Boolean specifying whether the DataFrame contains error data
    Returns:
        output_df   : Pandas DataFrame containing V-band magnitude and color terms
    """
    output_df = pd.DataFrame(index=input_df.index.values)

    output_df['U'] = add_series([input_df['V'], input_df['B-V'], input_df['U-B']], err=err)
    output_df['B'] = add_series([input_df['V'], input_df['B-V']], err=err)
    output_df['V'] = input_df['V']
    output_df['R'] = add_series([input_df['V'], input_df['V-R']], sub=True, err=err)
    output_df['I'] = add_series([input_df['V'], input_df['V-I']], sub=True, err=err)

    output_df = output_df.round(int(precision))
    output_df = output_df.replace(np.nan, 'INDEF', regex=True)

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Calculate (PSF Magnitudes, PSF Correction) & (PHOT Magnitudes, Aperture Correction)
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_psfmag(file_als, zpmag_df, zperr_df):
    """
    Calculates true magnitudes of the SN from Txs'dumped ALS & MAG files obtained after photometry.
    Arg
        file_als    : File containing Tx'dumped magnitudes obtained from ALS files generated by IRAF
        zpmag_df    : Pandas DataFrame containing zero-point magnitudes
        zperr_df    : Pandas DataFrame containing zero-point errors
    Returns:
        mag_df      : Pandas DataFrame containing broadband magnitudes
        err_df      : Pandas DataFrame containing errors in magnitudes
    """
    psf_df = pd.read_csv(filepath_or_buffer=file_als, sep='\s+', names=list_psfcol, index_col=0, engine='python')
    psf_df['FILTER'] = psf_df['IFILTER'].apply(lambda x: str(x)[-1])
    psf_df = psf_df[['FILTER', 'PSFMAG', 'PSFERR']]
    psf_df = psf_df.sort_index().sort_values(by='FILTER', kind='mergesort')

    indexes = psf_df.index.values
    star_count = len(set(indexes)) - 1

    psfstars_df = psf_df.loc[psf_df.index != star_count + 1].copy()
    psfsn_df = psf_df.loc[psf_df.index == star_count + 1].copy()

    file_mag = file_als[:-4] + 'mag4'
    mag_df = pd.read_csv(filepath_or_buffer=file_mag, sep='\s+', names=list_magcol, index_col=0, engine='python')

    psfstars_df['MAG_2'] = mag_df['MAG_2']
    psfstars_df['PSFCOR'] = psfstars_df['PSFMAG'] - psfstars_df['MAG_2']

    data_grouped = psfstars_df[['PSFCOR', 'FILTER']].groupby(['FILTER'])
    mean = {}
    stdev = {}

    for band in set(psfstars_df['FILTER'].values):
        temp_list = reject(data_grouped.get_group(name=band)['PSFCOR'].tolist(), iterations=int(star_count / 3) + 1)
        mean[band] = np.mean(temp_list)
        stdev[band] = np.std(temp_list)

    psfsn_df['COR_MEAN'] = psfsn_df['FILTER'].apply(lambda x: mean[x])
    psfsn_df['COR_STD'] = psfsn_df['FILTER'].apply(lambda x: stdev[x])

    date = file_als.split("_")[1]
    psfsn_df['ZP_MAG'] = psfsn_df['FILTER'].apply(lambda x: zpmag_df.loc[date, x])
    psfsn_df['ZP_ERR'] = psfsn_df['FILTER'].apply(lambda x: zperr_df.loc[date, x])

    psfsn_df['FMAG'] = psfsn_df['PSFMAG'] - psfsn_df['COR_MEAN'] - psfsn_df['ZP_MAG']
    psfsn_df['FERR'] = add_series([psfsn_df['PSFERR'], psfsn_df['COR_STD'], psfsn_df['ZP_ERR']], err=True)
    psfsn_df = psfsn_df.round(int(precision))
    psfsn_df = append_missing_data(psfsn_df[['FILTER', 'FMAG', 'FERR']]).set_index('FILTER')

    psfsn_df.to_csv("OUTPUT_SNMag_" + str(date), sep=" ", index=False)

    return date, psfsn_df


def calculate_tempmag(file_mag, zpmag_df, zperr_df):
    """
    Calculates true magnitudes of the SN from Tx'dumped MAG files obtained after template subtraction photometry.
    Arg
        file_mag    : File containing Tx'dumped magnitudes obtained from MAG files generated by IRAF
        zpmag_df    : Pandas DataFrame containing zero-point magnitudes
        zperr_df    : Pandas DataFrame containing zero-point errors
    Returns:
        mag_df      : Pandas DataFrame containing broadband magnitudes
        err_df      : Pandas DataFrame containing errors in magnitudes
    """
    temp_df = pd.read_csv(filepath_or_buffer=file_mag, sep="\s+", names=list_magcol, index_col=0, engine='python')
    temp_df['FILTER'] = temp_df['IFILTER'].apply(lambda x: str(x)[-1])
    temp_df = temp_df[['FILTER', 'MAG_1', 'MERR_1']]
    temp_df = temp_df.sort_index().sort_values(by='FILTER', kind='mergesort')

    temp_df = temp_df.set_index('FILTER', append=True).replace('INDEF', np.nan)
    temp_df = temp_df.astype('float64').round(3).reset_index(level=['FILTER'])

    indexes = temp_df.index.values
    star_count = len(set(indexes))

    file_mag4 = file_mag[:-4] + 'mag4'
    mag_df = pd.read_csv(filepath_or_buffer=file_mag4, sep="\s+", names=list_magcol, index_col=0, engine='python')
    mag_df['FILTER'] = mag_df['IFILTER'].apply(lambda x: str(x)[-1])
    mag_df['APCOR'] = mag_df['MAG_1'] - mag_df['MAG_2']

    data_grouped = mag_df[['APCOR', 'FILTER']].groupby(['FILTER'])
    mean = {}
    stdev = {}

    for band in set(mag_df['FILTER'].values):
        temp_list = reject(data_grouped.get_group(name=band)['APCOR'].tolist(), iterations=int(star_count / 3) + 1)
        mean[band] = np.mean(temp_list)
        stdev[band] = np.std(temp_list)

    temp_df['COR_MEAN'] = temp_df['FILTER'].apply(lambda x: mean[x])
    temp_df['COR_STD'] = temp_df['FILTER'].apply(lambda x: stdev[x])

    date = file_mag.split("_")[1]
    temp_df['ZP_MAG'] = temp_df['FILTER'].apply(lambda x: zpmag_df.loc[date, x])
    temp_df['ZP_ERR'] = temp_df['FILTER'].apply(lambda x: zperr_df.loc[date, x])

    temp_df['FMAG'] = temp_df['MAG_1'] - temp_df['COR_MEAN'] - temp_df['ZP_MAG']
    temp_df['FERR'] = add_series([temp_df['MERR_1'], temp_df['COR_STD'], temp_df['ZP_ERR']], err=True)
    temp_df = temp_df.round(int(precision))
    temp_df = append_missing_data(temp_df[['FILTER', 'FMAG', 'FERR']]).set_index('FILTER')

    temp_df.to_csv("OUTPUT_SNMagTemp_" + date, sep=" ", index=False)

    return date, temp_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Zero Point Correction For Different Epochs
# ------------------------------------------------------------------------------------------------------------------- #
file_mag = 'OUTPUT_truestdmag'
file_err = 'OUTPUT_truestderr'

truemag_df = pd.read_csv(filepath_or_buffer=file_mag, sep="\s+", index_col=0, engine='python')
trueerr_df = pd.read_csv(filepath_or_buffer=file_err, sep="\s+", index_col=0, engine='python')

dict_zpmag = {}
dict_zperr = {}
list_instrmag = group_similar_files("", "OUTPUT_instrmag_*")
list_instrerr = group_similar_files("", "OUTPUT_instrerr_*")

for index, file_name in enumerate(list_instrmag):
    obsmag_df = pd.read_csv(filepath_or_buffer=file_name, sep="\s+", index_col=0, engine='python')
    obserr_df = pd.read_csv(filepath_or_buffer=list_instrerr[index], sep="\s+", index_col=0, engine='python')

    obsmag_df = obsmag_df.replace('INDEF', np.nan).round(decimals=int(precision))
    obserr_df = obserr_df.replace('INDEF', np.nan).round(decimals=int(precision))

    dict_zpmag[file_name.split("_")[-1]] = {}
    dict_zperr[file_name.split("_")[-1]] = {}

    for column in truemag_df:
        temp_list = reject((obsmag_df[column] - truemag_df[column]).dropna().tolist(), iterations=4)

        if len(temp_list) != 0:
            zpmag = np.mean(temp_list, dtype='float64')
            zperr = (np.std(temp_list) ** 2 + trueerr_df[column].mean() ** 2 + obserr_df[column].mean() ** 2) ** 0.5
        else:
            zpmag = np.nan
            zperr = np.nan

        dict_zpmag[file_name.split("_")[-1]][column] = round(zpmag, int(precision))
        dict_zperr[file_name.split("_")[-1]][column] = round(zperr, int(precision))

zpmag_df = pd.DataFrame(data=dict_zpmag).T
zperr_df = pd.DataFrame(data=dict_zperr).T

display_text("Zero Point Correction For Different Epochs Of Observation Have Been Computed")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate DataFrame With Date Of Observation Related To JD And Filter
# ------------------------------------------------------------------------------------------------------------------- #
list_fits = group_similar_files("", "ca_*.fits", exceptions="ts_,psf,sub")
dict_JD = {}

for file_name in list_fits:
    JD, date, ifilter = calculate_jd(file_name)
    if date not in dict_JD:
        dict_JD[date] = {}
    dict_JD[date][ifilter] = round(float(JD), int(precision))

JD_df = pd.DataFrame(data=dict_JD)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Final Magnitudes For The SN During Different Epochs
# ------------------------------------------------------------------------------------------------------------------- #
list_psfmag = group_similar_files("", "output_*als1")

finalmag_df = pd.DataFrame()
for file_name in list_psfmag:
    date, psfdate_df = calculate_psfmag(file_name, zpmag_df, zperr_df)
    finalmag_df = pd.concat([finalmag_df, pd.concat([psfdate_df, JD_df.loc[:, date].rename("JD")], axis=1)])

finalmag_df = finalmag_df[['JD', 'FMAG', 'FERR']]
finalmag_df = finalmag_df.dropna().sort_values(by='JD').sort_index(kind='mergesort')
finalmag_df.to_csv("OUTPUT_FinalSNMag", sep=" ", index=True, index_label="FILTER")

display_text("Photometric Magnitudes For The Supernova Have Been Computed")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Final Magnitudes For The SN During Different Epochs (After Template Subtraction)
# ------------------------------------------------------------------------------------------------------------------- #
list_tempmag = group_similar_files("", "output_*mag1")

tempmag_df = pd.DataFrame()
for file_name in list_tempmag:
    date, mag_df = calculate_tempmag(file_name, zpmag_df, zperr_df)
    tempmag_df = pd.concat([tempmag_df, pd.concat([mag_df, JD_df.loc[:, date].rename("JD")], axis=1)])

tempmag_df = tempmag_df[['JD', 'FMAG', 'FERR']]
tempmag_df = tempmag_df.dropna().sort_values(by='JD').sort_index(kind='mergesort')
tempmag_df.to_csv("OUTPUT_FinalSNMagTemp", sep=" ", index=True, index_label="FILTER")

display_text("Template Subtracted Broadband Photometric Magnitudes For The Supernova Have Been Computed")
# ------------------------------------------------------------------------------------------------------------------- #
