#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxx---------Calculate True Magnitudes Of The Secondary Standards In The SN Field----------xxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import numpy as np
import pandas as pd
from pyraf import iraf
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In This Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 4
day_std = '2018-06-24'
object_name = '2018cow'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Site Extinction Coefficients With Errors In Different Bands
# ------------------------------------------------------------------------------------------------------------------- #
# FILTER  EXTINCTION_MEAN   EXTINCTION_ERROR
#   U         0.36              0.07
#   B         0.21              0.04
#   V         0.12              0.04
#   R         0.09              0.04
#   I         0.05              0.03

eeta = {'U': 0.36, 'B': 0.21, 'V': 0.12, 'R': 0.09, 'I': 0.05}
eeta_err = {'U': 0.07, 'B': 0.04, 'V': 0.04, 'R': 0.04, 'I': 0.03}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In Manipulating Pandas DataFrames
# ------------------------------------------------------------------------------------------------------------------- #

# For Generating DataFrame Containing V-Band Magnitude & Color Terms From Landolt Standard Field Matrix
# Standard_Matrix = [[V, B-V, U-B, V-R, R-I, V-I, V_Err, B-V_Err, U-B_Err, V-R_Err, R-I_Err, V-I_Err],
#                    [V, B-V, U-B, V-R, R-I, V-I, V_Err, B-V_Err, U-B_Err, V-R_Err, R-I_Err, V-I_Err],
#                    [V, B-V, U-B, V-R, R-I, V-I, V_Err, B-V_Err, U-B_Err, V-R_Err, R-I_Err, V-I_Err],
#                    [V, B-V, U-B, V-R, R-I, V-I, V_Err, B-V_Err, U-B_Err, V-R_Err, R-I_Err, V-I_Err]]]
# Mag/Color      Error        Indices
#   V            V_Err        (0, 6)
#   B-V          B-V_Err      (1, 7)
#   U-B          U-B_Err      (2, 8)
#   V-R          V-R_Err      (3, 9)
#   R-I          R-I_Err      (4, 10)
#   V-I          V-I_Err      (5, 11)

# For Creating A DataFrame Containing UBVRI Magnitudes From A DataFrame Containing V-Band Magnitudes & Color Terms
# DataFrame_Columns = ['V', 'B-V', 'U-B', 'V-R', 'R-I', 'V-I']
# Colors            Terms               Indexes
#   U           V + B-V + U-B         [0 + 1 + 2]
#   B           V + B-V               [0 + 1]
#   V           V                     [0]
#   R           V - (V-R)             [0 - 3]
#   I           V - (V-I)             [0 - 5]

# For Obtaining Beta Values
# Observed ColorMag DataFrame - ObsDF
# Standard ColorMag DataFrame - StdDF
# List of Standard Alpha Values - Alpha
# Use Standard Alpha Values & Obtain Beta Values

# Beta[0] = (B-V) - Alpha[0] * (B-V)obs
# Beta[1] = (U-B) - Alpha[1] * (U-B)obs
# Beta[2] = (V-R) - Alpha[2] * (V-R)obs
# Beta[3] = (R-I) - Alpha[3] * (R-I)obs
# Beta[4] = (V-I) - Alpha[4] * (V-I)obs
# Beta[5] = V - Vobs - Alpha[5] * (B-V)
# Beta[6] = V - Vobs - Alpha[6] * (V-R)

# Beta[0] = StdDF['B-V'] - Alpha[0] * ObsDF['B-V']
# Beta[1] = StdDF['U-B'] - Alpha[1] * ObsDF['U-B']
# Beta[2] = StdDF['V-R'] - Alpha[2] * ObsDF['V-R']
# Beta[3] = StdDF['R-I'] - Alpha[3] * ObsDF['R-I']
# Beta[4] = StdDF['V-I'] - Alpha[4] * ObsDF['V-I']
# Beta[5] = StdDF['V'] - ObsDF['V'] - Alpha[5] * StdDF['B-V']
# Beta[6] = StdDF['V'] - ObsDF['V'] - Alpha[6] * StdDF['V-R']

list_alpha = [0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]

filters = ['U', 'B', 'V', 'R', 'I']
colors = ['B-V', 'U-B', 'V-R', 'R-I', 'V-I']
color_indices = [tuple(value.split('-')) for value in colors]

mag_terms = ['V'] + colors
mag_indices = [(0, 6), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)]

list_magcol = ['ID', 'IMAGE', 'IFILTER', 'XCENTER', 'YCENTER', 'SKY_COUNTS', 'AIRMASS', 'APER_1', 'APER_2', 'MAG_1',
               'MAG_2', 'ERR_1', 'ERR_2']
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# V-Band Magnitude and Color Terms For Landolt Standards
# ------------------------------------------------------------------------------------------------------------------- #

PG0918 = [[13.327, -0.271, -1.081, -0.129, -0.159, -0.288, 0.0024, 0.0024, 0.0030, 0.0019, 0.0055, 0.0063],
          [14.490,  0.536, -0.032,  0.325,  0.336,  0.661, 0.0033, 0.0058, 0.0095, 0.0039, 0.0076, 0.0085],
          [13.963,  0.765,  0.366,  0.417,  0.370,  0.787, 0.0034, 0.0072, 0.0159, 0.0025, 0.0045, 0.0056],
          [13.537,  0.631,  0.087,  0.367,  0.357,  0.722, 0.0020, 0.0028, 0.0048, 0.0015, 0.0022, 0.0028],
          [12.272,  1.044,  0.821,  0.575,  0.535,  1.108, 0.0021, 0.0030, 0.0071, 0.0016, 0.0018, 0.0018]]

PG0231 = [[16.105, -0.329, -1.192, -0.162, -0.371, -0.534, 0.0068, 0.0083, 0.0045, 0.0276, 0.1066, 0.1221],
          [12.772,  0.710,  0.270,  0.405,  0.394,  0.799, 0.0008, 0.0015, 0.0030, 0.0011, 0.0030, 0.0030],
          [14.735,  1.448,  1.342,  0.954,  0.998,  1.951, 0.0030, 0.0072, 0.0178, 0.0034, 0.0026, 0.0057],
          [13.702,  0.671,  0.114,  0.399,  0.385,  0.783, 0.0014, 0.0078, 0.0148, 0.0028, 0.0064, 0.0085],
          [14.027,  1.088,  1.046,  0.675,  0.586,  1.256, 0.0029, 0.0075, 0.0312, 0.0081, 0.0064, 0.0110],
          [13.804,  0.677,  0.201,  0.390,  0.369,  0.757, 0.0046, 0.0040, 0.0075, 0.0035, 0.0017, 0.0023]]

PG0942 = [[14.004, -0.294, -1.175, -0.130, -0.149, -0.280, 0.0045, 0.0056, 0.0069, 0.0069, 0.0120, 0.0144],
          [14.731,  0.783,  0.339,  0.610,  0.477,  1.081, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
          [14.108,  0.525,  0.085,  0.368,  0.333,  0.697, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
          [14.989,  0.727,  0.369,  0.539,  0.376,  0.909, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
          [13.707,  0.564,  0.129,  0.348,  0.343,  0.687, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042]]

PG2213 = [[14.124, -0.217, -1.125, -0.092, -0.110, -0.203, 0.0022, 0.0028, 0.0063, 0.0044, 0.0085, 0.0092],
          [14.178,  0.673,  0.100,  0.406,  0.403,  0.808, 0.0050, 0.0033, 0.0057, 0.0030, 0.0050, 0.0060],
          [12.706,  0.749,  0.297,  0.427,  0.402,  0.829, 0.0011, 0.0023, 0.0026, 0.0008, 0.0015, 0.0015],
          [15.109,  0.721,  0.177,  0.426,  0.404,  0.830, 0.0045, 0.0057, 0.0068, 0.0023, 0.0068, 0.0064]]

PG1657 = [[15.015, -0.149, -0.940, -0.063, -0.033, -0.100, 0.0067, 0.0053, 0.0090, 0.0087, 0.0270, 0.0330],
          [14.033,  1.069,  0.730,  0.573,  0.539,  1.113, 0.0007, 0.0007, 0.0064, 0.0064, 0.0057, 0.0127],
          [14.721,  0.708,  0.065,  0.417,  0.420,  0.838, 0.0021, 0.0064, 0.0071, 0.0014, 0.0000, 0.0014],
          [15.225,  0.840,  0.385,  0.521,  0.444,  0.967, 0.0000, 0.0042, 0.0085, 0.0057, 0.0127, 0.0071]]

PG1633 = [[14.397, -0.192, -0.974, -0.093, -0.116, -0.212, 0.0025, 0.0022, 0.0047, 0.0033, 0.0089, 0.0111], 
          [15.256,  0.873,  0.320,  0.505,  0.511,  1.015, 0.0036, 0.0052, 0.0090, 0.0036, 0.0093, 0.0111],
          [12.969,  1.081,  1.007,  0.590,  0.502,  1.090, 0.0017, 0.0020, 0.0069, 0.0012, 0.0014, 0.0020],
          [13.229,  1.134,  1.138,  0.618,  0.523,  1.138, 0.0025, 0.0022, 0.0038, 0.0016, 0.0022, 0.0038],
          [13.691,  0.535, -0.025,  0.324,  0.327,  0.650, 0.0020, 0.0020, 0.0050, 0.0017, 0.0033, 0.0033]]

PG2331 = [[15.182, -0.066, -0.487, -0.012, -0.031, -0.044, 0.0057, 0.0071, 0.0035, 0.0078, 0.0057, 0.0127],
          [13.051,  0.741,  0.257,  0.419,  0.401,  0.821, 0.0021, 0.0014, 0.0014, 0.0014, 0.0014, 0.0014],
          [14.744,  0.819,  0.429,  0.481,  0.454,  0.935, 0.0035, 0.0007, 0.0014, 0.0035, 0.0064, 0.0021]]

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
    Removes the file 'file_name' in the constituent directory.
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
    Removes similar files based on the string 'common_text'.
    Args:
        common_text : String containing partial name of the files to be deleted
    Returns:
        None
    """
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)


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
    list_reject = map(float, list_reject)
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


def reject_series(input_series):
    """
    Rejects outliers from the input 'list_values'.
    Args:
        input_series  : Input Pandas Series of elements
    Returns:
        output_series : Modified Pandas Series after rejecting outliers from the input 'input_series'
    """
    input_series = input_series.replace('INDEF', np.nan).dropna().astype('float64')
    input_series = input_series.sort_values(ascending=True)

    sigma = 1
    output_series = pd.Series()
    while len(output_series.index) == 0:
        output_series = input_series[(input_series - input_series.median()).abs() < sigma * input_series.std()]
        sigma += 0.5

    return output_series


def file_is_empty(path):
    """
    Checks if a file is empty or not.
    Args:
        path      : Location of the file to be checked
    Returns:
        bool_file : Whether the file is empty or not
    """
    bool_file = (os.path.getsize(path) == 0)
    return bool_file


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
    task(textfile='@' + str(file_temp), fields=fields, expr='yes', Stdout=str(output_file))
    remove_file(file_temp)

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


def mul_series(input_series, mult_num):
    """
    Multiply the input series by a number (mult_num) if and only if the value is not an 'INDEF'.
    Args:
        input_series    : Input Pandas series on which multiplication operation is to be performed
        mult_num        : Number to multiplied to the input series
    Returns:
        output_series   : Output Pandas series which is a result of multiplication operation
    """
    output_series = input_series.apply(lambda val: val * float(mult_num) if val != 'INDEF' else val)

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
    output_df = output_df.replace('INDEF', np.nan)

    return output_df


def unorgmag_to_ubvriframe(input_df):
    """
    Creates a pandas DataFrame with broadband magnitudes from an input DataFrame with unorganised magnitudes.
    Args: 
        input_df    : Pandas DataFrame containing magnitudes and color terms
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes
    """
    dict_stars = {}
    for index, row in input_df.iterrows():
        if index not in dict_stars.keys():
            dict_stars[index] = {}
        dict_stars[index][row[0]] = row[1]

    output_df = pd.DataFrame(data=dict_stars).T[filters]
    output_df = output_df.apply(pd.to_numeric, errors='coerce').round(int(precision))
    output_df = output_df.replace(np.nan, 'INDEF')

    return output_df


def unorgmag_to_colormagframe(input_df, err=False):
    """
    Creates a pandas DataFrame with magnitudes and color terms from an input DataFrame with unorganised magnitudes.
    Args: 
        input_df    : Pandas DataFrame containing magnitudes and color terms
        err         : True, if the DataFrame contains error data
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes and color terms
    """
    dict_stars = {}
    for index, row in input_df.iterrows():
        if index not in dict_stars.keys():
            dict_stars[index] = {}
        dict_stars[index][row[0]] = row[1]

    for value in dict_stars.values():
        for index1, index2 in color_indices:
            if value[index1] != 'INDEF' and value[index2] != 'INDEF':
                if not err:
                    append_val = float(value[index1]) - float(value[index2])
                else:
                    append_val = (float(value[index1]) ** 2 + float(value[index2]) ** 2) ** 0.5
            else:
                append_val = 'INDEF'

            value[index1 + '-' + index2] = append_val

    output_df = pd.DataFrame(data=dict_stars).T[mag_terms]
    output_df = output_df.apply(pd.to_numeric, errors='coerce').round(int(precision))

    return output_df


def matrix_to_colormagframe(standard_matrix):
    """
    Creates a pair of pandas DataFrame with V-band magnitude and color terms from an input matrix containing 
    Landolt standard field data.
    Args:
        standard_matrix : Python list of lists which contain standard star magnitudes and color terms with errors
    Returns:
        mag_df          : Pandas DataFrame containing magnitude data
        err_df          : Pandas DataFrame containing error data
    """
    input_df = pd.DataFrame(data=standard_matrix, index=range(1, len(standard_matrix) + 1))
    mag_df = pd.DataFrame()
    err_df = pd.DataFrame()

    for (mag, err) in mag_indices:
        mag_df[mag_terms[mag]] = input_df[mag]
        err_df[mag_terms[mag]] = input_df[err]

    return mag_df, err_df


def ubvrimag_to_colormagframe(input_df, err=False):
    """
    Creates a pandas DataFrame with magnitudes and color terms from an input DataFrame with unorganised magnitudes.
    Args: 
        input_df    : Pandas DataFrame containing magnitudes and color terms
        err         : True, if the DataFrame contains error data
    Returns:
        output_df   : Pandas DataFrame containing broadband magnitudes and color terms
    """

    output_df = pd.DataFrame()

    output_df['V'] = input_df['V']
    output_df['B-V'] = add_series([input_df['B'], input_df['V']], sub=True, err=err)
    output_df['U-B'] = add_series([input_df['U'], input_df['B']], sub=True, err=err)
    output_df['V-R'] = add_series([input_df['V'], input_df['R']], sub=True, err=err)
    output_df['R-I'] = add_series([input_df['R'], input_df['I']], sub=True, err=err)
    output_df['V-I'] = add_series([input_df['V'], input_df['I']], sub=True, err=err)

    output_df = output_df.apply(pd.to_numeric, errors='coerce').round(int(precision))

    return output_df


def colormag_to_ubvriframe(input_df, err=False):
    """
    Creates a Pandas DataFrame with magnitudes and color terms from an input DataFrame with V-band magnitude and
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
    output_df = output_df.replace(np.nan, 'INDEF')

    return output_df


def calculate_stdmag(obs_df, alpha, beta, err=False):
    """
    Calculates standard color-magnitude values for the secondary standards from the observed data and 
    the alpha values given the DataFrame containing beta values.
    Args:
        obs_df      : Pandas DataFrame containing observed color-magnitude data
        alpha       : List containing standard extinction coefficients
        beta        : Pandas Series containing beta values for different color terms
        err         : True, if the DataFrames in the input contain error data
    Returns:
        std_df      : Pandas DataFrame containing standard color-magnitude data
    """

    def comp_series(obs_series, alpha_val, beta_val, err=err):
        output_series = mul_series(obs_series, alpha_val)
        # obs_series = obs_series.replace('INDEF', np.nan)
        # output_series = obs_series * alpha_val
        if not err:
            return output_series.apply(lambda x: x + beta_val if x != 'INDEF' else 'INDEF')
        else:
            return output_series.apply(lambda x: (x ** 2 + beta_val ** 2) ** 0.5 if x != 'INDEF' else 'INDEF')

    std_df = pd.DataFrame()

    std_df['B-V'] = comp_series(obs_df['B-V'], alpha[0], beta['Beta_0'], err=err)
    std_df['U-B'] = comp_series(obs_df['U-B'], alpha[1], beta['Beta_1'], err=err)
    std_df['V-R'] = comp_series(obs_df['V-R'], alpha[2], beta['Beta_2'], err=err)
    std_df['R-I'] = comp_series(obs_df['R-I'], alpha[3], beta['Beta_3'], err=err)
    std_df['V-I'] = comp_series(obs_df['V-I'], alpha[4], beta['Beta_4'], err=err)
    std_df['V'] = add_series([obs_df['V'], comp_series(std_df['B-V'], alpha[5], beta['Beta_5'], err=err)], err=err)
    std_df['V1'] = add_series([obs_df['V'], comp_series(std_df['V-R'], alpha[6], beta['Beta_6'], err=err)], err=err)

    return std_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Calculating Beta Values
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_betaframe(obs_df, std_df, alpha, err=False):
    """
    Calculates Beta values for the night from the observed data and the landolt standard data using site extinction
    coefficients.
    Args:
        obs_df      : Pandas DataFrame containing observed color-magnitude data
        std_df      : Pandas DataFrame containing standard color-magnitude data
        alpha       : List containing standard extinction coefficients
        err         : True, if the DataFrames in the input contain error data
    Returns:
        output_df   : Pandas DataFrame containing beta values for different color terms
    """
    output_df = pd.DataFrame()

    output_df['Beta_0'] = add_series([std_df['B-V'], mul_series(obs_df['B-V'], alpha[0])], sub=True, err=err)
    output_df['Beta_1'] = add_series([std_df['U-B'], mul_series(obs_df['U-B'], alpha[1])], sub=True, err=err)
    output_df['Beta_2'] = add_series([std_df['V-R'], mul_series(obs_df['V-R'], alpha[2])], sub=True, err=err)
    output_df['Beta_3'] = add_series([std_df['R-I'], mul_series(obs_df['R-I'], alpha[3])], sub=True, err=err)
    output_df['Beta_4'] = add_series([std_df['V-I'], mul_series(obs_df['V-I'], alpha[4])], sub=True, err=err)
    output_df['Beta_5'] = add_series([std_df['V'], obs_df['V'], mul_series(std_df['B-V'], alpha[5])], sub=True, err=err)
    output_df['Beta_6'] = add_series([std_df['V'], obs_df['V'], mul_series(std_df['V-R'], alpha[6])], sub=True, err=err)

    return output_df


def calculate_netbeta(file_standard, std_matrix, list_alpha):
    """
    Creates a pandas DataFrame & file with magnitudes and color terms from an input file with TX'Dumped magnitudes
    from Mag files.
    Args:
        file_standard   : Text file from which the broadband magnitudes have to be extracted
        std_matrix      : Standard matrix of the landolt field whose observed intrumental magnitudes are used
        list_alpha      : List containing standard extinction coefficients
    Returns:
        netbeta         : Pandas DataFrame containing Beta values
        neterr          : Pandas DataFrame containing error of Beta values

    """
    obsmag_df, obserr_df = calculate_colormag(file_mag=file_standard)
    stdmag_df, stderr_df = matrix_to_colormagframe(std_matrix)

    betamag_df = calculate_betaframe(obsmag_df, stdmag_df, list_alpha, err=False)
    betaerr_df = calculate_betaframe(obserr_df, stderr_df, list_alpha, err=True)

    betamag_df = betamag_df.apply(pd.to_numeric, errors='coerce').round(int(precision))
    betaerr_df = betaerr_df.apply(pd.to_numeric, errors='coerce').round(int(precision))

    return betamag_df, betaerr_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Obtain Instrumental Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_colormag(file_mag):
    """
    Calculates instrumental magnitudes from Tx'dumped mag files obtained after photometry.
    Arg
        file_mag    : File containing Tx'dumped magnitudes obtained from MAG files generated by IRAF
    Returns:
        mag_df      : Pandas DataFrame containing broadband magnitudes
        err_df      : Pandas DataFrame containing errors in magnitudes
    """
    data_df = pd.read_csv(filepath_or_buffer=file_mag, sep='\s+', names=list_magcol, index_col=0, engine='python')
    star_count = len(set(data_df.index.values))

    data_df['FILTER'] = data_df['IFILTER'].apply(lambda x: str(x)[-1])
    data_df['APCOR'] = data_df['MAG_1'] - data_df['MAG_2']
    data_df['EETA'] = data_df['FILTER'].apply(lambda x: eeta[x])
    data_df['EETAERR'] = data_df['FILTER'].apply(lambda x: eeta_err[x])

    data_grouped = data_df[['APCOR', 'FILTER']].groupby(['FILTER'])
    mean = {}
    stdev = {}

    for band in set(data_df['FILTER'].values):
        temp_list = reject(data_grouped.get_group(name=band)['APCOR'].tolist(), iterations=int(star_count / 3) + 1)
        mean[band] = np.mean(temp_list)
        stdev[band] = np.std(temp_list)

    data_df['COR_MEAN'] = data_df['FILTER'].apply(lambda x: mean[x])
    data_df['COR_STD'] = data_df['FILTER'].apply(lambda x: stdev[x])

    data_df['INSTR_MAG'] = data_df['MAG_1'] - data_df['COR_MEAN'] - data_df['AIRMASS'] * data_df['EETA']
    data_df['INSTR_ERR'] = add_series([data_df['COR_STD'], data_df['ERR_1']], err=True)

    data_df = data_df[['FILTER', 'INSTR_MAG', 'INSTR_ERR']].round(int(precision))
    data_df = append_missing_data(data_df)

    mag_df = ubvrimag_to_colormagframe(unorgmag_to_ubvriframe(data_df[['FILTER', 'INSTR_MAG']]), err=False)
    err_df = ubvrimag_to_colormagframe(unorgmag_to_ubvriframe(data_df[['FILTER', 'INSTR_ERR']]), err=True)

    return mag_df, err_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# MAG Files To Be Used In Determining True Secondary Standard Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
std_matrix = PG0918
ctext_magstd = 'ca_' + day_std + '_*PG*.mag.2'
file_standard = 'output_PG' + day_std + '_mag2'

ctext_SN = 'ca_' + day_std + '*' + object_name + '*.mag.4'
file_SNfield = 'output_' + day_std + '_mag4'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Obtain DataFrame Containing Beta Values For The Nights On Which Standards Are Observed
# ------------------------------------------------------------------------------------------------------------------- #
txdump(common_text=ctext_magstd, output_file=file_standard)
if file_is_empty(file_standard):
    print ("Magnitude Files For Standard Are Not Present In The Current Directory")
    sys.exit(1)

betamag_df, betaerr_df = calculate_netbeta(file_standard, std_matrix, list_alpha)
display_text("Beta Values Were Computed Using The Landolt Standard Field")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Writes Beta Values For The Nights On Which Standards Are Observed
# ------------------------------------------------------------------------------------------------------------------- #
betamag_df.to_csv('OUTPUT_betamag', sep=' ', index=True)
betaerr_df.to_csv('OUTPUT_betaerr', sep=' ', index=True)
display_text("Beta Values Were Written Onto Files 'OUTPUT_betamag' and 'OUTPUT_betaerr'")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates Mean And STDEV of Beta For V-Band Magnitude And Color Terms
# ------------------------------------------------------------------------------------------------------------------- #
net_beta = pd.Series()
net_err = pd.Series()

for column in betamag_df:
    tempmag = reject_series(betamag_df[column])
    temperr = betaerr_df[column].loc[tempmag.index.values]

    net_beta[column] = tempmag.mean()
    net_err[column] = ((tempmag.std() ** 2) + temperr.apply(lambda x: (x / len(temperr.index)) ** 2).sum()) ** 0.5

net_beta = net_beta.round(int(precision))
net_err = net_err.round(int(precision))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate True Magnitudes Of Secondary Standards By Calculating Instrumental Magnitudes & Beta Values
# ------------------------------------------------------------------------------------------------------------------- #
txdump(common_text=ctext_SN, output_file=file_SNfield)
if file_is_empty(file_SNfield):
    print ("SN Magnitude Files For {0} Are Not Present In The Current Directory".format(day_std))
    sys.exit(1)

instrmag_df, instrerr_df = calculate_colormag(file_SNfield)

stdmag_df = calculate_stdmag(instrmag_df, list_alpha, net_beta, err=False)
stderr_df = calculate_stdmag(instrerr_df, list_alpha, net_err, err=True)

truemag_df = colormag_to_ubvriframe(stdmag_df, err=False)
trueerr_df = colormag_to_ubvriframe(stderr_df, err=True)

truemag_df.to_csv('OUTPUT_truestdmag', sep=' ', index=True)
trueerr_df.to_csv('OUTPUT_truestderr', sep=' ', index=True)

display_text("True Secondary Standard Magnitudes Were Computed")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Output The Photometric Magnitudes Onto A Latex Table
# ------------------------------------------------------------------------------------------------------------------- #
comb_df = pd.DataFrame(index=truemag_df.index.values)
comb_df.index.name = 'ID'

for band in filters:
    comb_df[band] = truemag_df[band].apply(lambda x: "{:.2f}".format(x)) + r"$\pm$" + \
                    trueerr_df[band].apply(lambda x: "{:.2f}".format(x))

comb_df.reset_index().to_latex('_SecStdMag.tex', escape=False, index=False)
display_text("True Secondary Standard Magnitudes Have Been Written Onto A Latex File")
# ------------------------------------------------------------------------------------------------------------------- #
