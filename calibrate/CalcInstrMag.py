#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxx-----------Calculate Instrumental Magnitudes From TXDumped Photometry Files--------------xxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
filters = ['U', 'B', 'V', 'R', 'I']
list_magcol = ['ID', 'IMAGE', 'IFILTER', 'XCENTER', 'YCENTER', 'SKY_COUNTS', 'AIRMASS',
               'APER_1', 'APER_2', 'MAG_1', 'MAG_2', 'ERR_1', 'ERR_2']
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
    list_reject = filter(lambda x: x != np.nan, list_values)
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
# Function To Obtain Instrumental Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

def add_series(list_series, sub=False, err=False):
    """
    Adds multiple Pandas Series Column wise and obtains
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
            data_ext = [[band] + ['INDEF'] * (len(input_df.columns.values) - 1) for _ in range(0, len(star_id))]
            input_df = pd.concat([pd.DataFrame(data_ext, columns=input_df.columns.values, index=star_id), input_df])

    output_df = input_df.sort_values(by='FILTER').sort_index(kind='mergesort')
    output_df = output_df.replace('INDEF', np.nan, regex=True)

    return output_df


def unorgmag_to_ubvriframe(input_df):
    """
    Creates a Pandas DataFrame with broadband magnitudes from an input DataFrame with unorganised magnitudes.
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
    output_df = output_df.apply(pd.to_numeric, errors='coerce').round(decimals=int(precision))
    output_df = output_df.replace(np.nan, 'INDEF', regex=True)

    return output_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Obtain Instrumental Magnitudes From MAG Files
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_instrmag(file_mag):
    """
    Calculates instrumental magnitudes from Tx'dumped mag files obtained after photometry.
    Arg
        file_mag    : File containing Tx'dumped magnitudes obtained from MAG files generated by IRAF
    Returns:
        mag_df      : Pandas DataFrame containing broadband magnitudes
        err_df      : Pandas DataFrame containing errors in magnitudes
    """
    file_df = pd.read_csv(filepath_or_buffer=file_mag, sep="\s+", names=list_magcol, index_col=0, engine='python')
    file_df = file_df.replace('INDEF', np.nan)
    file_df[['MAG_1', 'MAG_2']] = file_df[['MAG_1', 'MAG_2']].astype('float64')

    indexes = file_df.index.values
    star_count = len(set(indexes))

    file_df['FILTER'] = file_df['IFILTER'].apply(lambda x: str(x)[-1])
    file_df['APCOR'] = file_df['MAG_1'] - file_df['MAG_2']

    data_grouped = file_df[['APCOR', 'FILTER']].groupby(['FILTER'])
    mean = {}
    stdev = {}

    for band in set(file_df['FILTER'].values):
        temp_list = reject(data_grouped.get_group(name=band)['APCOR'].tolist(), iterations=int(star_count / 3) + 1)
        mean[band] = np.mean(temp_list)
        stdev[band] = np.std(temp_list)

    file_df['COR_MEAN'] = file_df['FILTER'].apply(lambda x: mean[x])
    file_df['COR_STD'] = file_df['FILTER'].apply(lambda x: stdev[x])
    file_df['INSTR_MAG'] = file_df['MAG_1'] - file_df['COR_MEAN']
    file_df['INSTR_ERR'] = add_series([file_df['COR_STD'], file_df['ERR_1']], err=True)

    file_df = file_df.round(int(precision))
    file_df = file_df[['FILTER', 'INSTR_MAG', 'INSTR_ERR']]
    file_df = append_missing_data(file_df)

    mag_df = unorgmag_to_ubvriframe(file_df[['FILTER', 'INSTR_MAG']])
    err_df = unorgmag_to_ubvriframe(file_df[['FILTER', 'INSTR_ERR']])

    date = file_mag.split("_")[1]
    mag_df.to_csv("OUTPUT_instrmag_" + date, sep=" ", index=True)
    err_df.to_csv("OUTPUT_instrerr_" + date, sep=" ", index=True)

    net_df = pd.DataFrame()
    for column in mag_df:
        net_df[column] = mag_df[column].apply(lambda x: str(x) + "+/-") + err_df[column].apply(lambda x: str(x))

    net_df.to_csv("OUTPUT_instr_" + date, sep=" ", index=True)

    return mag_df, err_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Instrumental Magnitudes For Each Day From The Mag Files
# ------------------------------------------------------------------------------------------------------------------- #
list_mag = group_similar_files("", "*_mag4")

for file_name in list_mag:
    calculate_instrmag(file_name)

display_text("Instrumental Magnitudes Have Been Computed From MAG Files")
# ------------------------------------------------------------------------------------------------------------------- #
