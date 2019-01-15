#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxx-------------------PHOT0METRIC DATA REDUCTION PIPELINE-----------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import shutil
import easygui
import numpy as np
import matplotlib.pyplot as plt
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Location Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DEFAULT_DIR = "/home/avinash/Supernovae_Data/2016gfy/"
# DIR_SNE = easygui.enterbox(msg='Enter the Location Of The Directory:', title='Directory For Reduction',
#                                 default=DEFAULT_DIR)
DIR_SNE = "/home/avinash/Supernovae_Data/2016gfy/"
DIR_PHOT = "/home/avinash/Supernovae_Data/2016gfy/Photometry/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions In File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def remove_file(file_name):
    """
    Removes the file 'file_name' in the constituent directory.
    Args:
         file_name  : Name of the file to be removed
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
        common_text : String containing partial filename
    Returns:
        None
    """
    for residual_files in glob.glob(common_text):
        os.remove(residual_files)


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


def text_list_to_python_list(text_list):
    """
    Returns data in the file 'text_list' as a python_list.
    Args:
        text_list   : Input file containing filenames
    Returns:
        python_list : List of all the elements in the file 'text_list'
    Raises:
        Error : File Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("Error : File '" + text_list + "' Not Found")


def python_list_to_text_list(python_list, text_list):
    """
    Put the data from the input 'python_list' to a file 'text_list' line-wise.
    Args:
        python_list : Name of the python_list from which data has to be read
        text_list   : Name of the text file onto which data has to be appended
    Returns:
        None
    """
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(element + '\n')


def list_lists_to_list(list_lists, text_list):
    """
    Groups filenames from a list 'list_lists' onto a single file 'text_list'.
    Args:
        list_lists  : List containing the names of different lists
        text_list   : Name of the file onto which all the filenames from the 'list_lists' have to be appended
    Returns:
        list_files  : Python list containing the names of all the constituent files
    """
    list_name = []
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list = f.read().split()
            for element in file_list:
                list_name.append(element)

    python_list_to_text_list(list_files, text_list)

    return list_files


def append_date(file_name, date):
    """
    Appends date of observation mentioned in string 'date' to the file 'file_name' and creates a copy
    of the same file.
    Args:
        file_name       : File for which the date appended copy has to be created
        date            : Date of the observation of the file
    Returns:
        output_filename : Final name of the
    """
    output_filename = date + '_' + file_name
    remove_file(output_filename)
    shutil.copy(file_name, output_filename)

    return output_filename


def execute_script(exec_dir, python_script):
    """
    Executes the python script in the mentioned directory 'home_dir'. Removes the script
    from the 'home_dir' after execution.
    Args:
        exec_dir      : Directory where the python script is to be run
        python_script : Name of the python script that is to be executed
    Returns:
        None
    """
    os.system('cp ' + python_script + ' ' + exec_dir + ' 1')
    os.system('python ' + exec_dir + python_script + ' 1')
    os.system('rm ' + exec_dir + python_script + ' 1')

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Collect All Folders In The Supernova Directory
# ------------------------------------------------------------------------------------------------------------------- #
list_folders = glob.glob(DIR_SNE + '*/')
list_folders.sort()
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Pre-Processing Of Observational Data
# # ------------------------------------------------------------------------------------------------------------------- #
# for index in range(0, len(list_folders)):
#     exec_dir = DIR_SNE + list_folders[index]
#     os.chdir(exec_dir)
#     execute_script(exec_dir, 'pre_processing.py')
# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Prefix Photometric Object Files In A Folder With Date & Copy Photometric Files Into A New Directory 'Photometry'
# ------------------------------------------------------------------------------------------------------------------- #
# os.mkdir(str(DIR_PHOT))
for folder in list_folders:
    os.chdir(folder)
    date = folder[:-1].split("/")[-1]
    list_object = group_similar_files('', 'fbs_2016gfy*.fits')
    for file_name in list_object:
        shutil.copy(file_name, date + '_' + file_name)
        shutil.copy(date + '_' + file_name, DIR_PHOT)

os.chdir(DIR_SNE)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Collecting Final Spectras
# ------------------------------------------------------------------------------------------------------------------- #
#os.chdir(DIR_SPEC)
#list_raw_folders = glob.glob('*/')
#list_folders = [date[:-1] for date in list_raw_folders if len(date[:-1]) == 5]
#list_folders.sort()

#try:
#    os.mkdir(DEFAULT_DIR + 'Final_Spectra/')
#except OSError:
#    pass

#for date in list_folders:
#    exec_dir = DIR_SPEC + date
#    os.chdir(exec_dir)
#    list_cfwcbs = group_similar_files('', 'cfwcbs*.fits')
#    for file_name in list_cfwcbs:
#        output_filename = append_date(file_name, date)
#        shutil.copy(str(output_filename), DEFAULT_DIR + 'Final_Spectra/')
# ------------------------------------------------------------------------------------------------------------------- #

# # ------------------------------------------------------------------------------------------------------------------- #
# # Alignment Of Images
# # ------------------------------------------------------------------------------------------------------------------- #
# execute_script(DIR_PHOT, 'align_dft.py')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Photometry Of Object Frames
# # ------------------------------------------------------------------------------------------------------------------- #
# execute_script(DIR_PHOT, 'photometry.py')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Spectroscopy Of Object Frames
# # ------------------------------------------------------------------------------------------------------------------- #
# execute_script(DIR_SPEC, 'spectroscopy.py')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Flux Calibration Of Final Spectras
# # ------------------------------------------------------------------------------------------------------------------- #
# execute_script(DIR_SPEC, 'fluxcalib.py')
# # ------------------------------------------------------------------------------------------------------------------- #

