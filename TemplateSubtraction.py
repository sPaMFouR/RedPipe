#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxx--------------------------TEMPLATE SUBTRACTION------------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

"""
Assumes that you have aligned template photometric images and just need to carry out template subtraction on the
photometric images.
"""

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
# Details of Template Images
# ------------------------------------------------------------------------------------------------------------------- #
dict_temp = {'u': ['template_U.fits', 8.30], 'b': ['template_B.fits', 6.90], 'v': ['template_V.fits', 6.10],
             'r': ['template_R.fits', 5.90], 'i': ['template_I.fits', 5.30]}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
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


def text_list_to_python_list(text_list):
    """
    Returns data in the file 'text_list' as a python_list.
    Args:
        text_list   : Input file containing filenames
    Returns:
        python_list : List of all the elements in the file 'text_list'
    Raises:
        Error : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, "r+") as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("Error : File " + str(text_list) + " Not Found")
        sys.exit(1)


def list_statistics(python_list):
    """
    Returns the statistics of the list of elements in the input 'python_list'.
    Args:
        python_list  : Input list of elements
    Returns:
        value_mean  : Mean of the list of elements
        value_median: Median of the list of elements
        value_std   : Standard Deviation of the list of elements
    """
    value_mean = np.mean(python_list)
    value_median = np.median(python_list)
    value_std = np.std(python_list)

    return value_mean, value_median, value_std


def reject(python_list, iterations=2):
    """
    Rejects outliers from the input 'python_list'.
    Args:
        python_list : Input list of elements
        iterations  : No. of iterations of rejection to be run on the input list
    Returns:
        reject_list : Output list after rejecting outliers from the input 'python_list'
    """
    reject_list = [float(val) for val in python_list]
    reject_list.sort()

    for _ in range(0, iterations):
        if len(python_list) > 2:
            value_mean, value_median, value_std = list_statistics(reject_list)

            if abs(reject_list[0] - value_median) < abs(reject_list[-1] - value_median):
                remove_index = -1
            else:
                remove_index = 0

            if abs(reject_list[remove_index] - value_median) > value_std:
                reject_list.pop(remove_index)

    return reject_list


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
# IRAF Tasks
# ------------------------------------------------------------------------------------------------------------------- #

def imexam(file_name, def_key='a', coord_file='stars.coo', log_imexam='log_imexam'):
    """
    Examines the image 'file_name' at the coordinates mentioned in the file 'coord_file'
    and logs the output onto the file 'log_imexam'.
    Args:
        file_name    : Name of the FITS file
        def_key      : Default key for cursor input
        coord_file   : Text file listing the coordinates of selected stars in the field
        log_imexam   : Name of the text list to record log of IMEXAM
    Returns:
        None
    """
    remove_file(str(log_imexam))

    task = iraf.images.tv.imexam
    task.unlearn()

    task.logfile = str(log_imexam)                  # Log File To Record Output Of The Commands
    task.keeplog = 'yes'                            # Log Output Results?
    task.defkey = str(def_key)                      # Default Key For Cursor x-y Input List
    task.imagecur = str(coord_file)                 # Image Display Cursor Input
    task.use_display = 'no'                         # Use The Image Display?

    task(input=str(file_name), frame=1)


def imexam_multi(text_list, def_key='a', coord_file='stars.coo', log_imexam='log_imexam'):
    """
    Examines the images in the list 'python_list' at the coordinates mentioned in the file "stars.coo"
    and logs the output onto the file "log_imexam".
    Args:
        text_list    : Text list containing names of FITS files
        def_key      : Default key for cursor input
        coord_file   : Text file listing the coordinates of selected stars in the field
        log_imexam   : Name of the text list to record log of IMEXAM
    Returns:
        None
    """
    remove_file(str(log_imexam))
    list_files = text_list_to_python_list(str(text_list))

    task = iraf.images.tv.imexam
    task.unlearn()

    for file_name in list_files:
        task.logfile = str(log_imexam)              # Log File To Record Output Of The Commands
        task.keeplog = 'yes'                        # Log Output Results?
        task.defkey = str(def_key)                  # Default Key For Cursor x-y Input List
        task.imagecur = str(coord_file)             # Image Display Cursor Input
        task.use_display = 'no'                     # Use The Image Display?

        task(input=str(file_name), frame=1)


def imarith(operand1, operand2, operation='-', prefix_str='skysub_'):
    """
    Subtracts the master_bias image from the files in the file 'text_list_tbs'.
    Args:
        operand1      : Value of the 1st operand
        operand2      : Value of the 2nd operand
        operation     : Operation to be performed on the files
        prefix_str    : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    task = iraf.images.imutil.imarith
    task.unlearn()

    output_file_name = str(prefix_str) + str(operand1)
    remove_file(str(output_file_name))
    task(operand1=operand1, op=str(operation), operand2=operand2, result=str(output_file_name))


def gauss(file_name, sigma, prefix_str='conv_'):
    """
    Convolve the input image with an elliptical gaussian function.
    Args:
        file_name       : Input image to be fit
        sigma           : Sigma of gaussian along major axis of ellipse
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        output_filename : Name of the convolved FITS file
    """
    task = iraf.images.imfilter.gauss
    task.unlearn()

    if file_name.startswith('skysub_'):
        output_file_name = str(prefix_str) + file_name
    else:
        output_file_name = str(prefix_str) + file_name.lstrip('s_')

    remove_file(str(output_file_name))
    task(input=str(file_name), output=str(output_file_name), sigma=sigma)

    return output_file_name


def calculate_fwhm(text_list_files, coord_file='stars.coo', log_imexam='log_fwhm'):
    """
    Calculates the Mean FWHM of all the files in the list 'text_list_files'. It determines the FWHM
    using the IMEXAMINE task on the stars mentioned in the file "stars.coo".
    Args:
        text_list_files : Text list containing names of FITS files whose FWHM is to be determined
        coord_file      : Text file listing the coordinates of selected stars in the field
        log_imexam      : Name of the text list to record log of IMEXAM
    Returns:
        list_mean_fwhm  : Python list containing Mean FWHM of all the FITS files
    """
    imexam_multi(str(text_list_files), def_key='a', coord_file=str(coord_file), log_imexam=str(log_imexam))
    coord_df = pd.read_csv(coord_file, sep="\s+", header=None, engine='python')
    data_df = pd.read_csv(log_imexam, sep="\s+", comment="#", header=None, engine='python')
    count = coord_df.shape[0]
    rows, columns = data_df.shape
    col_moffat = columns - 2

    list_fwhm = [data_df.iloc[0 + count * i: (i + 1) * count, col_moffat].tolist() for i in range(0, rows / count)]

    list_mean_fwhm = []
    for fwhm_values in list_fwhm:
        mean = float(np.mean(a=reject(fwhm_values)))
        list_mean_fwhm.append(round(mean, 1))

    display_text("FWHM Of All The FITS Files Have Been Computed")

    return list_mean_fwhm


def calculate_skybg(text_list_files, coord_file='sky.coo', log_imexam='log_galaxy'):
    """
    Calculates the Mean sky background of all the files in the list 'text_list_files'. It determines the
    sky background using the IMEXAMINE task on the coordinates mentioned in the file "sky.coo".
    Args:
        text_list_files : Text list containing names of FITS files whose sky background is to be determined
        coord_file      : Text file listing the coordinates of open sky in the field
        log_imexam      : Name of the text list to record log of IMEXAM
    Returns:
        list_mean_skybg : Python list containing Mean sky background of all the FITS files
    """
    imexam_multi(str(text_list_files), def_key='m', coord_file=str(coord_file), log_imexam=str(log_imexam))
    coord_df = pd.read_csv(coord_file, sep="\s+", header=None, engine='python')
    data_df = pd.read_csv(log_imexam, sep="\s+", comment="#", header=None, engine='python')
    count = coord_df.shape[0]
    rows, columns = data_df.shape
    col_mean = 2

    list_skybg = [data_df.iloc[0 + count * i: (i + 1) * count, col_mean].tolist() for i in range(0, rows / count)]

    list_mean_skybg = []
    for skybg_values in list_skybg:
        mean = float(np.mean(a=reject(skybg_values)))
        list_mean_skybg.append(round(mean, 1))

    return list_mean_skybg


def calculate_flux(file_name, coord_file='galaxy.coo', log_imexam='log_galaxy'):
    """
    Calculates the Mean sky background of all the files in the list 'text_list_files'. It determines the
    sky background using the IMEXAMINE task on the coordinates mentioned in the file "sky.coo".
    Args:
        file_name   : Name of the FITS file whose galaxy count is to be determined
        coord_file  : Text file listing the coordinates of open sky in the field
        log_imexam  : Name of the text list to record log of IMEXAM
    Returns:
        list_flux   : Python list containing Mean sky background of all the FITS files
    """
    imexam(file_name, def_key='m', coord_file=str(coord_file), log_imexam=str(log_imexam))
    data_df = pd.read_csv(log_imexam, sep="\s+", comment="#", header=None, engine='python')
    list_galaxy = data_df[2]

    display_text("Galaxy Background Has Been Computed For - " + str(file_name))

    return list_galaxy

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of Photometry Tasks (PHOT, PSTSELECT, PSF, ALLSTAR)
# ------------------------------------------------------------------------------------------------------------------- #
remove_resfile = True
if remove_resfile:
    for text in ['skysub_*', 'conv_*', 'sc*', 'ts_*', 'list_*', 'log_*']:
        remove_similar_files(common_text=str(text))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Sky Background Of Object Frames And Subtract It From The Original Frame
# ------------------------------------------------------------------------------------------------------------------- #

for band in dict_temp.keys():
    temp_list = group_similar_files('list_object' + band, common_text='*ASASSN14dq-' + band + '*.fits',
                                    exceptions='psf,sub')
    list_skybg = calculate_skybg('list_object' + band)
    for index, file_name in enumerate(temp_list):
        imarith(operand1=file_name, operand2=list_skybg[index])
        display_text("Sky Background Has Been Subtracted From - " + file_name)

display_text("Sky Background Has Been Subtracted From All The Object Frames")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate FWHM From Sky Subtracted Object Frames And Convolve Frames Based On FWHM Values
# Also Performs Template Subtraction Of The Object Frames
# ------------------------------------------------------------------------------------------------------------------- #

for band in dict_temp.keys():
    obj_list = group_similar_files('list_skysub' + band, common_text='skysub_*-' + band + '*.fits')
    list_fwhm = calculate_fwhm('list_skysub' + band)
    for index, file_name in enumerate(obj_list):
        sigma_gauss = abs(list_fwhm[index] ** 2 - dict_temp[band][1] ** 2) ** 0.5 / 2.3548
        print("Sigma_Gaussian = {0}".format(sigma_gauss))

        if list_fwhm[index] > dict_temp[band][1]:
            print("Template Frame Will Be Convolved")
            temp_file = gauss(dict_temp[band][0], sigma=sigma_gauss)
            obj_file = file_name
        elif list_fwhm[index] < dict_temp[band][1]:
            print("Object Frame Will Be Convolved")
            temp_file = dict_temp[band][0]
            obj_file = gauss(file_name, sigma=sigma_gauss)
        else:
            temp_file = dict_temp[band][0]
            obj_file = file_name

        list_tempflux = calculate_flux(temp_file)
        list_objflux = calculate_flux(obj_file)

        scale_factor = np.mean([x / y for x, y in zip(list_objflux, list_tempflux)])
        imarith(operand1=temp_file, operation='*', operand2=scale_factor, prefix_str='sc')
        imarith(operand1=file_name, operation='-', operand2='sc' + temp_file, prefix_str='ts_')

        display_text("Template Subtraction Has Been Performed On - " + str(file_name))

remove_similar_files(common_text="conv_*")
remove_similar_files(common_text="skysub_*")

# ------------------------------------------------------------------------------------------------------------------- #
