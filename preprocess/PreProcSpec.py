#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxx-------------------------SPECTROSCOPIC DATA REDUCTION-----------------------xxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import numpy as np
from pyraf import iraf
from astropy.io import fits
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Telescope CCD Specifications
# ------------------------------------------------------------------------------------------------------------------- #
read_noise = 4.87
ccd_gain = 1.22
data_max = 55000
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
FILTER_keyword = 'FILTER'
GRISM_keyword = 'GRISM'
OBJECT_keyword = 'OBJECT'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_NAME = '2018cow'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Clipping Section (In Pixels) For The Spectroscopic Frames
# ------------------------------------------------------------------------------------------------------------------- #
clipsec_twod = '[*, 51:3200]'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.images(_doprint=0)
iraf.crutil(_doprint=0)
iraf.ccdred.instrument = 'ccddb$kpno/camera.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
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
    for residual_files in glob.glob(common_text):
        os.remove(residual_files)


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


def text_list_to_python_list(text_list):
    """
    Returns data in the file 'text_list' as a python_list.
    Args:
        text_list   : Input file containing filenames
    Returns:
        python_list : List of all the elements in the file 'text_list'
    Raises:
        Error       : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("Error : File {0} Not Found".format(text_list))


def python_list_to_text_list(python_list, text_list):
    """
    Put the data from the input 'python_list' to a file 'text_list' line-wise.
    Args:
        python_list : Python_list from which data has to be read
        text_list   : Name of the text file onto which data has to be appended
    Returns:
        None
    """
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element) + '\n')


def list_lists_to_list(list_lists, text_list):
    """
    Groups filenames from a list 'list_lists' onto a single file 'text_list'.
    Args:
        list_lists      : List containing the names of different lists
        text_list       : Name of the file onto which all the filenames from the 'list_lists' have to be appended
    Returns:
        list_filename   : Python list containing the names of all the constituent files
    """
    list_filename = []
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            list_contents = f.read().split()
            for element in list_contents:
                list_filename.append(element)
    python_list_to_text_list(list_filename, text_list)

    return list_filename


def check_fileempty(file_name):
    """
    Args:
        file_name   : Name of the file to be checked whether it's empty
    Returns:
        bool        : Returns true if file is empty
    """
    return not bool(len(open(file_name).readlines()))


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

def hedit(textlist_files, field_name, value, add_keyword='no'):
    """
    Edits the header key specified by the 'field_name' of all the FITS files in the file 'textlist_files'
    and substitutes it with 'value'.
    Args:
        textlist_files  : Text list containing the names of files whose header is to be edited
        add_keyword     : Whether header keyword is to be added to the files
        field_name      : The header keyword to be edited for the above files
        value           : The new value of the header keyword
    Returns:
        None
    """
    task = iraf.images.imutil.hedit
    task.unlearn()

    task.verify = 'no'  # Verify Each Edit Operation?
    task.add = add_keyword  # Add Rather Than Edit Fields?
    task.show = 'no'  # Print Record Of Each Edit Operation?
    task.update = 'yes'  # Enable Updating Of The Image Header?

    task(images='@' + textlist_files, fields=field_name, value=value)


def zero_combine(textlist_bias, master_bias):
    """
    Combines Bias images in the file 'textlist_bias' to form a master_bias image.
    Args:
        textlist_bias   : Text list containing the names of photometric bias images
        master_bias     : Name of the master bias file
    Returns:
        None
    """
    if check_fileempty(textlist_bias):
        display_text("Error: No BIAS Images To Combine In Directory {0}".format(os.getcwd()))
        sys.exit(1)

    task = iraf.noao.imred.ccdred.zerocombine
    task.unlearn()

    task.combine = 'median'  # Type Of Combine Operation
    task.reject = 'avsigclip'  # Type Of Rejection
    task.ccdtype = ''  # CCD Image Type To Combine
    task.process = 'no'  # Process Images Before Combining?
    task.delete = 'no'  # Delete Input Images After Combining?
    task.rdnoise = float(read_noise)  # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)  # CCD Gain (In e-/ADU)

    remove_file(master_bias)
    task(input='@' + textlist_bias, output=master_bias)


def bias_subtract(textlist_tbs, master_bias, prefix_str='bs_'):
    """
    Subtracts the master_bias image from the files in the file 'textlist_tbs'.
    Args:
        textlist_tbs    : Text list containing the filenames from which bias is to be subtracted
        master_bias     : Name of the master bias file
        prefix_str      : Prefix to distinguish the output FITS file from the original FITS file
    Returns:
        None
    """
    list_files = text_list_to_python_list(textlist_tbs)
    task = iraf.images.imutil.imarith
    task.unlearn()

    data = fits.getdata(filename=master_bias, ext=0)
    median = np.median(data[30:220, 100:3100])

    for file_name in list_files:

        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(operand1=file_name, op='-', operand2=median, result=output_filename)


def imcopy(textlist_clip, clip_section, prefix_str='c'):
    """
    Clips the file 'file_name' based on the string 'clip_section'.
    Args:
        textlist_clip   : Text list containing the names of the FITS files to be clipped
        clip_section    : Python list containing the section of the FITS file to be copied
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_clip = text_list_to_python_list(textlist_clip)

    task = iraf.images.imutil.imcopy
    task.unlearn()

    for file_name in list_clip:
        output_filename = prefix_str + file_name
        if prefix_str != '':
            remove_file(output_filename)
        task(input=file_name + clip_section, output=output_filename, verbose='no')


def crmedian(textlist_cosmic, clip_section, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        textlist_cosmic : Text list containing the names of FITS files to be corrected for cosmic rays
        clip_section    : Python list containing the section of the FITS file to be copied
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(textlist_cosmic)
    imcopy(textlist_clip=textlist_cosmic, clip_section=clip_section, prefix_str='')

    task = iraf.noao.imred.crutil.crmedian
    task.unlearn()

    task.lsigma = 25                                        # Low Clipping Sigma Factor
    task.ncsig = 10                                         # Column Box Size For Sigma Calculation

    for file_name in list_cosmic:
        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(input=file_name, output=output_filename, crmask='', median='', sigma='', residual='')


def cosmicray_check(textlist_cosmic, prefix_str='cr_'):
    """
    Subtracts the cosmic ray corrected image from the original image. Both the images are clipped.
    This is performed only after cosmic ray correction has been performed on images.
    Args:
        textlist_cosmic : Text list containing the names of FITS files to be corrected for cosmic rays
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(textlist_cosmic)

    task = iraf.images.imutil.imarith
    task.unlearn()

    for file_name in list_cosmic:
        output_filename = prefix_str + file_name[3:]
        remove_file(output_filename)
        task(operand1=file_name, op='-', operand2='c' + file_name, result=output_filename)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# remove_resfile = eg.boolbox(msg='Remove Residual Files From Previous Run Of This Script?',
#                             title='Remove Residual Files', choices=['Yes', 'No'])
remove_resfile = True
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of PreProcSpec.py
# ------------------------------------------------------------------------------------------------------------------- #
if remove_resfile:
    for text in ['*.ms.*', 'bs_*', '*cbs_*', 'list_*', '*tmp&']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Similar Type Of FITS Files
# Edits The 'OBJECT' Keyword In The Headers Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
dict_imgtype = {'biasspec*.fits': ['list_biasspec', 'BIAS_SPEC'], 'FeAr*.fits': ['list_FeAr', 'FeAr'],
                'FeNe*.fits': ['list_FeNe', 'FeNe'], OBJECT_NAME + '*.fits': ['list_object', OBJECT_NAME],
                'Feige34*.fits': ['list_Feige34', 'Feige34'], 'Feige66*.fits': ['list_Feige66', 'Feige66'],
                'Feige110*.fits': ['list_Feige110', 'Feige110'], 'Feige*.fits': ['list_Feige', ''],
                '*gr7*.fits': ['list_gr7', ''], '*gr8*.fits': ['list_gr8', '']}

for (text, [listname, objname]) in dict_imgtype.items():
    group_similar_files(listname, common_text=text)
    if objname != '':
        hedit(listname, OBJECT_keyword, objname)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Edits The 'FILTER' Keyword In The Headers Of FITS Files
# Edits The 'GRISM' Keyword In The Headers Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
hedit('list_spec', FILTER_keyword, '1Free')

hedit('list_gr7', GRISM_keyword, '4Grism7')
hedit('list_gr8', GRISM_keyword, '3Grism8')
hedit('list_FeAr', GRISM_keyword, '4Grism7')
hedit('list_FeNe', GRISM_keyword, '3Grism8')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Combines The BIAS Frames To Obtain MASTER_BIAS
# Applies BIAS Correction To The Files In The list 'list_tbs'
# ------------------------------------------------------------------------------------------------------------------- #
list_list_spec = ['list_biasspec', 'list_object', 'list_FeNe', 'list_FeAr', 'list_Feige']
list_lists_to_list(list_list_spec, 'list_spec')
list_lists_to_list(list_list_spec[1:], 'list_tbs_spec')

if not check_fileempty('list_biasspec'):
    zero_combine(textlist_bias='list_biasspec', master_bias='mbiasspec.fits')
    bias_subtract(textlist_tbs='list_tbs_spec', master_bias='mbiasspec.fits')
else:
    display_text("Error: No BIAS Images To Combine In Directory {0}".format(os.getcwd()))
    sys.exit(1)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups BIAS Subtracted Spectra Of Objects+Standards & Lamps
# Clips The BIAS Subtracted Spectroscopic Images And Corrects For Cosmic Ray Hits
# Computes The Subtracted Images Through Subtraction
# ------------------------------------------------------------------------------------------------------------------- #
group_similar_files('list_bs_gr', common_text='bs_*gr*.fits')
group_similar_files('list_bs_lamp', common_text='bs_Fe*.fits', exceptions='gr')

crmedian(textlist_cosmic='list_bs_gr', clip_section=clipsec_twod)
imcopy(textlist_clip='list_bs_lamp', clip_section=clipsec_twod)
cosmicray_check('list_bs_gr')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes All The Text LISTS Created In The Current Directory
# ------------------------------------------------------------------------------------------------------------------- #
remove_similar_files('*list*')
# ------------------------------------------------------------------------------------------------------------------- #
