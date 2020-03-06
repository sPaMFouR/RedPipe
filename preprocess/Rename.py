#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxx-------------------------RENAME RAW FILES-----------------------xxxxxxxxxxxxxxxxxxxxxxxx #
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
read_noise = 5.75
ccd_gain = 0.28
data_max = 700000
fields_extr = '$I, OBJECT, EXPTIME, FILTER, GRISM, COMMENT'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
iraf.imred(_doprint=0)
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

    task.verify = 'no'              # Verify Each Edit Operation?
    task.add = add_keyword          # Add Rather Than Edit Fields?
    task.show = 'no'                # Print Record Of Each Edit Operation?
    task.update = 'yes'             # Enable Updating Of The Image Header?

    task(images='@' + textlist_files, fields=field_name, value=value)


def hselect(ctext, fields=fields_extr):
    """
    Selects the header key specified by the 'fields' of all the FITS files having the common_text 'ctext'.
    Args:
        ctext       : Common text of the names of files whose header is to be extracted
        field_name  : The header keyword to be edited for the above files
    Returns:
        None
    """
    task = iraf.images.imutil.hselect
    task.unlearn()

    task.missing = 'INDEF'          # Value For Missing Keywords

    task(images=ctext, fields=fields, expr='yes', Stdout='list_summary.cl')


def rename(textlist_files, value='fits', field_name='extn'):
    """
    Add '.fits' extension to the Raw FITS files in the list 'textlist_files'.
    Args:
        textlist_files  : Text list containing the names of files whose header is to be edited
        value           : New File name or Field name
        field_name      : Field to be modified (all|dir|root|extn)
    Returns:
        None
    """
    task = iraf.system.rename
    task.unlearn()

    task(files='@' + textlist_files, newname=value, field=field_name)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# remove_resfile = eg.boolbox(msg='Remove Residual Files From Previous Run Of This Script?',
#                             title='Remove Residual Files', choices=['Yes', 'No'])
remove_resfile = True
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of Rename.py
# ------------------------------------------------------------------------------------------------------------------- #
if remove_resfile:
    for text in ['list_*', '*.cl']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Groups Similar Type Of FITS Files
# Edits The 'OBJECT' Keyword In The Headers Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
textlist_files = 'list_files'
group_similar_files(textlist_files, common_text='*', exceptions='.fits, *.py')

rename(textlist_files)
hselect(ctext='*.fits[0]', fields=fields_extr)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Removes All The Text LISTS Created In The Current Directory
# ------------------------------------------------------------------------------------------------------------------- #
remove_similar_files('*list*')
# ------------------------------------------------------------------------------------------------------------------- #
