#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx-------------------------COSMIC RAY CORRECTION--------------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import cosmics as cs
import easygui as eg
from pyraf import iraf
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
grism_keyword = 'GRISM'
filter_keyword = 'IFILTER'
object_keyword = 'OBJECT'
exptime_keyword = 'EXPTIME'
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
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("Error : File '{0}' Not Found".format(text_list))


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
            f.write(str(element) + "\n")


def list_lists_to_list(list_lists, text_list):
    """
    Groups filenames from a list 'list_lists' onto a single file 'text_list'.
    Args:
        list_lists  : List containing the names of different lists
        text_list   : Name of the file onto which all the filenames from the 'list_lists' have to be appended
    Returns:
        list_name   : Python list containing the names of all the constituent files
    """
    list_name = []
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list = f.read().split()
            for element in file_list:
                list_name.append(element)
    python_list_to_text_list(list_name, text_list)

    return list_name

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.crutil(_doprint=0)
iraf.images(_doprint=0)
iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def hedit(textlist_files, field_name, value, add_keyword='no'):
    """
    Edits the header key specified by the 'field_name' of all the FITS files in the file 'text_list_files'
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

    task.verify = 'no'                                      # Verify Each Edit Operation?
    task.add = add_keyword	                            # Add Rather Than Edit Fields?
    task.show = 'no'                                        # Print Record Of Each Edit Operation?
    task.update = 'yes'                                     # Enable Updating Of The Image Header?

    task(images='@' + textlist_files, fields=field_name, value=value)


def cosmicrays(ctext, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image.
    Args:
        ctext           : Common text of all the files to be corrected for cosmic rays
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_crreject = group_similar_files("", common_text=ctext)

    task = iraf.noao.imred.crutil.cosmicrays
    task.unlearn()

    for image in list_crreject:
        output_filename = prefix_str + image
        remove_file(output_filename)
        task(input=image, output=output_filename)


def crmedian(textlist_cosmic, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        textlist_cosmic     : Text list containing the names of FITS files to be corrected for cosmic rays
        prefix_str          : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(textlist_cosmic)

    task = iraf.noao.imred.crutil.crmedian
    task.crmask = ''
    task.median = ''
    task.sigma = ''
    task.resid = ''
    task.unlearn()

    task.lsigma = 25                                            # Low Clipping Sigma Factor
    task.ncsig = 10                                             # Column Box Size For Sigma Calculation

    for file_name in list_cosmic:
        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(input=file_name, output=output_filename)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Cosmic Ray Removal
# ------------------------------------------------------------------------------------------------------------------- #

def cosmicray_removal(textlist_cosmic, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        textlist_cosmic : Text list containing names of FITS file to be corrected for Cosmic rays
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(textlist_cosmic)

    for file_name in list_cosmic:
        input_array, input_header = cs.fromfits(file_name)
        input_object = cs.cosmicsimage(input_array, gain=float(ccd_gain), readnoise=float(read_noise),
                                       sigclip=15.0, sigfrac=0.5, objlim=5.0, satlevel=int(data_max), verbose=False)
        input_object.run(maxiter=2, verbose=False)

        output_filename = prefix_str + file_name
        remove_file(output_filename)
        cs.tofits(outfilename=output_filename, pixelarray=input_object.cleanarray, hdr=input_header)


def cosmicray_check(textlist_cosmic, prefix_str='cr_'):
    """
    Subtracts the cosmic ray corrected image from the original image. Both the images are clipped.
    This is performed only after cosmic ray correction has been performed on images.
    Args:
        textlist_cosmic	: Text list containing the names of FITS files to be corrected for cosmic rays
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(textlist_cosmic)

    task = iraf.images.imutil.imarith
    task.unlearn()

    for file_name in list_cosmic:
        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(operand1=file_name, op='-', operand2='c' + file_name, result=output_filename)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files', choices=['Yes', 'No'])
# ctext = eg.enterbox(msg='Enter The Common Text Of Files On Which Cosmic Ray Correction Has To Be Done?',
#                    title='Cosmic Ray Correction', default='a_*.fits')
rmv_files = True
ctext = '*.fits'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Resdiual Files From Previous Run Of This Script
# ------------------------------------------------------------------------------------------------------------------- #
if rmv_files:
    for text in ['cr_*.fits']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Applies Cosmic Ray Correction To OBJECT Images & Checks The Corrected Images For Cosmic Rays Through 'cr_*.fits'
# ------------------------------------------------------------------------------------------------------------------- #
list_crreject = 'list_crreject'
group_similar_files(list_crreject, ctext)
crmedian(list_crreject)
cosmicray_check(list_crreject)
# ------------------------------------------------------------------------------------------------------------------- #
