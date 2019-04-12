#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx-------------------------PHOTOMETRIC DATA PRE-PROCESSING----------------------xxxxxxxxxxxxxxxxxxxx #
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
object_name = '2018gj'
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
        print "Error : File " + str(text_list) + " Not Found"


def python_list_to_text_list(python_list, text_list):
    """
    Put the data from the input 'python_list' to a file 'text_list' line-wise.
    Args:
        python_list : Python_list from which data has to be read
        text_list   : Name of the text file onto which data has to be appended
    Returns:
        None
    """
    with open(str(text_list), "w") as f:
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
        with open(str(file_name), 'r') as f:
            file_list = f.read().split()
            for elements in file_list:
                list_name.append(str(elements))
    python_list_to_text_list(list_name, str(text_list))

    return list_name

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.crutil(_doprint=0)
iraf.images(_doprint=0)
iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def hedit(text_list_files, field_name, value, add_keyword='no'):
    """
    Edits the header key specified by the 'field_name' of all the FITS files in the file 'text_list_files'
    and substitutes it with 'value'.
    Args:
        text_list_files : Text list containing the names of files whose header is to be edited
        add_keyword     : Whether header keyword is to be added to the files
        field_name      : The header keyword to be edited for the above files
        value           : The new value of the header keyword
    Returns:
        None
    """
    task = iraf.images.hedit
    task.unlearn()

    task.verify = 'no'                                      # Verify Each Edit Operation?
    task.add = str(add_keyword)                             # Add Rather Than Edit Fields?
    task.show = 'no'                                        # Print Record Of Each Edit Operation?
    task.update = 'yes'                                     # Enable Updating Of The Image Header?

    task(images='@' + str(text_list_files), fields=field_name, value=value)


def ccdproc(text_list, prefix_str):
    """
    Corrects Input images in the file 'text_list'.
    Args:
        text_list  : Text list containing the names of photometric bias images
        prefix_str : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    task = iraf.noao.imred.ccdred.ccdproc
    task.unlearn()

    task.trim = 'no'	                                    # Trim The Image?
    task.fixpix = 'no'					    # Fix Bad CCD Lines & Columns?
    task.overscan = 'no'	                            # Apply Overscan Strip Correction?
    task.darkcor = 'no'                                     # Apply Dark Count Correction?
    task.flatcor = 'no'					    # Apply Flat Field Correction?
    task.ccdtype = ''                                       # CCD Image Type To Combine
    task.minreplace = 1					    # Minimum Flat-Field Value?




    remove_file(str(master_bias))
    task(input='@' + str(text_list_bias), output=str(master_bias))


def zero_combine(text_list_bias, master_bias):
    """
    Combines Bias images in the file 'text_list_bias' to form a master_bias image.
    Args:
        text_list_bias  : Text list containing the names of photometric bias images
        master_bias     : Name of the master bias file
    Returns:
        None
    """
    task = iraf.noao.imred.ccdred.zerocombine
    task.unlearn()

    task.combine = 'median'                                 # Type Of Combine Operation
    task.reject = 'avsigclip'                               # Type Of Rejection
    task.ccdtype = ''                                       # CCD Image Type To Combine
    task.process = 'no'                                     # Process Images Before Combining?
    task.delete = 'no'                                      # Delete Input Images After Combining?
    task.rdnoise = float(read_noise)                        # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                             # CCD Gain (In e-/ADU)

    remove_file(str(master_bias))
    task(input='@' + str(text_list_bias), output=str(master_bias))


def bias_subtract(text_list_tbs, master_bias, prefix_str='bs_'):
    """
    Subtracts the master_bias image from the files in the file 'text_list_tbs'.
    Args:
        text_list_tbs : Text list containing the filenames from which bias is to be subtracted
        master_bias   : Name of the master bias file
        prefix_str    : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_files = text_list_to_python_list(str(text_list_tbs))

    task = iraf.images.imutil.imarith
    task.unlearn()

    for file_name in list_files:
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(operand1=str(file_name), op='-', operand2=str(master_bias), result=str(output_file_name))


def image_statistics(ctext, text_list_bsflat_mean):
    """
    Gives the statistics of the images to be analysed from the file 'text_list_stat'.
    Args:
        ctext                  : Common text of the bias subtracted flat images for which mean is to be computed
        text_list_bsflat_mean  : Text list in which mean values for the FLAT frames are to be printed out
    Returns:
        list_bsflat_mean       : Python list of all the mean values for the FLAT frames
    """
    text_list_bsflat = "list_bsflat"
    group_similar_files(str(text_list_bsflat), common_text=str(ctext))

    task = iraf.images.imutil.imstat
    task.unlearn()

    task.lower = 'INDEF'                                # Lower Limit For Pixel Values
    task.upper = 'INDEF'                                # Upper Limit For Pixel Values
    task.cache = 'no'                                   # Cache Image In Memory?

    task(format='no', images='@' + str(text_list_bsflat), field='mean', Stdout=str(text_list_bsflat_mean))
    task(format='yes', images='@' + str(text_list_bsflat), field='image, mean, stddev, min, max',
         Stdout=str(text_list_bsflat_mean).rstrip('mean') + 'stat')
    remove_file(str(text_list_bsflat))
    list_bsflat_mean = text_list_to_python_list(str(text_list_bsflat_mean))

    return list_bsflat_mean


def flat_normalise(ctext, text_list_bsflat_mean, prefix_str='n'):
    """
    Normalises the flat images based on the mean value mentioned in the file 'text_list_bsflat_mean'.
    Args:
        ctext                  : Common text of the FLAT frames which have to be normalised
        text_list_bsflat_mean  : Text list of the mean counts of the above files
        prefix_str             : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_bsflat = group_similar_files("", common_text=str(ctext))
    list_bsflat_mean = image_statistics(str(ctext), text_list_bsflat_mean=str(text_list_bsflat_mean))

    task = iraf.images.imutil.imarith
    task.unlearn()

    for i in range(0, len(list_bsflat)):
        output_file_name = str(prefix_str) + str(list_bsflat[i])
        remove_file(str(output_file_name))
        task(operand1=str(list_bsflat[i]), op='/', operand2=float(list_bsflat_mean[i]), result=str(output_file_name))


def flat_combine(ctext):
    """
    Combines flat images from the file 'text_list_nbsflat'.
    Args:
        ctext    : Common text of the normalised bias subtracted flat images of the same filter
    Returns:
        None
    """
    group_similar_files("list_nbsflat", common_text=str(ctext))
    filter_name = str(ctext)[-7]

    task = iraf.noao.imred.ccdred.flatcombine
    task.unlearn()

    task.combine = 'median'                             # Type Of Combine Operation
    task.reject = 'avsigclip'                           # Type Of Rejection
    task.ccdtype = ''                                   # CCD Image Type To Combine
    task.process = 'no'                                 # Process Images Before Combining?
    task.delete = 'no'                                  # Delete Input Images After Combining?
    task.rdnoise = float(read_noise)                    # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                         # CCD Gain (In e-/ADU)
    task.mclip = 'yes'                                  # Use Median In Sigma Clipping Algorithms?
    task.blank = 1                                      # Value If There Are No Pixels

    output_file_name = 'mflat' + str(filter_name) + '.fits'
    remove_file(str(output_file_name))
    task(input="@list_nbsflat", output=str(output_file_name))
    remove_file("list_nbsflat")


def flat_correction(ctext, master_flat, prefix_str='f'):
    """
    Corrects for flat fielding errors in the OBJECT image.
    Args:
        ctext         : Common text of the bias subtracted object images of the same filter
        master_flat   : Name of the master flat image needed for flat correction
        prefix_str    : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_bsobject = group_similar_files("", common_text=str(ctext))

    task = iraf.images.imutil.imarith
    task.unlearn()

    for image in list_bsobject:
        output_file_name = str(prefix_str) + str(image)
        remove_file(str(output_file_name))
        task(operand1=str(image), op='/', operand2=str(master_flat), result=str(output_file_name))


def cosmic_rays(ctext, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image.
    Args:
        ctext       : Common text of all the files to be corrected for cosmic rays
        prefix_str  : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_crreject = group_similar_files("", common_text=str(ctext))

    task = iraf.noao.imred.crutil.cosmicrays
    task.unlearn()

    for image in list_crreject:
        output_file_name = str(prefix_str) + str(image)
        remove_file(str(output_file_name))
        task(input=str(image), output=str(output_file_name))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Cosmic Ray Removal
# ------------------------------------------------------------------------------------------------------------------- #

def cosmicray_removal(text_list_cosmic, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        text_list_cosmic : Text list containing names of FITS file to be corrected for Cosmic rays
        prefix_str    : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(str(text_list_cosmic))

    for file_name in list_cosmic:
        input_array, input_header = cs.fromfits(str(file_name))
        input_object = cs.cosmicsimage(input_array, gain=float(ccd_gain), readnoise=float(read_noise),
                                       sigclip=15.0, sigfrac=0.5, objlim=5.0, satlevel=int(data_max), verbose=False)
        input_object.run(maxiter=2, verbose=False)

        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        cs.tofits(outfilename=str(output_file_name), pixelarray=input_object.cleanarray, hdr=input_header)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files', choices=['Yes', 'No'])
rmv_files = True
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Resdiual Files From Previous Run Of This Script
# ------------------------------------------------------------------------------------------------------------------- #
if rmv_files:
    for text in ['tmp*', 'fbs_*', 'nbs_*', 'bs_', 'mbias*', 'mflat*']:
        remove_similar_files(common_text=str(text))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Group Similar Type Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
group_similar_files("list_bias", 'bias*.fits', exceptions='spec')
group_similar_files("list_flat", 'flat*.fits')
group_similar_files("list_object", object_name + '-*.fits')
group_similar_files("list_standard", 'pg*.fits')

group_similar_files("list_u", '*-u*.fits')
group_similar_files("list_b", '*-b*.fits')
group_similar_files("list_v", '*-v*.fits')
group_similar_files("list_r", '*-r*.fits')
group_similar_files("list_i", '*-i*.fits')

list_list_phot = ["list_bias", "list_flat", "list_object", "list_standard"]
list_lists_to_list(list_list_phot[1:], "list_phot")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Edit The 'OBJECT' Keyword In The Header Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
hedit("list_flat", 'OBJECT', 'FLAT')
hedit("list_bias", 'OBJECT', 'BIAS')
hedit("list_object", 'OBJECT', str(object_name))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Edit The 'FILTER' Keyword In The Header Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
hedit("list_u", 'FILTER', '7BesU')
hedit("list_b", 'FILTER', '6BesB')
hedit("list_v", 'FILTER', '5BesV')
hedit("list_r", 'FILTER', '4BesR')
hedit("list_i", 'FILTER', '3BesI')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Edit The 'GRISM' Keyword In The Header Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
hedit("list_phot", 'GRISM', '8Free')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Combines The BIAS Frames To Obtain MASTER_BIAS And Applies BIAS Correction To The Files In The list 'list_tbs'
# ------------------------------------------------------------------------------------------------------------------- #
zero_combine(text_list_bias="list_bias", master_bias='mbias.fits')
bias_subtract(text_list_tbs="list_phot", master_bias='mbias.fits')
bias_subtract(text_list_tbs="list_standard", master_bias='mbias.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Normalises The FLAT Frames
# ------------------------------------------------------------------------------------------------------------------- #
flat_normalise(ctext='bs_flat*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Combines Normalised FLAT Frames To Form Master FLATS In Each Of The FILTERs
# ------------------------------------------------------------------------------------------------------------------- #
flat_combine(ctext='nbs_flat-u*.fits')
flat_combine(ctext='nbs_flat-b*.fits')
flat_combine(ctext='nbs_flat-v*.fits')
flat_combine(ctext='nbs_flat-r*.fits')
flat_combine(ctext='nbs_flat-i*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Applies FLAT Correction To OBJECT Images In Each Of The FILTERs
# -------------------------------------------------------------------------------------------------------------------
flat_correction(ctext="bs_*-u*.fits", master_flat='mflatu.fits')
flat_correction(ctext="bs_*-b*.fits", master_flat='mflatb.fits')
flat_correction(ctext="bs_*-v*.fits", master_flat='mflatv.fits')
flat_correction(ctext="bs_*-r*.fits", master_flat='mflatr.fits')
flat_correction(ctext="bs_*-i*.fits", master_flat='mflati.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Applies Cosmic Ray Correction To OBJECT Images
# # ------------------------------------------------------------------------------------------------------------------- #
# cosmic_rays(ctext="fbs_*.fits")
# # ------------------------------------------------------------------------------------------------------------------- #&
