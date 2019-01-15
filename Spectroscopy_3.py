#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx-------------------------SPECTROSCOPIC DATA PRE-PROCESSING----------------------xxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import math
import ephem
import shutil
import datetime
import cosmics as cs
from pyraf import iraf
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = "Indian Astronomical Observatory, Hanle"
OBS_LONG = '78:57:51'
OBS_LAT = '32:46:46'
OBS_ALT = 4500
OBS_TIMEZONE = +5.5
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
RA_keyword = 'RA'
DEC_keyword = 'DEC'
date_keyword = 'DATE-OBS'
grism_keyword = 'GRISM'
filter_keyword = 'IFILTER'
object_keyword = 'OBJECT'
airmass_keyword = 'AIRMASS'
exptime_keyword = 'EXPTIME'
time_start_keyword = 'TM_START'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
name_object = 'ASASSN14dq'
RA_object = '21:57:59.9'
DEC_object = '+24:16:08.1'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Standard Stars Details
# ------------------------------------------------------------------------------------------------------------------- #
mag_Feige110 = 11.83
RA_Feige110 = '23:19:58.4'
DEC_Feige110 = '-05:09:56.2'
mag_Feige34 = 11.18
RA_Feige34 = '10:39:36.7'
DEC_Feige34 = '+43:06:09.3'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Location Of AIRMASS File, Extinction File, Calibration Directory
# ------------------------------------------------------------------------------------------------------------------- #
FILE_AIRMASS = '/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/Airmass.asc'
FILE_EXTINCTION = '/iraf/iraf/noao/lib/onedstds/iaoextinct.dat'
DIR_CALIBRATION = '/iraf/iraf/noao/lib/onedstds/spec50cal/'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Clipping Wavelengths(In Angstroms) For Different Grisms
# ------------------------------------------------------------------------------------------------------------------- #
lower_clip_gr7 = 3750
upper_clip_gr7 = 6200
lower_clip_gr8 = 6100
upper_clip_gr8 = 9150
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.images(_doprint=0)
iraf.astutil(_doprint=0)
iraf.crutil(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.apextract(_doprint=0)
iraf.onedspec(_doprint=0)
iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
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


def check_ifexists(path):
    """
    Checks if a file or directory exists.
    Args:
        path        : Path whose existence is to be checked
    Returns:
        True        : Returns True only when the path exists
    """
    if path[-1] == '/':
        if os.path.exists(str(path)):
            return True
            pass
        else:
            print ("Error : Directory '" + str(path) + "' Does Not Exist")

    else:
        if os.path.isfile(str(path)):
            return True
            pass
        else:
            print ("Error : File '" + str(path) + "' Does Not Exist")

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

    task.verify = 'no'  # Verify Each Edit Operation?
    task.add = str(add_keyword)  # Add Rather Than Edit Fields?
    task.show = 'no'  # Print Record Of Each Edit Operation?
    task.update = 'yes'  # Enable Updating Of The Image Header?

    task(images='@' + str(text_list_files), fields=field_name, value=value)


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

    task.combine = 'median'                                     # Type Of Combine Operation
    task.reject = 'avsigclip'                                   # Type Of Rejection
    task.ccdtype = ''                                           # CCD Image Type To Combine
    task.process = 'no'                                         # Process Images Before Combining?
    task.delete = 'no'                                          # Delete Input Images After Combining?
    task.rdnoise = float(read_noise)                            # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                                 # CCD Gain (In e-/ADU)

    remove_file(str(master_bias))
    task(input='@' + str(text_list_bias), output=str(master_bias))


def bias_subtract(text_list_tbs, master_bias, prefix_str='bs_'):
    """
    Subtracts the master_bias image from the files in the file 'text_list_tbs'.
    Args:
        text_list_tbs : Text list containing the filenames from which bias is to be subtracted
        master_bias   : Name of the master bias file
        prefix_str : Prefix to distinguish the aligned FITS file from the original FITS file
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


def imcopy_clip(text_list_clip, clip_section, prefix_str='c'):
    """
    Clips the file 'file_name' based on the string 'clip_section'.
    Args:
        text_list_clip : Text list containing the names of the FITS files to be clipped
        clip_section   : Python list containing the section of the FITS file to be copied
        prefix_str  : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_clip = text_list_to_python_list(str(text_list_clip))

    task = iraf.images.imutil.imcopy
    task.unlearn()

    for file_name in list_clip:
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(input=str(file_name) + str(clip_section), output=str(output_file_name), verbose='no')


def crmedian(text_list_cosmic, clip_section, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        text_list_cosmic  : Text list containing the names of FITS files to be corrected for cosmic rays
        clip_section      : Python list containing the section of the FITS file to be copied
        prefix_str     : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(str(text_list_cosmic))
    imcopy_clip(text_list_clip=str(text_list_cosmic), clip_section=str(clip_section), prefix_str='')

    task = iraf.noao.imred.crutil.crmedian
    task.unlearn()

    task.lsigma = 25                                            # Low Clipping Sigma Factor
    task.ncsig = 10                                             # Column Box Size For Sigma Calculation

    for file_name in list_cosmic:
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(input=str(file_name), output=str(output_file_name))


def cosmicray_check(text_list_cosmic, prefix_str='cr_'):
    """
    Subtracts the cosmic ray corrected image from the original image. Both the images are clipped.
    This is performed only after cosmic ray correction has been performed on images.
    Args:
        text_list_cosmic  : Text list containing the names of FITS files to be corrected for cosmic rays
        prefix_str     : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_cosmic = text_list_to_python_list(str(text_list_cosmic))

    task = iraf.images.imutil.imarith
    task.unlearn()

    for file_name in list_cosmic:
        output_file_name = str(prefix_str) + str(file_name)[3:]
        remove_file(str(output_file_name))
        task(operand1=str(file_name), op='-', operand2='c' + str(file_name), result=str(output_file_name))


def apall_object(text_list_extract):
    """
    Extracts 1-D Spectra for file 'file_name' from given 2-D Spectra. Basically, this runs APALL task in IRAF.
    Args:
        text_list_extract : Text list containing the names of FITS files from which 1-D Spectra is to be extracted
    Returns:
        None
    """
    list_extract = text_list_to_python_list(str(text_list_extract))

    task = iraf.noao.twodspec.apextract.apall
    task.unlearn()

    task.format = 'multispec'                                   # Extracted Spectra Format
    task.interactive = 'yes'                                    # Run Task Interactively?
    task.extras = 'yes'                                         # Extract Sky, Sigma, etc.?
    task.trace = 'yes'                                          # Trace Apertures?
    task.t_function = 'spline3'                                 # Trace Fitting Function
    task.t_order = 4                                            # Trace Fitting Function Order
    task.background = 'median'                                  # Background To Subtract
    task.clean = 'yes'                                          # Detect And Replace Bad Pixels?
    task.nfind = 1                                              # Box Car Smoothing Length For Sky
    task.pfit = 'fit1d'                                         # Background Profile Fitting Type
    task.saturation = int(data_max)                             # Saturation Level For The CCD
    task.readnoise = float(read_noise)                          # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                                 # CCD Gain (In e-/ADU)
    task.lsigma = 4                                             # Lower Rejection Threshold
    task.usigma = 4                                             # Higher Rejection Threshold

    if len(list_extract) != 0:
        for file_name in list_extract:
            output_file_name = str(file_name).rstrip('.fits') + '.ms.fits'
            remove_file(str(output_file_name))
            task(input=str(file_name))
    else:
        print ("Error : No Object Spectra In '{0}' Available To Extract".format(str(text_list_extract)))


def apall_lamp(text_list_extract, reference_spec):
    """
    Extracts 1-D Spectra for lamps from given 2-D Spectra using reference 1-D Spectra 'reference_name'.
    Basically, this runs APALL task in IRAF.
    Args:
        text_list_extract : Text list containing the names of FITS files from which 1-D Spectra is to be extracted
        reference_spec    : Reference file to be used in extracting 1-D Spectra
    Returns:
        None
    """
    list_extract = text_list_to_python_list(str(text_list_extract))

    task = iraf.noao.twodspec.apextract.apall
    task.unlearn()

    task.format = 'multispec'                                   # Extracted Spectra Format
    task.interactive = 'no'                                     # Run Task Interactively?
    task.extras = 'no'                                          # Extract Sky, Sigma, etc.?
    task.trace = 'no'                                           # Trace Apertures?
    task.background = 'none'                                    # Background To Subtract
    task.clean = 'no'                                           # Detect And Replace Bad Pixels?
    task.nfind = 1                                              # Box Car Smoothing Length For Sky
    task.pfit = 'fit1d'                                         # Profile Fitting Type
    task.saturation = int(data_max)                             # Saturation Level For The CCD
    task.readnoise = float(read_noise)                          # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                                 # CCD Gain (In e-/ADU)
    task.lsigma = 4                                             # Lower Rejection Threshold
    task.usigma = 4                                             # Higher Rejection Threshold

    if len(list_extract) != 0:
        for file_name in list_extract:
            output_file_name = str(file_name).rstrip('.fits') + '.ms.fits'
            remove_file(str(output_file_name))
            task(input=str(file_name), references=str(reference_spec))
    else:
        print ("Error : No Lamp Spectra In '{0}' Available To Extract".format(str(text_list_extract)))


def dispcor(text_list_wcalib, prefix_str='w'):
    """
    Corrects for dispersion and resamples Spectra. Basically, this task perform wavelength calibration.
    Args:
        text_list_wcalib : Text list containing the names of FITS files to be wavelength calibrated
        prefix_str       : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_wcalib = text_list_to_python_list(str(text_list_wcalib))

    task = iraf.noao.onedspec.dispcor
    task.unlearn()

    for file_name in list_wcalib:
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(input=str(file_name), output=str(output_file_name))


def standard(file_name, prefix_str='std_'):
    """
    Adds standard stars to sensitivity file. Basically, runs STANDARD task in IRAF.
    Args:
        file_name     : Input 1-D standard star Spectra FITS file
        prefix_str    : Prefix to distinguish the standard star FITS file from the original FITS file
    Returns:
        output_file_name : Name of the output flux file
    """
    file_header = fits.getheader(file_name)
    airmass = file_header[str(airmass_keyword)]
    exptime = file_header[str(exptime_keyword)]
    object_name = file_header[str(object_keyword)]
    grism_name = file_header[str(grism_keyword)]

    if object_name == 'Feige110':
        star_name = 'feige110'
        mag = mag_Feige110
    elif object_name == 'Feige34':
        star_name = 'feige34'
        mag = mag_Feige34
    else:
        print ("Error : Unknown Standard Star '{0}' Chosen For Flux Calibration".format(str(object_name)))
        sys.exit(1)

    if grism_name == '4Grism7':
        grism_suffix = '-gr7'
    elif grism_name == '3Grism8':
        grism_suffix = '-gr8'
    else:
        print ("Error : Grism Keyword '{0}' Not Identified".format(str(grism_name)))
        sys.exit(1)

    task = iraf.noao.onedspec.standard
    task.unlearn()

    task.extinction = FILE_EXTINCTION                           # Location Of Extinction File
    task.caldir = DIR_CALIBRATION                               # Directory Containing Calibration Data
    task.observatory = 'iao'                                    # Observatory Of Data
    task.interact = 'no'                                        # Graphic Interaction To Define New Bandpasses?
    task.airmass = float(airmass)                               # Airmass Value
    task.exptime = float(exptime)                               # Exposure Time In Seconds
    task.mag = float(mag)                                       # Magnitude Of Star
    task.magband = 'V'                                          # Magnitude Band

    output_file_name = str(prefix_str) + str(star_name) + str(grism_suffix)
    remove_file(str(output_file_name))
    task(input=str(file_name), star_name=str(star_name), output=str(output_file_name))

    return output_file_name


def sensfunc(file_name, prefix_str='sens_'):
    """
    Determines sensitivity and extinction functions. Basically, runs SENSFUNC task in IRAF.
    Args:
        file_name    : Input 1-D standard star Spectra FITS file
        prefix_str   : Prefix to distinguish the sensitivity function FITS file from the original FITS file
    Returns:
        None
    """
    task = iraf.onedspec.sensfunc
    task.unlearn()

    task.extinction = FILE_EXTINCTION                           # Extinction File
    task.newextinction = "extinct.dat"                          # Output Revised Extinction File
    task.observatory = 'iao'                                    # Observatory Of Data
    task.function = 'spline3'                                   # Fitting Function
    task.order = 6                                              # Order Of Fit
    task.interactive = 'no'                                     # Switch On Interactive Mode?

    output_file_name = str(prefix_str) + str(file_name)[4:]
    remove_similar_files(str(output_file_name) + '*')
    task(standards=str(file_name), sensitivity=str(output_file_name))


def calibrate(ctext, prefix_str='f'):
    """
    Applies extinction corrections and flux calibrations. Basically, runs CALIBRATE task in IRAF.
    Args:
        ctext       : Common text of the FITS files to be flux calibrated
        prefix_str  : Prefix to distinguish the flux calibrated FITS file from the original FITS file
    Returns:
        None
    """
    list_files = group_similar_files("", common_text=ctext)
    
    for file_name in list_files:
        file_header = fits.getheader(file_name)
        airmass = file_header[str(airmass_keyword)]
        exptime = file_header[str(exptime_keyword)]
        grism_name = file_header[str(grism_keyword)]
    
        sens_func = ''
        if grism_name == '4Grism7':
            sens_func = glob.glob('sens_feige*gr7*')[0].rstrip('.0001.fits')
        elif grism_name == '3Grism8':
            sens_func = glob.glob('sens_feige*gr8*')[0].rstrip('.0001.fits')
        else:
            print ("Error : Grism Keyword '{0}' Not Identified".format(str(grism_name)))
    
        task = iraf.onedspec.calibrate
        task.unlearn()
    
        task.extinction = FILE_EXTINCTION                           # Extinction File
        task.observatory = 'iao'                                    # Observatory Of Data
        task.extinct = 'yes'                                        # Apply Extinction Correction?
        task.flux = 'yes'                                           # Apply Flux Correction?
        task.airmass = float(airmass)                               # Airmass Value
        task.exptime = float(exptime)                               # Exposure Time In Seconds
    
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(input=str(file_name), sens=str(sens_func), output=str(output_file_name))


def scopy_sens(ctext):
    """
    Copies the clipped version (based on wavelength) of the file 'file_name' onto the same file.
    To be done on wavelength calibrated Spectra.
    Args:
        ctext   : Common text of the FITS files to be clipped and copied
    Returns:
        None
    """
    list_copy = group_similar_files("", common_text=str(ctext))

    for file_name in list_copy:
        file_header = fits.getheader(file_name)
        grism_name = file_header[str(grism_keyword)]
    
        if grism_name == '4Grism7':
            lower_limit = lower_clip_gr7
            upper_limit = 7800
        elif grism_name == '3Grism8':
            lower_limit = 5900
            upper_limit = upper_clip_gr8
        else:
            print ("Error : Grism Keyword '{0}' Not Identified".format(str(grism_name)))
            sys.exit(1)
    
        task = iraf.onedspec.scopy
        task.unlearn()
    
        task.clobber = 'yes'                                    # Modify existing output images?
        task.format = 'multispec'                               # Output Spectra format
        task.w1 = int(lower_limit)                              # Starting Wavelength
        task.w2 = int(upper_limit)                              # Ending Wavelength
        task.verbose = 'no'                                     # Print Operations
    
        task(input=str(file_name), output=str(file_name))


def scopy_fcb(ctext, prefix_str='c'):
    """
    Copies the clipped version (based on wavelength) of the file 'file_name' onto a new file.
    Args:
        ctext       : Common text of the FITS files to be clipped and copied
        prefix_str  : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_copy = group_similar_files("", common_text=str(ctext))

    for file_name in list_copy:
        file_header = fits.getheader(file_name)
        grism_name = file_header[str(grism_keyword)]
    
        if grism_name == '4Grism7':
            lower_limit = lower_clip_gr7
            upper_limit = upper_clip_gr7
        elif grism_name == '3Grism8':
            lower_limit = lower_clip_gr8
            upper_limit = upper_clip_gr8
        else:
            print ("Error : Grism Keyword '{0}' Not Identified".format(str(grism_name)))
            sys.exit(1)
    
        task = iraf.onedspec.scopy
        task.unlearn()
    
        task.clobber = 'yes'                                        # Modify existing output images?
        task.format = 'multispec'                                   # Output Spectra format
        task.w1 = int(lower_limit)                                  # Starting Wavelength
        task.w2 = int(upper_limit)                                  # Ending Wavelength
        task.verbose = 'no'                                         # Print Operations
    
        output_file_name = str(prefix_str) + str(file_name)
        remove_file(str(output_file_name))
        task(input=str(file_name), output=str(output_file_name))


def scombine(ctext):
    """
    Combines the Spectra in the two grisms (Grism7 & Grism8) based on a common sampling region.
    Args:
        ctext   : Common text of 1-D Spectra to be combined
    Returns:
        None
    """
    text_list_comb = "list_cfwcbs"
    list_comb = group_similar_files(str(text_list_comb), common_text=str(ctext))
    sampling_region = str(lower_clip_gr8 + 25) + ':' + str(upper_clip_gr7 - 25)

    task = iraf.onedspec.scombine
    task.unlearn()

    task.combine = 'median'                                     # Type of combine operation
    task.reject = 'avsigclip'                                   # Type of rejection
    task.sample = sampling_region                               # Print Operations
    task.rdnoise = float(read_noise)                            # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                                 # CCD Gain (In e-/ADU)

    output_file_name = list_comb[0].split('-')[0] + '.ms.fits'
    remove_file(str(output_file_name))
    task(input='@' + str(text_list_comb), output=str(output_file_name))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Cosmic Ray Removal
# ------------------------------------------------------------------------------------------------------------------- #

def cosmicray_removal(text_list_cosmic, clip_section, prefix_str='c'):
    """
    Corrects for cosmic rays in the OBJECT image after clipping based on the string 'clip_section'
    Args:
        text_list_cosmic : Text list containing names of FITS file to be corrected for Cosmic rays
        clip_section     : String list containing the section of the FITS file to be copied
        prefix_str    : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    imcopy_clip(text_list_clip=str(text_list_cosmic), clip_section=clip_section, prefix_str='')
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
# Functions For Airmass Correction
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_airmass(ctext):
    """
    Calculates AIRMASS for the list of all FITS files and appends respective details in the headers.
    Args:
        ctext   : Common text of all FITS files whose headers have to be edited
    Returns:
        None
    """
    list_files = group_similar_files("", common_text=ctext)

    for file_name in list_files:
        hdulist = fits.open(file_name, mode='update')
        file_header = hdulist[0].header

        object_string = file_header[str(object_keyword)]
        date_obs = file_header[str(date_keyword)]
        time_start = file_header[str(time_start_keyword)]

        object_ra = RA_object
        object_dec = DEC_object

        if object_string == 'Feige110':
            object_ra = RA_Feige110
            object_dec = DEC_Feige110
        elif object_string == 'Feige34':
            object_ra = RA_Feige34
            object_dec = DEC_Feige34

        if str(RA_keyword) in file_header:
            file_header.set(str(RA_keyword), object_ra)
        else:
            file_header.append((str(RA_keyword), object_ra))

        if str(DEC_keyword) in file_header:
            file_header.set(str(DEC_keyword), object_dec)
        else:
            file_header.append((str(DEC_keyword), object_dec))

        time_obs = str(datetime.timedelta(seconds=int(time_start)))
        time_utc = str(date_obs) + ' ' + str(time_obs)
        julian_day = ephem.julian_date(time_utc)

        telescope = ephem.Observer()
        telescope.lon = OBS_LONG
        telescope.lat = OBS_LAT
        telescope.elevation = OBS_ALT
        telescope.pressure = 0
        telescope.epoch = ephem.J2000
        telescope.date = time_utc

        obj_pos = ephem.FixedBody()
        obj_pos._ra = object_ra
        obj_pos._dec = object_dec
        obj_pos._epoch = ephem.J2000
        obj_pos.compute(telescope)

        time_sidereal = telescope.sidereal_time()
        object_alt = Angle(str(obj_pos.alt) + ' degrees').degree
        airmass = 1 / math.cos(math.radians(90 - object_alt))

        list_keywords = ['OBSERVAT', 'OBS_LAT', 'OBS_LONG', 'OBS_ALT', 'TIMEZONE', 'DATE_OBS', 'UT', 'JD', 'ST', 'RA',
                         'DEC', 'ALT', 'AZ', 'AIRMASS']

        dict_header = {'OBSERVAT': str(OBS_NAME), 'OBS_LAT': str(OBS_LAT), 'OBS_LONG': str(OBS_LONG),
                       'OBS_ALT': str(OBS_ALT), 'TIMEZONE': str(OBS_TIMEZONE), 'DATE_OBS': str(date_obs),
                       'UT': str(time_utc), 'JD': str(julian_day), 'ST': str(time_sidereal), 'RA': str(object_ra),
                       'DEC': str(object_dec), 'ALT': str(obj_pos.alt), 'AZ': str(obj_pos.az), 'AIRMASS': str(airmass)}

        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(str(keyword), remove_all=True)
            file_header.append(card=(keyword, dict_header[keyword]))

        hdulist.flush()
        hdulist.close()


def edit_header(file_name):
    """
    Edits the header of the file 'file_name' based on the object_keyword.
    Args:
        file_name : FITS file whose header has to be edited
    Returns:
        None
    """
    (file_data, file_header) = fits.getdata(filename=str(file_name), header=True, ext=0)
    object_value = file_header[str(object_keyword)]
    date_obs = file_header[str(date_keyword)]

    object_ra = RA_object
    object_dec = DEC_object

    if object_value == 'Feige110':
        object_ra = RA_Feige110
        object_dec = DEC_Feige110
    elif object_value == 'Feige34':
        object_ra = RA_Feige34
        object_dec = DEC_Feige34

    file_header.set('DATE_OBS', str(date_obs))
    file_header.set(str(RA_keyword), object_ra)
    file_header.set(str(DEC_keyword), object_dec)

    fits.writeto(filename=str(file_name), data=file_data, header=file_header, clobber=True)

# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Checks If Files & Directories To Be Used Exist Or Not
# # ------------------------------------------------------------------------------------------------------------------- #
# list_paths = [FILE_AIRMASS, FILE_EXTINCTION, DIR_CALIBRATION]
#
# for path in list_paths:
#     check_ifexists(path)
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Removes Residual Files From Previous Run Of Pre-Processing
# # ------------------------------------------------------------------------------------------------------------------- #
# remove_similar_files("*.ms.*")
# remove_similar_files("bs_*")
# remove_similar_files("cbs_*")
# remove_similar_files("list_*")
# remove_similar_files('*tmp*')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Groups Similar Type Of FITS Files
# # ------------------------------------------------------------------------------------------------------------------- #
# group_similar_files("list_object", "ASASSN14dq-*.fits")
# group_similar_files("list_biasspec", "biasspec*.fits")
#
# group_similar_files("list_gr7", "*gr7*.fits")
# group_similar_files("list_gr8", "*gr8*.fits")
#
# group_similar_files("list_FeNe", "FeNe*.fits")
# group_similar_files("list_FeAr", "FeAr*.fits")
# group_similar_files("list_Feige110", "Feige110*.fits")
# group_similar_files("list_Feige34", "Feige34*.fits")
#
# list_list_spec = ["list_biasspec", "list_object", "list_FeNe", "list_FeAr", "list_Feige110", "list_Feige34"]
# list_spec = list_lists_to_list(list_list_spec, "list_spec")
# list_tbs_spec = list_lists_to_list(list_list_spec[1:], "list_tbs_spec")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Edits The 'OBJECT' Keyword In The Headers Of FITS Files
# # ------------------------------------------------------------------------------------------------------------------- #
# hedit("list_biasspec", 'OBJECT', 'BIAS_SPEC')
# hedit("list_FeNe", 'OBJECT', 'FeNe')
# hedit("list_FeAr", 'OBJECT', 'FeAr')
# hedit("list_Feige110", 'OBJECT', 'Feige110')
# hedit("list_Feige34", 'OBJECT', 'Feige34')
# hedit("list_Halogen", 'OBJECT', 'Halogen')
# hedit("list_object", 'OBJECT', str(name_object))
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Edits The 'FILTER' Keyword In The Headers Of FITS Files
# # ------------------------------------------------------------------------------------------------------------------- #
# hedit("list_spec", 'FILTER', '1Free')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Edits The 'GRISM' Keyword In The Headers Of FITS Files
# # ------------------------------------------------------------------------------------------------------------------- #
# hedit("list_gr7", 'GRISM', '4Grism7')
# hedit("list_gr8", 'GRISM', '3Grism8')
# hedit("list_FeAr", 'GRISM', '4Grism7')
# hedit("list_FeNe", 'GRISM', '3Grism8')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Combines The BIAS Frames To Obtain MASTER_BIAS And Applies BIAS Correction To The Files In The list 'list_tbs'
# # ------------------------------------------------------------------------------------------------------------------- #
# zero_combine(text_list_bias="list_biasspec", master_bias="mbiasspec.fits")
# bias_subtract(text_list_tbs="list_tbs_spec", master_bias="mbiasspec.fits")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Groups BIAS Subtracted Spectra Of Objects+Standards & Lamps
# # ------------------------------------------------------------------------------------------------------------------- #
# group_similar_files("list_bs_gr", common_text="bs_*gr*.fits")
# group_similar_files("list_bs_lamp", common_text="bs_Fe*.fits", exceptions='gr')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Clips The BIAS Subtracted Spectroscopic Images And Corrects For Cosmic Ray Hits
# # ------------------------------------------------------------------------------------------------------------------- #
# # cosmicray_removal(text_list_cosmic="list_bs_gr", clip_section='[*, 51:3200]')
#
# crmedian(text_list_cosmic="list_bs_gr", clip_section='[*, 51:3200]')
#
# imcopy_clip(text_list_clip="list_bs_lamp", clip_section='[*, 51:3200]')
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Checks The Corrected Images For Cosmic Rays Through 'cr_*.fits'
# # ------------------------------------------------------------------------------------------------------------------- #
# cosmicray_check("list_bs_gr")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Groups Extracted 1-D Spectra Of Objects & Standards (Based On GRISM)
# # ------------------------------------------------------------------------------------------------------------------- #
# list_cbs_gr7_ms = group_similar_files("list_cbs_gr7_ms", common_text="cbs_*gr7*.ms.fits")
# list_cbs_gr8_ms = group_similar_files("list_cbs_gr8_ms", common_text="cbs_*gr8*.ms.fits")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Extracts 1-D Spectra (From 2-D Spectra) For Lamps
# # ------------------------------------------------------------------------------------------------------------------- #
# group_similar_files("list_cbs_FeAr", common_text="cbs_FeAr*.fits", exceptions='ms')
# group_similar_files("list_cbs_FeNe", common_text="cbs_FeNe*.fits", exceptions='ms')
#
# if len(list_cbs_gr7_ms) != 0:
#     apall_lamp("list_cbs_FeAr", reference_spec=list_cbs_gr7_ms[0].rstrip('.ms.fits'))
# else:
#     print ("Error : Reference Spectra Not Available For Grism 7")
#
# if len(list_cbs_gr8_ms) != 0:
#     apall_lamp("list_cbs_FeNe", reference_spec=list_cbs_gr8_ms[0].rstrip('.ms.fits'))
# else:
#     print ("Error : Reference Spectra Not Available For Grism 8")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Groups Extracted 1-D Spectra Of Lamps
# # ------------------------------------------------------------------------------------------------------------------- #
# list_cbs_FeAr_ms = group_similar_files("", common_text="cbs_FeAr*.ms.fits")
# list_cbs_FeNe_ms = group_similar_files("", common_text="cbs_FeNe*.ms.fits")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Appends Header Keyword 'REFSPEC1' To Spectra Of Objects & Standards With Extracted Lamp Spectra
# # ------------------------------------------------------------------------------------------------------------------- #
# if len(list_cbs_FeAr_ms) != 0:
#     hedit("list_cbs_gr7_ms", field_name='REFSPEC1', value=list_cbs_FeAr_ms[0], add_keyword='yes')
# else:
#     print ("ERROR: Header Keyword 'REFSPEC1' Cannot Be Appended -> 1-D Lamp Spectra For Grism7 Not Available")
#
# if len(list_cbs_FeNe_ms) != 0:
#     hedit("list_cbs_gr8_ms", field_name='REFSPEC1', value=list_cbs_FeNe_ms[0], add_keyword='yes')
# else:
#     print ("ERROR: Header Keyword 'REFSPEC1' Cannot Be Appended -> 1-D Lamp Spectra For Grism8 Not Available")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Performs Wavelength Calibration On 1-D Spectra Of Objects & Standards
# # ------------------------------------------------------------------------------------------------------------------- #
# group_similar_files("list_cbs_gr_ms", "cbs_*gr*.ms.fits")
# dispcor("list_cbs_gr_ms")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Calculates AIRMASS And Clips Wavelength Calibrated 1-D Spectra Of Objects & Standards
# # ------------------------------------------------------------------------------------------------------------------- #
# calculate_airmass(ctext="wcbs_*.fits")
# # scopy_sens(ctext="wcbs_*.fits")
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Determines Sensitivity Function Of The CCD Using Wavelength Calibrated Standard Star Data
# # ------------------------------------------------------------------------------------------------------------------- #
# list_wcbs_standard = group_similar_files("", common_text="wcbs_*Feige*")
#
# for file_name in list_wcbs_standard:
#     sensfunc(standard(file_name))
# # ------------------------------------------------------------------------------------------------------------------- #
#
#
# # ------------------------------------------------------------------------------------------------------------------- #
# # Performs Flux Calibration On Object Spectra Using CCD's Sensitivity Function
# # ------------------------------------------------------------------------------------------------------------------- #
# calibrate(ctext="wcbs_*ASASSN*.fits")
# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Combine Operation On Clipped Flux Calibrated Object Spectra (If Both Grism7 & Grism8 Data Exist)
# ------------------------------------------------------------------------------------------------------------------- #
list_fwcbs = group_similar_files("", common_text="fwcbs_*gr*.fits")

if len(list_fwcbs) >= 2:
    scopy_fcb(ctext="fwcbs_*gr*.fits")
    scombine(ctext="cfwcbs_*gr*.fits")
elif len(list_fwcbs) == 1:
    output_file_name = 'c' + list_fwcbs[0]
    scopy_sens(ctext="fwcbs_*.fits")
    shutil.copy(list_fwcbs[0], str(output_file_name))
else:
    print ("Error : Needs Both Grism7 & Grism8 Clipped Flux Calibrated 1-D Spectra To Combine")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes All The Text LISTS Created In The Current Directory
# ------------------------------------------------------------------------------------------------------------------- #
remove_similar_files("*list*")
# ------------------------------------------------------------------------------------------------------------------- #
