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
import math
import ephem
import shutil
import datetime
from pyraf import iraf
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = 'Indian Astronomical Observatory, Hanle'
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
DATE_keyword = 'DATE-OBS'
GRISM_keyword = 'GRISM'
OBJECT_keyword = 'OBJECT'
AIRMASS_keyword = 'AIRMASS'
EXPTIME_keyword = 'EXPTIME'
TIMESTART_keyword = 'TM_START'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_NAME = '2018cow'
OBJECT_RA = '16:16:00.22'
OBJECT_DEC = '+22:16:04.83'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Standard Stars Details
# ------------------------------------------------------------------------------------------------------------------- #
mag_Feige34 = 11.18
RA_Feige34 = '10:39:36.7'
DEC_Feige34 = '+43:06:09.3'

mag_Feige66 = 10.50
RA_Feige66 = '12:37:23.5'
DEC_Feige66 = '+25:03:59.9'

mag_Feige67 = 11.63
RA_Feige67 = '12:41:51.79'
DEC_Feige67 = '+17:31:19.75'

mag_Feige110 = 11.82
RA_Feige110 = '23:19:58.4'
DEC_Feige110 = '-05:09:56.2'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Location Of AIRMASS File, Extinction File, Calibration Directory
# ------------------------------------------------------------------------------------------------------------------- #
FILE_EXTINCTION = '/iraf/iraf/noao/lib/onedstds/iaoextinct.dat'
DIR_CALIBRATION = '/iraf/iraf/noao/lib/onedstds/spec50cal/'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Clipping Wavelengths (In Angstroms) For Different Grisms (For Copying & Combining Respectively)
# ------------------------------------------------------------------------------------------------------------------- #
lowerclip_gr7 = 3500
upperclip_gr7 = 7720
lowerclip_gr8 = 5900
upperclip_gr8 = 9150

combclip_gr7 = 6200
combclip_gr8 = 6100
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.apextract(_doprint=0)
iraf.onedspec(_doprint=0)
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
            f.write(element + '\n')


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


def check_ifexists(path):
    """
    Checks if a file or directory exists.
    Args:
        path    : Path whose existence is to be checked
    Returns:
        True    : Returns True only when the path exists
    Raises:
        Error   : Directory/File 'path' Not Found
    """
    if path[-1] == '/':
        if os.path.exists(path):
            return True
            pass
        else:
            print ("ERROR : Directory '{0}' Does Not Exist".format(path))

    else:
        if os.path.isfile(path):
            return True
            pass
        else:
            print ("ERROR : File '{0}' Does Not Exist".format(path))


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
    task = iraf.images.hedit
    task.unlearn()

    task.verify = 'no'                                      # Verify Each Edit Operation?
    task.add = add_keyword                                  # Add Rather Than Edit Fields?
    task.show = 'no'                                        # Print Record Of Each Edit Operation?
    task.update = 'yes'                                     # Enable Updating Of The Image Header?

    task(images='@' + textlist_files, fields=field_name, value=value)


def apall_object(textlist_extract):
    """
    Extracts 1-D Spectra for file 'file_name' from given 2-D Spectra. Basically, this runs APALL task in IRAF.
    Args:
        textlist_extract    : Text list containing the names of FITS files from which 1-D Spectra is to be extracted
    Returns:
        None
    Raises:
        Error               : No Object Spectra In 'textlist_extract' Available To Extract
    """
    list_extract = text_list_to_python_list(textlist_extract)

    task = iraf.noao.twodspec.apextract.apall
    task.unlearn()

    task.format = 'multispec'                               # Extracted Spectra Format
    task.interactive = 'yes'                                # Run Task Interactively?
    task.extras = 'yes'                                     # Extract Sky, Sigma, etc.?
    task.trace = 'yes'                                      # Trace Apertures?
    task.t_function = 'spline3'                             # Trace Fitting Function
    task.t_order = 4                                        # Trace Fitting Function Order
    task.background = 'median'                              # Background To Subtract
    task.clean = 'yes'                                      # Detect And Replace Bad Pixels?
    task.weights = 'variance'                               # Extraction Weights
    task.nfind = 1                                          # Box Car Smoothing Length For Sky
    task.pfit = 'fit1d'                                     # Background Profile Fitting Type
    task.saturation = int(data_max)                         # Saturation Level For The CCD
    task.readnoise = float(read_noise)                      # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                             # CCD Gain (In e-/ADU)
    task.lsigma = 4                                         # Lower Rejection Threshold
    task.usigma = 4                                         # Higher Rejection Threshold

    if len(list_extract) != 0:
        for file_name in list_extract:
            output_filename = file_name.rstrip('.fits') + '.ms.fits'
            remove_file(output_filename)
            task(input=file_name)
    else:
        print ("ERROR : No Object Spectra In '{0}' Available To Extract".format(textlist_extract))


def apall_lamp(textlist_extract, ref_spec):
    """
    Extracts 1-D Spectra for lamps from given 2-D Spectra using reference 1-D Spectra 'reference_name'.
    Basically, this runs APALL task in IRAF.
    Args:
        textlist_extract    : Text list containing the names of FITS files from which 1-D Spectra is to be extracted
        ref_spec            : Reference file to be used in extracting 1-D Spectra
    Returns:
        None
    Raises:
        Error               : No Lamp Spectra In 'textlist_extract' Available To Extract
    """
    list_extract = text_list_to_python_list(textlist_extract)

    task = iraf.noao.twodspec.apextract.apall
    task.unlearn()

    task.format = 'multispec'                               # Extracted Spectra Format
    task.interactive = 'no'                                 # Run Task Interactively?
    task.extras = 'no'                                      # Extract Sky, Sigma, etc.?
    task.trace = 'no'                                       # Trace Apertures?
    task.background = 'none'                                # Background To Subtract
    task.clean = 'no'                                       # Detect And Replace Bad Pixels?
    task.nfind = 1                                          # Box Car Smoothing Length For Sky
    task.pfit = 'fit1d'                                     # Background Profile Fitting Type
    task.weights = 'variance'                               # Extraction Weights
    task.saturation = int(data_max)                         # Saturation Level For The CCD
    task.readnoise = float(read_noise)                      # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                             # CCD Gain (In e-/ADU)
    task.lsigma = 4                                         # Lower Rejection Threshold
    task.usigma = 4                                         # Higher Rejection Threshold

    if len(list_extract) != 0:
        for file_name in list_extract:
            output_filename = file_name.rstrip('.fits') + '.ms.fits'
            remove_file(output_filename)
            task(input=file_name, references=ref_spec)
    else:
        print ("ERROR : No Lamp Spectra In '{0}' Available To Extract".format(textlist_extract))


def dispcor(ctext, prefix_str='w'):
    """
    Corrects for dispersion and resamples Spectra. Basically, this task perform wavelength calibration.
    Args:
        ctext        : Common text of the FITS files to be wavelength calibrated
        prefix_str   : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_wcalib = group_similar_files('', common_text=ctext)

    task = iraf.noao.onedspec.dispcor
    task.unlearn()

    for file_name in list_wcalib:
        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(input=file_name, output=output_filename)


def standard(file_name, prefix_str='std'):
    """
    Adds standard stars to sensitivity file. Basically, runs STANDARD task in IRAF.
    Args:
        file_name       : Input 1-D standard star Spectra FITS file
        prefix_str      : Prefix to distinguish the standard star FITS file from the original FITS file
    Returns:
        output_filename : Name of the output flux file
    Raises:
        Error           : Unknown Standard Star 'object_name' Chosen For Flux Calibration
    """
    file_header = fits.getheader(file_name)
    airmass = file_header[AIRMASS_keyword]
    exptime = file_header[EXPTIME_keyword]
    object_name = file_header[OBJECT_keyword]
    grism = file_header[GRISM_keyword]

    if object_name == 'Feige110':
        mag = mag_Feige110
    elif object_name == 'Feige34':
        mag = mag_Feige34
    elif object_name == 'Feige66':
        mag = mag_Feige66
    elif object_name == 'Feige67':
        mag = mag_Feige67
    else:
        print ("ERROR : Unknown Standard Star '{0}' Chosen For Flux Calibration".format(object_name))
        sys.exit(1)

    task = iraf.noao.onedspec.standard
    task.unlearn()

    task.extinction = FILE_EXTINCTION                       # Location Of Extinction File
    task.caldir = DIR_CALIBRATION                           # Directory Containing Calibration Data
    task.observatory = 'iao'                                # Observatory Of Data
    task.interact = 'no'                                    # Graphic Interaction To Define New Bandpasses?
    task.airmass = float(airmass)                           # Airmass Value
    task.exptime = float(exptime)                           # Exposure Time In Seconds
    task.mag = float(mag)                                   # Magnitude Of The Standard Star
    task.magband = 'V'                                      # Magnitude Band

    output_filename = prefix_str + file_name.lstrip('wcbs').rstrip('.ms.fits')
    remove_file(output_filename)
    task(input=file_name, star_name=object_name.lower(), output=output_filename)

    return output_filename


def sensfunc(file_name, prefix_str='sens'):
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

    task.extinction = FILE_EXTINCTION                       # Extinction File
    task.newextinction = 'extinct.dat'                      # Output Revised Extinction File
    task.observatory = 'iao'                                # Observatory Of Data
    task.function = 'spline3'                               # Fitting Function
    task.order = 6                                          # Order Of Fit
    task.interactive = 'no'                                 # Switch On Interactive Mode?

    output_filename = prefix_str + file_name.lstrip('std')
    remove_similar_files(output_filename + '*')
    task(standards=file_name, sensitivity=output_filename)


def calibrate(ctext, prefix_str='f'):
    """
    Applies extinction corrections and flux calibrations. Basically, runs CALIBRATE task in IRAF.
    Args:
        ctext       : Common text of the FITS files to be flux calibrated
        prefix_str  : Prefix to distinguish the flux calibrated FITS file from the original FITS file
    Returns:
        None
    Raises:
        Error       : Grism Keyword 'grism' Not Identified For File Name 'file_name'
    """
    list_files = group_similar_files('', common_text=ctext)

    for file_name in list_files:
        file_header = fits.getheader(file_name)
        airmass = file_header[AIRMASS_keyword]
        exptime = file_header[EXPTIME_keyword]
        grism = file_header[GRISM_keyword]

        if grism == '4Grism7':
            sens_func = glob.glob('sens_Feige*gr7*')[0].rstrip('.0001.fits')
        elif grism == '3Grism8':
            sens_func = glob.glob('sens_Feige*gr8*')[0].rstrip('.0001.fits')
        else:
            print ("ERROR : Grism Keyword '{0}' Not Identified For File Name '{1}'".format(grism, file_name))
            sys.exit(1)

        task = iraf.onedspec.calibrate
        task.unlearn()

        task.extinction = FILE_EXTINCTION                   # Extinction File
        task.observatory = 'iao'                            # Observatory Of Data
        task.extinct = 'yes'                                # Apply Extinction Correction?
        task.flux = 'yes'                                   # Apply Flux Correction?
        task.airmass = float(airmass)                       # Airmass Value
        task.exptime = float(exptime)                       # Exposure Time In Seconds

        output_filename = prefix_str + file_name
        remove_file(output_filename)
        task(input=file_name, sens=sens_func, output=output_filename)


def scopy(ctext, for_comb=False, prefix_str=''):
    """
    Copies the clipped version (based on wavelength) of the file 'file_name' onto the same file.
    To be done on wavelength calibrated Spectra.
    Args:
        ctext       : Common text of the FITS files to be clipped and copied
        prefix_str  : Prefix to distinguish the aligned FITS file from the original FITS file
        for_comb    : Boolean describing whether to just make copy of the spectra or copy it for combining
    Returns:
        None
    Raises:
        Error       : Grism Keyword 'grism' Not Identified For File Name 'file_name'
    """
    list_copy = group_similar_files('', common_text=ctext)

    for file_name in list_copy:
        file_header = fits.getheader(file_name)
        grism = file_header[GRISM_keyword]

        if grism == '4Grism7':
            lower_limit = lowerclip_gr7
            if for_comb:
                upper_limit = combclip_gr7
            else:
                upper_limit = upperclip_gr7
        elif grism == '3Grism8':
            upper_limit = upperclip_gr8
            if for_comb:
                lower_limit = combclip_gr8
            else:
                lower_limit = lowerclip_gr8
        else:
            print ("ERROR : Grism Keyword '{0}' Not Identified For File Name '{1}'".format(grism, file_name))
            sys.exit(1)

        task = iraf.onedspec.scopy
        task.unlearn()

        task.clobber = 'yes'                                # Modify existing output images?
        task.format = 'multispec'                           # Output Spectra format
        task.w1 = int(lower_limit)                          # Starting Wavelength
        task.w2 = int(upper_limit)                          # Ending Wavelength
        task.verbose = 'no'                                 # Print Operations

        if prefix_str != '':
            output_filename = prefix_str + file_name
            remove_file(output_filename)
        else:
            output_filename = file_name

        task(input=file_name, output=output_filename)


def scombine(ctext):
    """
    Combines the Spectra in the two grisms (Grism7 & Grism8) based on a common sampling region.
    Args:
        ctext   : Common text of 1-D Spectra to be combined
    Returns:
        None
    """
    textlist_comb = 'list_cfwcbs'
    list_comb = group_similar_files(textlist_comb, common_text=ctext)
    sampling_region = str(lowerclip_gr8 + 25) + ':' + str(upperclip_gr7 - 25)

    task = iraf.onedspec.scombine
    task.unlearn()

    task.combine = 'median'                                 # Type of combine operation
    task.reject = 'avsigclip'                               # Type of rejection
    task.sample = sampling_region                           # Wavelength Sampling Region
    task.rdnoise = float(read_noise)                        # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                             # CCD Gain (In e-/ADU)

    output_filename = list_comb[0].split('-')[0] + '.ms.fits'
    remove_file(output_filename)
    task(input='@' + textlist_comb, output=output_filename)

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
    list_files = group_similar_files('', common_text=ctext)

    for file_name in list_files:
        hdulist = fits.open(file_name, mode='update')
        file_header = hdulist[0].header

        object_str = file_header[OBJECT_keyword]
        date_obs = file_header[DATE_keyword]
        time_start = file_header[TIMESTART_keyword]

        object_ra = OBJECT_RA
        object_dec = OBJECT_DEC

        if object_str == 'Feige110':
            object_ra = RA_Feige110
            object_dec = DEC_Feige110
        elif object_str == 'Feige34':
            object_ra = RA_Feige34
            object_dec = DEC_Feige34
        elif object_str == 'Feige66':
            object_ra = RA_Feige66
            object_dec = DEC_Feige66

        if RA_keyword in file_header:
            file_header.set(RA_keyword, object_ra)
        else:
            file_header.append((RA_keyword, object_ra))

        if DEC_keyword in file_header:
            file_header.set(DEC_keyword, object_dec)
        else:
            file_header.append((DEC_keyword, object_dec))

        time_obs = str(datetime.timedelta(seconds=int(time_start)))
        time_utc = date_obs + ' ' + time_obs
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

        list_keywords = ['OBSERVAT', 'OBS_LAT', 'OBS_LONG', 'OBS_ALT', 'TIMEZONE', 'DATE_OBS', 'UT', 'JD', 'ST',
                         'RA', 'DEC', 'ALT', 'AZ', 'AIRMASS']

        dict_header = {'OBSERVAT': OBS_NAME, 'OBS_LAT': OBS_LAT, 'OBS_LONG': OBS_LONG, 'OBS_ALT': OBS_ALT,
                       'TIMEZONE': OBS_TIMEZONE, 'DATE_OBS': date_obs, 'UT': time_utc, 'JD': julian_day,
                       'ST': time_sidereal, 'RA': object_ra, 'DEC': object_dec, 'ALT': obj_pos.alt,
                       'AZ': obj_pos.az, 'AIRMASS': airmass}

        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all=True)
            file_header.append(card=(keyword, str(dict_header[keyword])))

        hdulist.flush()
        hdulist.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# remove_resfile = eg.boolbox(msg='Remove Residual Files From Previous Run Of This Script?',
#                             title='Remove Residual Files', choices=['Yes', 'No'])
remove_resfile = True
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of PostProcSpec.py
# ------------------------------------------------------------------------------------------------------------------- #
if remove_resfile:
    for text in ['fwcbs_*', 'cfwcbs_*', 'sens_*', 'std_*', 'list_*', '*tmp&']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Checks If Files & Directories To Be Used Exist Or Not
# ------------------------------------------------------------------------------------------------------------------- #
list_paths = [FILE_EXTINCTION, DIR_CALIBRATION]

for path in list_paths:
    check_ifexists(path)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Extracted 1-D Spectra Of Objects & Standards (Based On GRISM)
# Extracts 1-D Spectra (From 2-D Spectra) For Lamps
# ------------------------------------------------------------------------------------------------------------------- #
list_cbs_gr7_ms = group_similar_files('list_cbs_gr7_ms', common_text='cbs_*gr7*.ms.fits')
list_cbs_gr8_ms = group_similar_files('list_cbs_gr8_ms', common_text='cbs_*gr8*.ms.fits')

group_similar_files('list_cbs_FeAr', common_text='cbs_FeAr*.fits', exceptions='ms')
group_similar_files('list_cbs_FeNe', common_text='cbs_FeNe*.fits', exceptions='ms')

if len(list_cbs_gr7_ms) != 0:
    apall_lamp('list_cbs_FeAr', ref_spec=list_cbs_gr7_ms[0].rstrip('.ms.fits'))
else:
    print ("ERROR : Reference Object Spectra Not Available For Grism 7")

if len(list_cbs_gr8_ms) != 0:
    apall_lamp('list_cbs_FeNe', ref_spec=list_cbs_gr8_ms[0].rstrip('.ms.fits'))
else:
    print ("ERROR : Reference Object Spectra Not Available For Grism 8")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Extracted 1-D Spectra Of Lamps
# Appends Header Keyword 'REFSPEC1' To Spectra Of Objects & Standards With Extracted Lamp Spectra
# ------------------------------------------------------------------------------------------------------------------- #
list_cbs_FeAr_ms = group_similar_files('', common_text='cbs_FeAr*.ms.fits')
list_cbs_FeNe_ms = group_similar_files('', common_text='cbs_FeNe*.ms.fits')

if len(list_cbs_gr7_ms) != 0:
    if len(list_cbs_FeAr_ms) != 0:
        hedit('list_cbs_gr7_ms', field_name='REFSPEC1', value=list_cbs_FeAr_ms[0], add_keyword='yes')
    else:
        print ("ERROR : Header Keyword 'REFSPEC1' Cannot Be Appended -> 1-D Lamp Spectra For Grism7 Not Available")
        sys.exit(1)

if len(list_cbs_gr8_ms) != 0:
    if len(list_cbs_FeNe_ms) != 0:
        hedit('list_cbs_gr8_ms', field_name='REFSPEC1', value=list_cbs_FeNe_ms[0], add_keyword='yes')
    else:
        print ("ERROR : Header Keyword 'REFSPEC1' Cannot Be Appended -> 1-D Lamp Spectra For Grism8 Not Available")
        sys.exit(1)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Wavelength Calibration On 1-D Spectra Of Objects & Standards
# Calculates AIRMASS And Clips Wavelength Calibrated 1-D Spectra Of Objects & Standards
# ------------------------------------------------------------------------------------------------------------------- #
dispcor(ctext='cbs_*gr*.ms.fits')
calculate_airmass(ctext='wcbs_*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Clip The Wavelength Calibrated Spectra Before Flux Calibration (Both Standard & The Object)
# Determines Sensitivity Function Of The CCD Using Wavelength Calibrated Standard Star Data
# Performs Flux Calibration On Object Spectra Using CCD's Sensitivity Function
# ------------------------------------------------------------------------------------------------------------------- #
scopy(ctext='wcbs_*Feige*')
scopy(ctext='wcbs_' + OBJECT_NAME + '*.fits')

for file_name in group_similar_files('', common_text='wcbs_*Feige*'):
    sensfunc(standard(file_name))

calibrate(ctext='wcbs_*' + OBJECT_NAME + '*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Combine Operation On Clipped Flux Calibrated Object Spectra (If Both Grism7 & Grism8 Data Exist)
# Append The Name Of The Final 1-D Spectra With The Date Of Observation
# ------------------------------------------------------------------------------------------------------------------- #
list_fwcbs = group_similar_files('', common_text='fwcbs_*.fits')

if len(list_fwcbs) >= 2:
    scopy(ctext='fwcbs_*.fits', for_comb=True, prefix_str='c')
    scombine(ctext='cfwcbs_*gr*.fits')
elif len(list_fwcbs) == 1:
    scopy(ctext='fwcbs_*.fits')
    shutil.copy(list_fwcbs[0], list_fwcbs[0].split('-')[0] + '.ms.fits')
    shutil.copy(list_fwcbs[0], 'c' + list_fwcbs[0])
else:
    print ("ERROR : Needs Both Grism7 & Grism8 Flux Calibrated 1-D Spectra To Combine")

for file_name in group_similar_files('', common_text='cfwcbs*.fits'):
    header = fits.getheader(file_name, ext=0)
    date_obs = header[DATE_keyword]
    shutil.copy(file_name, date_obs + '_' + file_name)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes All The Text LISTS Created In The Current Directory
# ------------------------------------------------------------------------------------------------------------------- #
remove_similar_files('*list*')
# ------------------------------------------------------------------------------------------------------------------- #
