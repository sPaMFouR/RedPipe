#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxx-------------------PHOTOMETRY OF OBJECT FRAMES-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import numpy as np
import easygui as eg
from pyraf import iraf
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Home Directory & Its Component Folders
# ------------------------------------------------------------------------------------------------------------------- #
DIR_HOME = os.getcwd()
dir_text = eg.enterbox(msg='Enter The Home Directory Containing Component Folders: ', title='Home Directory',
                       default=DIR_HOME)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Telescope CCD Specifications
# ------------------------------------------------------------------------------------------------------------------- #
gain = 11.98
read_noise = 7
data_max = 50000

BackSub = False         # 'True' To Switch On Background Subtraction
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
UT_keyword = 'UT'
FILTER_keyword = 'FILTER'
AIRMASS_keyword = 'AIRMASS'
EXPTIME_keyword = 'EXPTIME'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.daophot(_doprint=0)
iraf.ptools(_doprint=0)
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
                test = re.search(text, file_name)
                if test:
                    try:
                        list_files.remove(file_name)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for index in range(0, len(list_files)):
                f.write(list_files[index] + '\n')

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
        print ("Error : File " + text_list + " Not Found")
        sys.exit(1)


def list_statistics(python_list):
    """
    Returns the statistics of the list of elements in the input 'python_list'.
    Args:
        python_list  : Input list of elements
    Returns:
        value_mean   : Mean of the list of elements
        value_median : Median of the list of elements
        value_std    : Standard Deviation of the list of elements
    """
    value_mean = np.mean(python_list)
    value_median = np.median(python_list)
    value_std = np.std(python_list)

    return value_mean, value_median, value_std


def reject(python_list):
    """
    Rejects outliers from the input 'python_list'.
    Args:
        python_list : Input list of elements
    Returns:
        reject_list : Output list after rejecting outliers from the input 'python_list'
    """
    reject_list = []
    pop = False
    for index in range(0, len(python_list)):
        reject_list.append(float(python_list[index]))

    reject_list.sort()
    value_mean, value_median, value_std = list_statistics(reject_list)

    if abs(reject_list[0] - value_median) < abs(reject_list[-1] - value_median):
        remove_index = -1
    else:
        remove_index = 0

    if abs(reject_list[remove_index] - value_median) > value_std:
        reject_list.pop(remove_index)
        pop = True

    if pop:
        value_mean, value_median, value_std = list_statistics(reject_list)
        if abs(reject_list[0] - value_median) < abs(reject_list[-1] - value_median):
            remove_index = -1
        else:
            remove_index = 0
        if abs(reject_list[remove_index] - value_median) > 2 * value_std:
            reject_list.pop(remove_index)

    return reject_list


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


def display_text(text_to_display):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text_to_display : Text to be displayed
    Returns:
        None
    """
    print ("\n" + "# " + "-" * (12 + len(text_to_display)) + " #")
    print ("# " + "-" * 5 + " " + text_to_display + " " + "-" * 5 + " #")
    print ("# " + "-" * (12 + len(text_to_display)) + " #" + "\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def imexam_fwhm(text_list, coord_file, log_imexam='log_imexam'):
    """
    Examines the images in the list 'python_list' at the coordinates mentioned in the file "stars.coo"
    and logs the output onto the file "log_imexam".
    Args:
        text_list    : Text list containing names of FITS files
        coord_file   : Text file listing the coordinates of selected stars in the field
        log_imexam   : Name of the text list to record log of IMEXAM
    Returns:
        None
    """
    remove_file(log_imexam)
    list_files = text_list_to_python_list(text_list)

    task = iraf.images.tv.imexam
    task.unlearn()

    for file_name in list_files:
        task.logfile = str(log_imexam)              # Log File To Record Output Of The Commands
        task.keeplog = 'yes'                        # Log Output Results?
        task.defkey = 'a'                           # Default Key For Cursor x-y Input List
        task.imagecur = str(coord_file)             # Image Display Cursor Input
        task.use_display = 'no'                     # Use The Image Display?

        task(input=str(file_name), frame=1)


def data_pars(fwhm_value, data_max):
    """
    Edits the data dependent parameters(DATAPARS) required by the DAOPHOT tasks.
    Args:
        fwhm_value  : Mean FWHM value for the image file
        data_max    : Maximum good pixel value
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.datapars
    task.unlearn()

    task.scale = 1.0                                # Scale Of The Image In Arcseconds Per Pixel
    task.fwhmpsf = float(fwhm_value)                # FWHM Of The PSF In Scale Units
    task.emission = 'yes'                           # All Features Are Considered To Be Emission Features
    task.datamin = 'INDEF'                          # Minimum Good Pixel Value
    task.datamax = data_max                         # Maximum Good Pixel Value
    task.noise = 'poisson'                          # Noise Model Used To Estimate Uncertainties In APPHOT Magnitudes
    task.sigma = 'INDEF'                            # Standard Deviation Of The Sky Pixels
    task.readnoise = read_noise                     # Readout Noise Of The CCD In Electrons
    task.epadu = gain                               # Gain Of The CCD In Electrons Per ADU
    task.exposure = EXPTIME_keyword                 # Exposure Time Keyword In Image Header
    task.airmass = AIRMASS_keyword                  # Airmass Keyword In Image Header
    task.filter = FILTER_keyword                    # Filter Keyword In Image Header
    task.obstime = UT_keyword                       # UT Keyword In Image Header


def center_pars():
    """
    Edits the centering algorthm parameters(CENTERPARS) required by the DAOPHOT tasks.
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.centerpars
    task.unlearn()

    task.calgorithm = 'centroid'                    # Centering Algorithm
    task.cbox = 5                                   # Centering Box Width In Scale Units
    task.cthreshold = 0                             # Centering Threshold In Sigma Above Background


def fitsky_pars(fwhm_value):
    """
    Edits the sky fitting algorithm parameters(FITSKYPARS) requred by the DAOPHOT tasks.
    Args:
        fwhm_value  : Mean FWHM value for the image file
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.fitskypars
    task.unlearn()

    if BackSub:
        task.salgorithm = 'centroid'                        # Sky Fitting Algorithm
        task.annulus = 5 * float(fwhm_value)            # Inner Radius Of Sky Annulus In Scale Units
        task.dannulus = 3                               # Width Of Sky Annulus In Scale Units
    else:
        task.salgorithm = 'constant'                    # Sky Fitting Algorithm
        task.skyvalue = 0                               # User Supplied Sky Value


def phot_pars():
    """
    Edits the photometry parameters(PHOTPARS) required by the DAOPHOT tasks.
    Args:
        None
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.photpars
    task.unlearn()

    task.weighting = 'constant'                             # Photometric Weighting Scheme
    task.aperture = "12,15,18,20,22,25,28,30,32,35,40"   # List Of Aperture Radii In Scale Units
    task.zmag = 25                                          # Zero Point Of Magnitude Scale


def phot(file_name, coord_file):
    """
    Performs PHOT task on the file 'file_name. Selects candidate stars from coordinate file 'coord_file'.
    Args:
        file_name    : FITS file on which aperture photometry is to be performed
        coord_file   : Name of the coordinate file containing candidate star
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.phot
    task.unlearn()

    task.interactive = 'no'                         # Interactive Mode?
    task.radplot = 'no'                             # Plot The Radial Profiles?
    task.verbose = 'no'                             # Print Messages About Progress Of The Task?
    task.verify = 'no'                              # Verify Critical Parameters?
    task.update = 'no'                              # Update Critical Parameters(If Verify Is Yes)?

    task(image=file_name, coords=coord_file, output='default')


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
    group_similar_files(file_temp, common_text=common_text)
    task(textfile='@' + file_temp, fields=fields, expr="yes", Stdout=output_file)
    remove_file(file_temp)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Performing Photometry
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_fwhm(text_list_files, coord_file='stars.coo', log_imexam='log_imexam'):
    """
    Calculates the Mean FWHM of all the files in the list 'list_files'. It determines the FWHM using
    IMEXAMINE task on the stars mentioned in the file "stars.coo".
    Args:
        text_list_files : Text list containing names of FITS files whose FWHM is to be determined
        coord_file      : Text file listing the coordinates of selected stars in the field
        log_imexam      : Name of the text list to record log of IMEXAM
    Returns:
        list_mean_fwhm  : Python list containing Mean FWHM of all the FITS files
    """
    imexam_fwhm(text_list_files, coord_file=coord_file, log_imexam=log_imexam)

    data_file = []
    with open(log_imexam, 'r') as fin:
        for line in fin:
            if not line.startswith('#'):
                data_file += line.rstrip().split()

    columns = 15
    rows = len(data_file) / columns
    no_of_stars = len(text_list_to_python_list(coord_file)) / 2

    list_mean_fwhm = []
    temp_list = []
    for index in range(0, rows):
        temp_list.append(data_file[14 + index * columns])
        if len(temp_list) % no_of_stars == 0:
            mean = float(np.mean(np.array(reject(temp_list))))
            list_mean_fwhm.append(round(mean, 1))
            temp_list = []

    return list_mean_fwhm


def aper_phot(text_list_files, text_list_fwhm, coord_file, data_max='INDEF'):
    """
    Performs aperture photometry (PHOT task) on the files in the list 'list_files'. Selects candidate
    stars from the coordinate file 'coord_file'.
    Args:
        text_list_files : List of all FITS files on which aperture photometry is to be performed
        text_list_fwhm  : List of Mean FWHM values of all the FITS files
        coord_file      : Name of the coordinate file containing candidate star
        data_max        : Maximum good pixel value
    Returns:
        None
    """
    list_files = text_list_to_python_list(text_list_files)
    list_fwhm = text_list_to_python_list(text_list_fwhm)

    for index in range(0, len(list_files)):
        data_pars(list_fwhm[index], data_max)
        center_pars()
        fitsky_pars(list_fwhm[index])
        phot_pars()
        phot(file_name=list_files[index], coord_file=coord_file)

    display_text('Aperture Photometry Is Completed')

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# remove_resfile = eg.boolbox(msg='Remove Residual Files From Previous Run Of This Script?',
#                             title='Remove Residual Files', choices=['Yes', 'No'])
# ctext = eg.enterbox(msg='Enter The Common Text Of Files On Which Photometry Is To Be Done?',
#                     title='Photometry Using IRAF', default='*.fits')
# coord_file = eg.enterbox(msg='Enter The File With Coordinates Of Field Stars:', title='Field Star Coordinates',
#                          default='stars.coo')

remove_resfile = True
ctext = 'HD*.fits'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of Photometry Tasks (PHOT, PSTSELECT, PSF, ALLSTAR)
# ------------------------------------------------------------------------------------------------------------------- #
if remove_resfile:
    for text in ['tmp*', '*.mag.*', 'list_*', 'log_*']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups FITS Files On Which Photometry Is To Be Performed
# ------------------------------------------------------------------------------------------------------------------- #
text_list = 'list_files'
text_list_fwhm = 'list_fwhm'

list_files = group_similar_files(text_list, common_text=ctext)
list_names = set([filename.split('_')[0] for filename in list_files])
dict_coord = {name: name + '.coo' for name in list_names}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Determines The FWHM Of The Stars In The Field (stars.coo)
# Performs Photometry On Images
# ------------------------------------------------------------------------------------------------------------------- #

for name in list_names:
    group_similar_files(text_list + name, common_text=name + '*.fits')
    list_objfwhm = calculate_fwhm(text_list_files=text_list + name, coord_file=dict_coord[name])
    python_list_to_text_list(python_list=list_objfwhm, text_list=text_list_fwhm + name)
    aper_phot(text_list + name, text_list_fwhm + name, dict_coord[name], data_max='INDEF')
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of Photometry Tasks (PHOT, PSTSELECT, PSF, ALLSTAR)
# ------------------------------------------------------------------------------------------------------------------- #
remove_similar_files(common_text='list_')
# ------------------------------------------------------------------------------------------------------------------- #
