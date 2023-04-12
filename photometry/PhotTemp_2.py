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
import pandas as pd
import easygui as eg
from pyraf import iraf
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = 'Indian Astronomical Observatory, Hanle'
OBS_LONG = '78:57:51'
OBS_LAT = '32:46:46'
OBS_ALT = 4486
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
# Extinction Coefficients (In Magnitudes) For Hanle In Different Photometric  Bands
# ------------------------------------------------------------------------------------------------------------------- #
eeta = {'7BesU': 0.36, '6BesB': 0.21, '5BesV': 0.12, '4BesR': 0.09, '3BesI': 0.05}
eeta_err = {'7BesU': 0.07, '6BesB': 0.04, '5BesV': 0.04, '4BesR': 0.04, '3BesI': 0.03}
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
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_NAME = '2022jli'
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
iraf.ccdred.instrument = "ccddb$kpno/camera.dat"
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
    """s
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
        Error : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("Error : File '{0}' Not Found".format(text_list))
        sys.exit(1)


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


def zip_list(list_lists):
    """
    Combines all the lists in a list to a single list element-wise.
    Args:
        list_lists  : List of all the lists which needs to be zipped
    Returns:
        new_list    : List with the combined elements
    """
    new_list = []
    for index in range(0, len(list_lists[0])):
        value = ''
        for val in range(0, len(list_lists)):
            value += list_lists[val][index]
        new_list.append(str(value))

    return new_list


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
        task.logfile = log_imexam	            # Log File To Record Output Of The Commands
        task.keeplog = 'yes'                        # Log Output Results?
        task.defkey = 'a'                           # Default Key For Cursor x-y Input List
        task.imagecur = coord_file	            # Image Display Cursor Input
        task.use_display = 'no'                     # Use The Image Display?

        display_text("CURRENT FILE - " + file_name)
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
    task.epadu = ccd_gain                           # Gain Of The CCD In Electrons Per ADU
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

    task.unlearn()
    task.salgorithm = 'mode'                        # Sky Fitting Algorithm
    task.annulus = 5 * float(fwhm_value)            # Inner Radius Of Sky Annulus In Scale Units
    task.dannulus = 3                               # Width Of Sky Annulus In Scale Units


def phot_pars(aperture_values):
    """
    Edits the photometry parameters(PHOTPARS) required by the DAOPHOT tasks.
    Args:
        aperture_values : Mean FWHM value for the image file
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.photpars
    task.unlearn()

    task.weighting = 'constant'                     # Photometric Weighting Scheme
    task.aperture = aperture_values                 # List Of Aperture Radii In Scale Units
    task.zmag = 25                                  # Zero Point Of Magnitude Scale


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
    task(textfile='@' + file_temp, fields=fields, expr='yes', Stdout=output_file)
    remove_file(file_temp)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Accessing & Manipulating Text File Data
# ------------------------------------------------------------------------------------------------------------------- #

def read_column(file_name, col_index, title_rows=0):
    """
    Extracts the specified column as a list from a text file.
    Args:
        file_name   : Text file from which the specified column has to be extracted
        col_index   : Index of the column to be extracted
        title_rows  : No. of rows used for title description
    Returns:
        data_column : List of all the elements extracted from the column
    """
    file_df = pd.read_csv(filepath_or_buffer=file_name, sep='\s+', header=None, skiprows=title_rows)
    data_column = file_df.iloc[:, col_index].tolist()

    return data_column


def read_file(file_name, title_rows=0):
    """
    Extracts the file data as a list of columns from a text file.
    Args:
        file_name   : Text file from which file data has to be extracted
        title_rows  : No. of rows used for title description
    Returns:
        data_file   : List of all columns extracted from the text file
    """
    file_df = pd.read_csv(filepath_or_buffer=file_name, sep='\s+', header=None, skiprows=title_rows)
    data_file = [file_df.iloc[:, index].tolist() for index in range(0, file_df.shape[1])]

    return data_file


def read_magfile(file_name, col_nos, fmt='{:>8}', title_rows=0):
    """
    Extracts the mag file data as a list of columns specified by 'col_nos' from a text file and formats
    the list according to the format specified in 'fmt'.
    Args:
        file_name   : Text file from which file data has to be extracted
        col_nos     : Indexes of columns to be read from the file ('7:9' or '7,8,9')
        fmt         : Format of storing data in the list
        title_rows  : No. of rows used for title description
    Returns:
        data_file   : List of all columns extracted from the text file
    """
    file_df = pd.read_csv(filepath_or_buffer=file_name, sep='\s+', header=None, skiprows=title_rows)
    rows, columns = file_df.shape

    if re.search(':', col_nos):
        col_indexes = range(int(col_nos.split(':')[0]), int(col_nos.split(':')[-1]), 1)
    elif re.search(',', col_nos):
        col_indexes = col_nos.split(',')
    else:
        print ("Invalid Format Of Entering Column Numbers")
        sys.exit(1)
        
    data_file = [file_df.iloc[:, index].tolist() for index in col_indexes]

    for col_index in range(0, len(col_indexes)):
        for row_index in range(0, rows):
            try:
                float(data_file[col_index][row_index])
            except ValueError:
                new_fmt = fmt[0:3] + 's}'
                data_file[col_index][row_index] = new_fmt.format(str(data_file[col_index][row_index]))
            else:
                data_file[col_index][row_index] = fmt.format(float(data_file[col_index][row_index]))

    return data_file

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Performing Photometry
# ------------------------------------------------------------------------------------------------------------------- #

def calculate_fwhm(textlist_files, coord_file='stars.coo', log_imexam='log_imexam'):
    """
    Calculates the Mean FWHM of all the files in the list 'list_files'. It determines the FWHM using
    IMEXAMINE task on the stars mentioned in the file "stars.coo".
    Args:
        textlist_files  : Text list containing names of FITS files whose FWHM is to be determined
        coord_file      : Text file listing the coordinates of selected stars in the field
        log_imexam      : Name of the text list to record log of IMEXAM
    Returns:
        list_mean_fwhm  : Python list containing Mean FWHM of all the FITS files
    """
    imexam_fwhm(textlist_files, coord_file=coord_file, log_imexam=log_imexam)
    coord_df = pd.read_csv(coord_file, sep='\s+', header=None)
    data_df = pd.read_csv(log_imexam, sep='\s+', comment="#", header=None)
    count = coord_df.shape[0]
    rows, columns = data_df.shape
    col_moff = columns - 2

    list_fwhm = [data_df.iloc[0 + count * i: (i + 1) * count, col_moff].tolist() for i in range(0, rows / count)]

    list_mean_fwhm = []
    for index, fwhm_values in enumerate(list_fwhm):
        mean = float(np.mean(a=reject(fwhm_values)))
        list_mean_fwhm.append(round(mean, 1))

    return list_mean_fwhm


def extract_val(aper_string, fwhm):
    """
    Calculates apertures to be calculated in terms of 'Pixels' from a string supplying apertures
    in terms of FWHM value of the image.
    Args:
        aper_string : String specifing apertures in terms of FWHM of the image
        fwhm        : FWHM of the image to which photometry is being done
    Returns:
        aper_values : String containing apertures to be used for photometry
    """
    if re.search(':', aper_string):
        list_aper = aper_string.split(':')
        if len(list_aper) == 2:
            list_aper = np.arange(float(list_aper[0]), 1 + float(list_aper[1]), 1)
        elif len(aper_string.split(':')) == 3:
            list_aper = np.arange(float(list_aper[0]), float(list_aper[2]) + float(list_aper[1]), float(list_aper[2]))
    else:
        list_aper = aper_string.split(',')

    aper_values = ''
    for value in list_aper:
        aper_values += str(float(value) * float(fwhm)) + ','

    return aper_values[:-1]


def aper_phot(textlist_files, textlist_fwhm, coord_file, phot_radius='1', data_max='INDEF'):
    """
    Performs aperture photometry (PHOT task) on the files in the list 'list_files'. Selects candidate
    stars from the coordinate file 'coord_file'.
    Args:
        textlist_files  : List of all FITS files on which aperture photometry is to be performed
        textlist_fwhm   : List of Mean FWHM values of all the FITS files
        coord_file      : Name of the coordinate file containing candidate star
        phot_radius     : String containing the apertures at which photometry is to be done("1,4")
        data_max        : Maximum good pixel value
    Returns:
        None
    """
    list_files = text_list_to_python_list(textlist_files)
    list_fwhm = text_list_to_python_list(textlist_fwhm)

    for index in range(0, len(list_files)):
        aperture_values = extract_val(str(phot_radius), list_fwhm[index])

        data_pars(list_fwhm[index], data_max)
        center_pars()
        fitsky_pars(list_fwhm[index])
        phot_pars(aperture_values)
        phot(file_name=list_files[index], coord_file=coord_file)

    display_text('Aperture Photometry Is Completed For Aperture Values (x FWHM): ' + str(phot_radius))


def tabular_mag(input_file, output_file):
    """
    Takes the output from a MAG or ALS file and computes a tabular magnitude file.
    Args:
        input_file  : Input MAG or ALS file
        output_file : Name of the output file onto which the tabular magnitudes are to be written
    Returns:
        None
    """
    data_file = read_file(input_file)
    [star_id, _, ifilter, xcenter, ycenter, sky_counts, airmass] = data_file[0:7]
    list_images = ['ca_' + input_file.split('_')[1] + '_fbs_' + OBJECT_NAME + '-' + band[-1].lower() + '.fits'
                   for band in ifilter]

    col_data = 7
    star_count = max(map(int, star_id))
    columns = len(data_file)
    rows = len(data_file[0])
    apertures = (columns - col_data) / 3

    list_col = []
    for index in range(0, 3):
        list_col.append(str(col_data + index * apertures) + ':' + str(col_data + (index + 1) * apertures))

    list_aper = zip_list(read_magfile(str(input_file), col_nos=list_col[0], fmt="{:8.2f}"))
    list_mag = zip_list(read_magfile(str(input_file), col_nos=list_col[1], fmt="{:10.4f}"))
    list_err = zip_list(read_magfile(str(input_file), col_nos=list_col[2], fmt="{:8.3f}"))

    aper_names = ""
    mag_names = ""
    err_names = ""
    for value in range(1, apertures + 1):
        aper_names += "{:8s}".format("AP_" + str(value))
        mag_names += "{:10s}".format("MAG_" + str(value))
        err_names += "{:8s}".format("ERR_" + str(value))

    with open(output_file, 'w') as fout:
        fout.write("{0:>4s}{1:>10s}{2:>11s}{3:>13s}  {4}{5}{6}\n\n".format
                   ("ID", "XCENTER", "YCENTER", "SKY_COUNTS", aper_names, mag_names, err_names))

        for index in range(0, rows):
            if index % star_count == 0:
                fout.write("# IMAGE #{0} - {1}, FILTER - {2}, AIRMASS - {3}\n".format
                           ((index / star_count) + 1, list_images[index], ifilter[index], airmass[index]))

            fout.write("{0:>3.0f}{1:>11.3f}{2:>11.3f}{3:>11.4f}{4}{5} {6}\n".format
                       (float(star_id[index]), float(xcenter[index]), float(ycenter[index]), float(sky_counts[index]),
                        list_aper[index], list_mag[index], list_err[index]))

            if index % star_count == star_count - 1:
                fout.write("\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# remove_resfile = eg.boolbox(msg='Remove Residual Files From Previous Run Of This Script?',
#                             title='Remove Residual Files', choices=['Yes', 'No'])
# ctext = eg.enterbox(msg='Enter The Common Text Of Raw Photometric Files:', title='Raw Files', default='*.fits')
# ctext_temp = eg.enterbox(msg='Enter The Common Text Of Template Subtracted Files:', title='Template Subtracted Files',
#                          default='ts_*.fits')
# coord_stars = eg.enterbox(msg='Enter The File With Coordinates Of Field Stars:', title='Field Star Coordinates',
#                           default='stars.coo')
# coord_sn = eg.enterbox(msg='Enter The File With Coordinates Of Supernova:', title='Supernova Coordinates',
#                        default='sn.coo')
# aperture_values = eg.enterbox(msg='Enter The Apertures At Which Photometry Is To Be Performed:',
#                               title='Apertures For Performing Photometry', default='1')

remove_resfile = True
ctext = '*.fits'
ctext_temp = 'ts_*.fits'
coord_stars = 'stars.coo'
coord_sn = 'sn.coo'
aperture_values = '1,4'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Remove Residual Files From Previous Run Of Photometry Tasks (PHOT)
# ------------------------------------------------------------------------------------------------------------------- #
if remove_resfile:
    for text in ['tmp*', 's_*', 'conv_*', 'sc*', 'list_*', 'log*', 'ts_*.mag.*']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Determines The FWHM Of The Image Using Stars Mentioned In The File "stars.coo"
# ------------------------------------------------------------------------------------------------------------------- #
textlist_raw = 'list_files'
textlist_fwhm = 'list_fwhm'

group_similar_files(textlist_raw, common_text=ctext, exceptions='ts_,template,psf,sub')
list_fwhm = calculate_fwhm(textlist_files=textlist_raw, coord_file=coord_stars)
python_list_to_text_list(python_list=list_fwhm, text_list=textlist_fwhm)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups FITS Files On Which Photometry Is To Be Performed & Performs Photometry On SN
# ------------------------------------------------------------------------------------------------------------------- #
textlist_temp = 'list_temp'
group_similar_files(textlist_temp, common_text=ctext_temp)
aper_phot(textlist_temp, textlist_fwhm, coord_sn, phot_radius=aperture_values, data_max='INDEF')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Mag Files From PHOT Task Into A Separate List & Makes A List Of Dates On Which Observation Was Done
# ------------------------------------------------------------------------------------------------------------------- #
mag_suffix = 1
list_mag = group_similar_files('', common_text='ts_*.mag.1')
list_dates = set([file_name[6:16] for file_name in list_mag])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups MAG Files From A Date Into A File 'list_date_mag1', Applies TXDUMP Task To Obtain 'output_date_mag1'
# Also Computes Tabular Magnitude Files From IRAF MAG Files Generated Through Photometry
# ------------------------------------------------------------------------------------------------------------------- #
for date in list_dates:
    txdump(common_text='ts_*' + date + '*.mag.1', output_file='output_' + date + '_mag1')
    tabular_mag(input_file='output_' + date + '_mag1', output_file='OUTPUT_tabulartemp_' + date)

display_text("Tabular Magnitudes Have Been Computed For MAG Files")
# ------------------------------------------------------------------------------------------------------------------- #

