#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxx-------------------PHOTOMETRY OF OBJECT FRAMES-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
from pyraf import iraf
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Telescope Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = 'AstroSat_UVIT'
platescale = 0.41625
fwhm_arcsec = 1.31
fwhm_pixels = fwhm_arcsec / platescale
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
EXPTIME_keyword = 'EXP_TIME'
FILTER_keyword = 'FILTER'
DATE_keyword = 'DATE-OBS'
TIME_keyword = 'TIME-OBS'
TSTART_keyword = 'TSTART'
TSTOP_keyword = 'TSTOP'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_NAME = '2018cow'
OBJECT_RA = '16:16:00.22'
OBJECT_DEC = '+22:16:04.83'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Data Containing Information On FILTERS (Extinction Coefficients For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS_UVIT.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Name']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
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
iraf.ccdred.instrument = 'ccddb$kpno/camera.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Handling Files & Lists
# ------------------------------------------------------------------------------------------------------------------- #

def remove_file(file_name):
    """
    Removes the file 'file_name' in the constituent directory.
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
    Removes similar files based on the string 'common_text'.
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

def data_pars(data_max='INDEF'):
    """
    Edits the data dependent parameters(DATAPARS) required by the DAOPHOT tasks.
    Args:
        data_max    : Maximum good pixel value
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.datapars
    task.unlearn()

    task.scale = 1.0  # Scale Of The Image In Arcseconds Per Pixel
    task.fwhmpsf = fwhm_pixels  # FWHM Of The PSF In Scale Units
    task.emission = 'yes'  # All Features Are Considered To Be Emission Features
    task.datamin = 'INDEF'  # Minimum Good Pixel Value
    task.datamax = data_max  # Maximum Good Pixel Value
    task.noise = 'poisson'  # Noise Model Used To Estimate Uncertainties In APPHOT Magnitudes
    task.sigma = 'INDEF'  # Standard Deviation Of The Sky Pixels
    task.readnoise = 0  # Readout Noise Of The CCD In Electrons
    task.epadu = 1  # Gain Of The CCD In Electrons Per ADU
    task.exposure = EXPTIME_keyword  # Exposure Time Keyword In Image Header
    task.airmass = ''  # Airmass Keyword In Image Header
    task.filter = ''  # Filter Keyword In Image Header
    task.obstime = TIME_keyword  # UT Keyword In Image Header


def center_pars():
    """
    Edits the centering algorthm parameters(CENTERPARS) required by the DAOPHOT tasks.
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.centerpars
    task.unlearn()

    task.calgorithm = 'centroid'  # Centering Algorithm
    task.cbox = 5  # Centering Box Width In Scale Units
    task.cthreshold = 0  # Centering Threshold In Sigma Above Background


def fitsky_pars(annulus=100, dannulus=5):
    """
    Edits the sky fitting algorithm parameters(FITSKYPARS) requred by the DAOPHOT tasks.
    Args:
        annulus     : Inner radius of the sky annnulus from the centroid
        dannulus    : Width of the sky annulus
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.fitskypars
    task.unlearn()

    task.unlearn()
    task.salgorithm = 'mode'  # Sky Fitting Algorithm
    task.annulus = int(annulus)  # Inner Radius Of Sky Annulus In Scale Units
    task.dannulus = int(dannulus)  # Width Of Sky Annulus In Scale Units


def phot_pars(aperture_values, zpmag=0):
    """
    Edits the photometry parameters(PHOTPARS) required by the DAOPHOT tasks.
    Args:
        aperture_values : Mean FWHM value for the image file
        zpmag           : Zero Point of the filter in which the observation is performed
    Returns:
        None
    """
    task = iraf.noao.digiphot.daophot.photpars
    task.unlearn()

    task.weighting = 'constant'  # Photometric Weighting Scheme
    task.aperture = aperture_values  # List Of Aperture Radii In Scale Units
    task.zmag = zpmag  # Zero Point Of Magnitude Scale


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

    task.interactive = 'no'  # Interactive Mode?
    task.radplot = 'no'  # Plot The Radial Profiles?
    task.verbose = 'no'  # Print Messages About Progress Of The Task?
    task.verify = 'no'  # Verify Critical Parameters?
    task.update = 'no'  # Update Critical Parameters(If Verify Is Yes)?

    task(image=file_name, coords=coord_file, output='default')


def txdump(file_name, output_file):
    """
    Performs TXDUMP task on the MAG file generated by photometry tasks. This extracts
    useful data from magnitude files.
    Args:
        file_name   : Name of the MAG file from which data is to be extracted
        output_file : Output file where data from the list of input files is to be written
    Returns:
        None
    """
    fields = "RAPERT, SUM, FLUX, MAG, MERR"

    task = iraf.noao.digiphot.ptools.txdump
    task.unlearn()

    task(textfile=file_name, fields=fields, expr='yes', Stdout=output_file)


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Performing Photometry
# ------------------------------------------------------------------------------------------------------------------- #

def extract_val(aper_string):
    """
    Calculates apertures to be calculated in terms of 'Pixels' from a string supplying apertures
    in terms of FWHM value of the image.
    Args:
        aper_string : String specifing apertures in terms of FWHM of the image
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
        aper_values += str(float(value)) + ','

    return aper_values[:-1]


def aper_phot(file_name, coord_file, phot_radius, zpmag=0, annulus=100, dannulus=5, data_max='INDEF'):
    """
    Performs aperture photometry (PHOT task) on the files in the list 'textlist_files'. Selects candidate
    stars from the coordinate file 'coord_file'.
    Args:
        file_name   :   FITS file on which aperture photometry is to be performed
        coord_file  :   Name of the coordinate file containing candidate star
        phot_radius :   String containing the apertures at which photometry is to be done("1,4")
        zpmag       :   Zero point magnitude of the filter in which observation is performed
        annulus     :   Inner radius of the sky annnulus from the centroid
        dannulus    :   Width of the sky annulus
        data_max    :   Maximum good pixel value
    Returns:
        None
    """
    data_pars(data_max)
    center_pars()
    fitsky_pars(annulus, dannulus)
    phot_pars(extract_val(phot_radius), zpmag)

    phot(file_name=file_name + '[0]', coord_file=coord_file)
    display_text("Aperture Photometry Is Completed For Aperture Values (x FWHM): {0}".format(phot_radius))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
remove_resfile = True
aperture_values = '3.2:37.2:0.2'
ctext = '*Sig*.fits'
ctext_coords = 'stars*.coo'

if remove_resfile:
    for text in ['*.mag.*', '*.eps']:
        remove_similar_files(text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Perform Photometry On Images
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files('', common_text=ctext)
list_coord = group_similar_files('', common_text=ctext_coords)

for index, file_name in enumerate(list_files):
    aper_phot(file_name, list_coord[index], phot_radius=aperture_values)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Light Curve For One Set of Observation
# ------------------------------------------------------------------------------------------------------------------- #
name_cols = ['APER', 'SUM', 'AREA', 'FLUX']
output_cols = ['PHASE', 'FILTER', 'EXPTIME', 'COUNTS', 'MAG', 'ERR']

list_magwithgal = group_similar_files('', common_text='*.mag.1')

dict_magwithgal = {}
for file_name in list_magwithgal:
    header = fits.getheader(file_name[:-7], ext=0)
    exptime = float(header[EXPTIME_keyword])
    band = header[FILTER_keyword]
    zpmag = filter_df.loc[band, 'ZeroPoint']
    zperr = filter_df.loc[band, 'ZPErr']

    data_uv = pd.read_csv(file_name, sep='\s+', usecols=[0, 1, 2, 3], names=name_cols, skiprows=79)
    data_uv = data_uv.replace('INDEF', np.nan).astype('float64')
    counts = data_uv.loc[data_uv['APER'] == 18., 'FLUX'].values[0]

    obsepoch = header[DATE_keyword] + 'T' + header[TIME_keyword][:-6]
    dict_magwithgal[obsepoch] = {}
    dict_magwithgal[obsepoch]['FILTER'] = band + '-' + filter_df.loc[band, 'Name']
    dict_magwithgal[obsepoch]['PHASE'] = round((header[TSTART_keyword] + header[TSTOP_keyword]) / 2., 3)
    dict_magwithgal[obsepoch]['EXPTIME'] = round(exptime, 3)
    dict_magwithgal[obsepoch]['COUNTS'] = round(counts, 7)
    dict_magwithgal[obsepoch]['MAG'] = round(-2.5 * np.log10(counts) + zpmag, 3)
    dict_magwithgal[obsepoch]['ERR'] = round(((1.087 / ((counts * exptime) ** 0.5)) ** 2 + zperr ** 2) ** 0.5, 3)

data_magwithgal = pd.DataFrame(dict_magwithgal).T
data_magwithgal.index.name = 'OBSTIME'
data_magwithgal = data_magwithgal[output_cols]
data_magwithgal.to_csv('OUTPUT_MagWithGal', sep=' ', index=True, header=True)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The LC From AstroSat UVIT
# ------------------------------------------------------------------------------------------------------------------- #

for band, band_df in data_magwithgal.groupby('FILTER'):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)

    ax.errorbar((band_df['PHASE'] - band_df['PHASE'].min()) / 3600, band_df['MAG'], yerr=band_df['ERR'], marker='*',
                color='k', markersize=12, linestyle='--', capthick=1, capsize=5, elinewidth=1)

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.tick_params(which='minor', direction='in', length=3, width=1, labelsize=14)
    ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=14)

    ax.invert_yaxis()
    ax.set_xlabel('Time Since First Observation [In Hours]', fontsize=16)
    ax.set_ylabel('FUV Magnitude [mag]', fontsize=16)

    fig.savefig('PLOT_' + band + '.eps', format='eps', dpi=1000, bbox_inches='tight')
    # plt.show()
    plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
