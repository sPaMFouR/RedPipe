#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx-------------------------FLUX CALIBRATION OF 1-D Spectra----------------------xxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import shutil
import numpy as np
import pandas as pd
import easygui as eg
from pyraf import iraf
from astropy.io import fits
import matplotlib.pyplot as plt
import specutils.io.read_fits as spec
from scipy.interpolate import CubicSpline, Rbf
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
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
JD_keyword = 'JD'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Files & Directories To Be Used
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_HOME = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
DIR_PHOT = "/home/avinash/Supernovae_Data/2016gfy/Photometry/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
FILE_BANDPASS = os.path.join(DIR_SPEC, 'Filter_BPF.asc')
list_paths = [DIR_HOME, DIR_PHOT, DIR_SPEC, FILE_BANDPASS]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
iraf.twodspec(_doprint=0)
iraf.onedspec(_doprint=0)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_HOME + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.set_index('FILTER')
list_filters = filter_df.index.tolist()

for index, row in filter_df.iterrows():
    if row['Offset'] > 0:
        filter_df.loc[index, 'Label'] = index + ' + ' + str(row['Offset'])
    elif row['Offset'] == 0:
        filter_df.loc[index, 'Label'] = index
    else:
        filter_df.loc[index, 'Label'] = index + ' - ' + str(abs(row['Offset']))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
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


def copy_files(inp_path, out_path, common_text, exceptions=''):
    """
    Copies similar files based on the string 'common_text' from the directory specified by 'inp_path'
    onto the directory specified by 'out_path'.
    Args:
        inp_path    : Path of the directory from which files are to be copied
        out_path    : Path of the directory to which files are to be copied
        common_text : String containing partial name of the files to be copied
        exceptions  : String containing the partial name of the files that need to be excluded
    Returns:
        None
    """
    os.chdir(inp_path)

    list_copy = group_similar_files('', common_text=common_text, exceptions=exceptions)
    for file_name in list_copy:
        shutil.copy(os.path.join(inp_path, file_name), out_path)

    os.chdir(DIR_CURNT)


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
        print ("\nError : File '{0}' Not Found\n".format(text_list))
        sys.exit(1)


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


def check_ifexists(path):
    """
    Checks if a file or directory exists.
    Args:
        path        : Path whose existence is to be checked
    Returns:
        True        : Returns True only when the path exists
    """
    if path[-1] == '/':
        if os.path.exists(path):
            pass
        else:
            print ("\nError : Directory '{0}' Does Not Exist\n".format(path))

    else:
        if os.path.isfile(str(path)):
            pass
        else:
            print ("\nError : File '{0}' Does Not Exist\n".format(path))


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

def sbands(common_text, output_file='OUTPUT_sbands'):
    """
    Performs bandpass spectrophotometry of 1-D spectra.
    Args:
        common_text : Common text of 1-D spectra files whose flux in different bands is to be calculated
        output_file : Name of the output file to record spectroscopic fluxes
    Returns:
        None
    """
    textlist_files = 'list_smspec'
    group_similar_files(textlist_files, common_text=common_text)

    task = iraf.noao.onedspec.sbands
    task.unlearn()

    task.normalize = 'yes'              # Normalize The Bandpasss Response?
    task.mag = 'no'                     # Output Results In Magnitudes?
    task.verbose = 'no'                 # Verbose Header?
    task.magzero = '0'                  # Magnitude Zero Point

    if os.path.isfile(FILE_BANDPASS):
        remove_file(output_file)
        task(input='@' + textlist_files, output=output_file, bands=FILE_BANDPASS)
        remove_file(textlist_files)
    else:
        print ("\nError : '{0}' Does Not Exist\n".format(FILE_BANDPASS))


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Reading & Writing Data From FITS Files
# ------------------------------------------------------------------------------------------------------------------- #

def read_jd(file_name):
    """
    Reads JD of observation of the file "file_name".
    Args:
        file_name   : Name of the 1-D Spectra whose JD of observation is to be found out
    Returns:
        julian_day  : Julian day of the 1-D spectra
    """
    julian_day = fits.getval(filename=file_name, keyword=JD_keyword)

    return julian_day


def read_1dspec(file_name):
    """
    Reads 1-D Spectra from a FITS file and returns wavelength and flux arrays.
    Args:
        file_name    : FITS file from which data has to be extracted
    Returns:
        wave_array   : Array containing wavelength values extracted from the 1-D Spectra
        flux_array   : Array containing flux values extracted from the 1-D Spectra
    """
    with fits.open(file_name) as hdulist:
        axis = int(hdulist[0].header['NAXIS'])
        if axis == 1:
            flux_array = hdulist[0].data
            wave_array = spec.read_fits_spectrum1d(file_name).dispersion
        else:
            flux_array = hdulist[0].data[0][0]
            wave_array = spec.read_fits_spectrum1d(file_name)[0].dispersion

    return wave_array, flux_array


def write_1dspec(ref_filename, flux_array, prefix_str):
    """
    Writes 1-D Spectra onto a FITS file.
    Args:
        ref_filename : FITS file from which header has to be extracted
        flux_array   : Array containing flux values
        prefix_str   : Prefix to distinguish the smoothened 1-D spectra from the original
    Returns:
        None
    """
    with fits.open(ref_filename) as hdulist:
        file_header = hdulist[0].header

    output_filename = prefix_str + ref_filename
    remove_file(output_filename)
    fits.writeto(output_filename, data=flux_array, header=file_header)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Accessing & Manipulating Text File Data
# ------------------------------------------------------------------------------------------------------------------- #

def read_specflux(file_name, file_specflux='OUTPUT_sbands'):
    """
    Reads spectroscopic fluxes from a text file 'file_specflux'. The fluxes are determined
    for the spectra specified by the file 'file_name'.
    Args:
        file_name     : 1-D Spectra whose spectroscopic fluxes are to be extracted
        file_specflux : Text file from which spectroscopic fluxes are to be extracted
    Returns:
        dict_specflux : Dictionary of spectroscopic fluxes in different bands
    """
    sbands_columns = ['FileName', 'FILTER', 'Flux']
    flux_df = pd.read_csv(file_specflux, header=None, names=sbands_columns, sep='\s+', engine='python')
    flux_df['Flux'] = flux_df['Flux'].apply(lambda x: '{0:6.4e}'.format(x))

    dict_specflux = {}
    for index, row in flux_df.iterrows():
        if re.search(file_name, row['FileName']):
            dict_specflux[filter_df.loc[row['FILTER'], 'CentreWave']] = row['Flux']

    return dict_specflux


def read_photflux(list_photfiles, julian_day, flux=True):
    """
    Reads broadband photometric magnitudes from a text list containing names of files containing
    photometric magnitudes. The magnitudes are determined for the epoch specified by 'julian_day'.
    Args:
        list_photfiles  : Text list of files containing different broadband photometric magnitudes
        julian_day      : Julian day close to which the photometric magnitude has to be extracted
        flux            : Boolean describing whether flux has to be returned (or magnitude)
    Returns:
        dict_photmag    : Dictionary of photometric magnitudes in different bands
        dict_photflux   : Dictionary of photometric fluxes in different bands
    """
    list_photfiles = text_list_to_python_list(list_photfiles)

    dict_photmag = {}
    dict_photflux = {}
    for file_photmag in list_photfiles:
        file_df = pd.read_csv(file_photmag, sep='\s+', engine='python').astype('float64')
        file_df = file_df[abs(file_df['JD'] - float(julian_day)) <= 0.25]

        if file_df.shape[0] == 1:
            dict_photmag[filter_df.loc[file_photmag[-1], 'CentreWave']] = file_df['Mag'].values[0]
            dict_photflux[filter_df.loc[file_photmag[-1], 'CentreWave']] = mag_to_flux(file_df['Mag'].values[0],
                                                                                       file_band=file_photmag)
        else:
            dict_photmag[filter_df.loc[file_photmag[-1], 'CentreWave']] = file_df['Mag'].mean()
            dict_photflux[filter_df.loc[file_photmag[-1], 'CentreWave']] = mag_to_flux(file_df['Mag'].mean(),
                                                                                       file_band=file_photmag)
    if flux:
        return dict_photflux
    else:
        return dict_photmag


def mag_to_flux(mag, file_band):
    """
    Converts magnitudes to flux values (Not extinction corrected).
    Args:
        mag         : Magnitude value to be converted to flux
        file_band   : Text list of files containing broadband photometric magnitudes
    Returns:
        flux        : Flux value corresponding to the input magnitude
    """
    if file_band[-1] in filter_df.index.values.tolist():
        flux = '{0:6.4e}'.format(10 ** (-0.4 * (mag + float(filter_df.loc[file_band[-1], 'ZeroPoint']) + 21.100)))
    else:
        print ("Error: Band Of Observation '{0}' Not Recognised".format(file_band[-1]))
        sys.exit(1)

    return flux


def get_zflux(dict_phot, cntrl_wav=7500):
    """
    Obtains Z-band (narrow band, 7500 Angstroms) flux value.
    Args:
        dict_phot   : Dictionary containing broadband photometric flux values
        cntrl_wav   : Central wavelength of the Z-band
    Returns:
        dict_phot   : Modified dictionary with Z-band flux value included
    """
    data_series = pd.Series(dict_phot, dtype='float64').dropna()
    spline = CubicSpline(data_series.index.values, data_series.values, bc_type='natural', extrapolate=True)
    dict_phot[int(cntrl_wav)] = '{0:6.4e}'.format(float(spline(int(cntrl_wav))))

    return dict_phot

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Smoothening Of 1-D Spectra
# ------------------------------------------------------------------------------------------------------------------- #

def smooth_1dspec(common_text, sp=1.2, kernel='gaussian', prefix_str='z_', plot=False):
    """
    Smoothens a 1-D spectra based on the smoothening parameter. Smoothening parameter
    is 'std.dev.' in case of isotropic Gaussian filter and is 'width' in the case of the
    non-isotropic box filter.
    Args:
        common_text : Common text of 1-D spectra files which have to be smoothened
        sp          : Smoothening parameter
        kernel      : Convolution Kernel used for smoothening (Gaussian or Box)
        prefix_str  : Prefix to distinguish the smoothened 1-D spectra from the original
        plot        : Boolean describing whether the smoothened spectra has to be plotted
    Returns:
        None
    """
    list_files = group_similar_files('', common_text=common_text)

    for file_name in list_files:
        wav_data, flux_data = read_1dspec(file_name)
        usable_kernel = Gaussian1DKernel(int(sp))

        if kernel.lower() != 'gaussian':
            if kernel.lower() == 'box':
                usable_kernel = Box1DKernel(int(sp))
            else:
                print ("Error: Kernel '{0}' Not Recognised".format(kernel))
                sys.exit(1)

        smoothed_data = convolve(flux_data, usable_kernel)
        write_1dspec(ref_filename=file_name, flux_array=smoothed_data, prefix_str=prefix_str)

        if plot:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            ax.plot(wav_data, flux_data, 'g', label='Original Spectrum')
            ax.plot(wav_data, smoothed_data, 'r', label='Smooth Spectrum')

            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Scaling Spectra With The Help Of Photometric Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

def scale_spectra(common_text, list_photfiles, prefix_str='f', plot=False):
    """
    Scales spectra acccording to the values in the array 'scale_array'. Basically, this step applies
    flux calibration on the spectra.
    Args:
        common_text     : Common text of 1-D Spectra files whose spectroscopic fluxes are to be scaled
        list_photfiles  : Text list of files containing different broadband photometric magnitudes
        prefix_str      : Prefix to distinguish the scaled 1-D spectra from the original
        plot            : Boolean describing whether the scaled spectra has to be plotted
    Returns:
        None
    """
    list_files = group_similar_files('', common_text=common_text)

    for file_name in list_files:
        dict_spec = read_specflux(file_name)
        dict_phot = read_photflux(list_photfiles=list_photfiles, julian_day=read_jd(file_name))
        dict_phot = get_zflux(dict_phot, cntrl_wav=7500)
        dict_scale = dict((key, str(float(dict_phot[key]) / float(dict_spec[key]))) for key in dict_phot.keys() if
                          key in dict_spec.keys())

#         if len(dict_scale.keys()) > 3:
#             del dict_scale[7500]

        # if len(dict_scale.keys()) > 4:
        #     order = 4
        # else:
        #     order = len(dict_scale.keys()) - 1

        series = pd.Series(dict_scale, dtype='float64').dropna()
        spline = CubicSpline(series.index.values, series.values, bc_type='natural', extrapolate=True)
        spline2 = Rbf(series.index.values, series.values)

        wave_data, flux_data = read_1dspec(file_name)
        scale_data = spline(wave_data)
        scale_data[scale_data < 0] = 0

        flux_moddata = np.multiply(np.asarray(flux_data), scale_data)
        write_1dspec(ref_filename=file_name, flux_array=flux_moddata, prefix_str=prefix_str)

        if plot:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111)

            wavenew = np.linspace(float(wave_data[0]), float(wave_data[-1]), 10000)
            ax.plot(series.index.values, series.values, 'o', label='Data Points')
            ax.plot(wavenew, spline(wavenew), 'k', label='CubicSpline')
            ax.plot(wavenew, spline2(wavenew), 'r', label='Rbf')

            ax.legend()
            ax.grid()
            plt.show()
            plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Checks If Directories To Be Used Exist Or Not
# ------------------------------------------------------------------------------------------------------------------- #
for path in list_paths:
    check_ifexists(path)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copies Files Needed For Flux Calibration (sp_*.asc) And Supernova Magnitudes To The Current Working Directory
# ------------------------------------------------------------------------------------------------------------------- #
# copy_files(inp_path=DIR_HOME, out_path=DIR_SPEC, common_text='Filter_*.asc', exceptions='BPF,SDSS')
copy_files(inp_path=DIR_PHOT, out_path=DIR_SPEC, common_text='OUTPUT_InterpSNMag*', exceptions='Row,Col')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files',
#                        choices=['Yes', 'No'])
# ctext = eg.enterbox('Enter The Common Text Of Files To Be Flux Calibrated?', title='Flux Calibration',
#                     default='*cfwcbs_*.ms.fits')
# bool_smooth = eg.boolbox('Perform Smoothening Of Spectra?', title='Smoothening 1-D Spectra', choices=['Yes', 'No'])

rmv_files = True
ctext = '*.fits'
bool_smooth = True
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes Residual Files From Previous Run Of Flux Calibration
# ------------------------------------------------------------------------------------------------------------------- #
os.chdir(DIR_SPEC)
if rmv_files:
    for text in ['z_*.fits', 'fz_*.fits', 'rfz_*.fits', 'list_smspec']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Smoothening, SBANDS Task & Finally Flux Calibration On 1-D Spectra
# ------------------------------------------------------------------------------------------------------------------- #
if bool_smooth:
    smooth_1dspec(common_text=ctext, sp=2, kernel='gaussian')
    ctext = 'z_' + ctext
    display_text('Smoothening of 1-D Spectra Have Been Performed')

sbands(common_text=ctext)
group_similar_files("list_interpmag", common_text='OUTPUT_InterpSNMag*')

scale_spectra(common_text=ctext, list_photfiles="list_interpmag", plot=False)
display_text('Flux Calibration Has Been Performed On All Of 1-D Spectra')
# ------------------------------------------------------------------------------------------------------------------- #
