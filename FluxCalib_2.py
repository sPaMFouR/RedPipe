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
import easygui as eg
from pyraf import iraf
from astropy.io import fits
import matplotlib.pyplot as plt
import specutils.io.read_fits as spec
from scipy.interpolate import UnivariateSpline
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
# Bessell Filter Central Wavelengths
# ------------------------------------------------------------------------------------------------------------------- #
dict_centralwav = {'U': 3700, 'B': 4200, 'V': 5300, 'R': 6000, 'I': 8050, 'Z': 7000, 'Y': 8500}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Zero Point Correction For Different Photometric Bands
# ------------------------------------------------------------------------------------------------------------------- #
zp_u = -0.152
zp_b = -0.602
zp_v = 0.000
zp_r = 0.555
zp_i = 1.271
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
JD_keyword = 'JD'
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
# Paths Of Files & Directories To Be Used
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_PHOT = "/home/avinash/Supernovae_Data/Photometry/"
DIR_SPECS = "/home/avinash/Supernovae_Data/Final_Spectra/"
FILE_BANDPASS = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/Filter_BPF.asc"
list_paths = [DIR_PHOT, DIR_SPECS, FILE_BANDPASS]
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


def copy_files(in_path, out_path, common_text, exceptions=''):
    """
    Copies similar files based on the string "common_text" from the directory specified by "in_path"
    onto the directory specified by "out_path".
    Args:
        in_path     : Path of the directory from which files are to be copied
        out_path    : Path of the directory to which files are to be copied
        common_text : String containing partial name of the files to be copied
        exceptions  : String containing the partial name of the files that need to be excluded
    Returns:
        None
    """
    owd = os.getcwd()
    os.chdir(str(in_path))

    list_copy = group_similar_files("", common_text=str(common_text), exceptions=str(exceptions))
    for file_name in list_copy:
        shutil.copy(os.path.join(str(in_path), str(file_name)), str(out_path))

    os.chdir(owd)


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
            pass
        else:
            print ("\nError : Directory '{0}' Does Not Exist\n".format(str(path)))

    else:
        if os.path.isfile(str(path)):
            pass
        else:
            print ("\nError : File '{0}' Does Not Exist\n".format(str(path)))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def sbands(ctext, output_file="output_sbands"):
    """
    Performs bandpass spectrophotometry of 1-D spectra.
    Args:
        ctext       : Common text of 1-D spectra files whose flux in different bands is to be calculated
        output_file : Name of the output file to record spectroscopic fluxes
    Returns:
        None
    """
    text_list_files = "list_smspec"
    group_similar_files(str(text_list_files), common_text=str(ctext))

    task = iraf.noao.onedspec.sbands
    task.unlearn()

    task.normalize = 'yes'                              # Normalize The Bandpasss Response?
    task.mag = 'no'                                     # Output Results In Magnitudes?
    task.verbose = 'no'                                 # Verbose Header?
    task.magzero = '0'                                  # Magnitude Zero Point

    if os.path.isfile(FILE_BANDPASS):
        remove_file(str(output_file))
        task(input='@' + str(text_list_files), output=str(output_file), bands=str(FILE_BANDPASS))
        remove_file(str(text_list_files))
    else:
        print ("\nError : '{0}' Does Not Exist\n".format(str(FILE_BANDPASS)))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Manipulating 1-D Spectra (Read, Write, Smoothen etc.)
# ------------------------------------------------------------------------------------------------------------------- #

def read_1dspec(file_name):
    """
    Reads 1-D Spectra from a FITS file and returns wavelength and flux arrays.
    Args:
        file_name    : FITS file from which data has to be extracted
    Returns:
        wave_array   : Array containing wavelength values extracted from the 1-D Spectra
        flux_array   : Array containing flux values extracted from the 1-D Spectra
    """
    with fits.open(str(file_name)) as hdulist:
        axis = int(hdulist[0].header['NAXIS'])
        if axis == 1:
            flux_array = hdulist[0].data
            wave_array = spec.read_fits_spectrum1d(str(file_name)).dispersion
        else:
            flux_array = hdulist[0].data[0][0]
            wave_array = spec.read_fits_spectrum1d(str(file_name))[0].dispersion

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
    with fits.open(str(ref_filename)) as hdulist:
        file_header = hdulist[0].header

    output_filename = str(prefix_str) + str(ref_filename)
    remove_file(str(output_filename))
    fits.writeto(str(output_filename), data=flux_array, header=file_header)


def smooth_1dspec(ctext, sp=1.2, kernel="gaussian", prefix_str='z_', plot=False):
    """
    Smoothens a 1-D spectra based on the smoothening parameter. Smoothening parameter
    is 'std.dev.' in case of isotropic Gaussian filter and is 'width' in the case of the
    non-isotropic box filter.
    Args:
        ctext       : Common text of 1-D spectra files which have to be smoothened
        sp          : Smoothening parameter
        kernel      : Convolution Kernel used for smoothening (Gaussian or Box)
        prefix_str  : Prefix to distinguish the smoothened 1-D spectra from the original
        plot        : Boolean describing whether the smoothened spectra has to be plotted
    Returns:
        None
    """
    list_files = group_similar_files("", common_text=str(ctext))

    for file_name in list_files:
        wav_data, flux_data = read_1dspec(str(file_name))
        usable_kernel = Gaussian1DKernel(int(sp))

        if kernel.lower() != "gaussian":
            if kernel.lower() == 'box':
                usable_kernel = Box1DKernel(int(sp))
            else:
                print ("Error: Kernel '{0}' Not Recognised".format(str(kernel)))
                sys.exit(1)

        smoothed_data = convolve(flux_data, usable_kernel)
        write_1dspec(ref_filename=str(file_name), flux_array=smoothed_data, prefix_str=str(prefix_str))

        if plot:
            plt.plot(wav_data, flux_data, 'g', label="Original Spectrum")
            plt.plot(wav_data, smoothed_data, 'r', label="Smooth Spectrum")
            plt.legend()
            plt.show()
            plt.close()


def scale_spectra(ctext, list_photfiles, prefix_str='f', plot=False):
    """
    Scales spectra acccording to the values in the array "scale_array". Basically, this step applies
    flux calibration on the spectra.
    Args:
        ctext           : Common text of 1-D Spectra files whose spectroscopic fluxes are to be scaled
        list_photfiles  : Text list of files containing different broadband photometric magnitudes
        prefix_str      : Prefix to distinguish the scaled 1-D spectra from the original
        plot            : Boolean describing whether the scaled spectra has to be plotted
    Returns:
        None
    """
    list_files = group_similar_files("", common_text=str(ctext))

    for file_name in list_files:
        dict_spec = read_specflux(str(file_name))
        dict_phot = read_photflux(list_photfiles=str(list_photfiles), julian_day=read_jd(str(file_name)))
        dict_modphot = dict((dict_centralwav[key], value) for (key, value) in dict_phot.items())
        dict_modspec = dict((dict_centralwav[key], value) for (key, value) in dict_spec.items())
        dict_modphot = get_zflux(dict_modphot, cntrl_wav=7000)

        common_keys = [key for key in dict_modphot.keys() if key in dict_modspec.keys()]
        dict_scale = dict((key, str(float(dict_modphot[key]) / float(dict_modspec[key]))) for key in common_keys)

        print (file_name, dict_modphot[7000], dict_modspec[7000])

        if len(dict_scale.keys()) > 2:
            order = len(dict_scale.keys()) - 1
            spline = UnivariateSpline(dict_scale.keys(), dict_scale.values(), k=int(order))

            with fits.open(str(file_name)) as hdulist:
                file_header = hdulist[0].header

            wav_data, flux_data = read_1dspec(str(file_name))
            scale_array = spline(wav_data)

            flux_moddata = np.multiply(np.asarray(flux_data), scale_array)
            output_filename = str(prefix_str) + str(file_name)
            fits.writeto(str(output_filename), data=flux_moddata, header=file_header, overwrite=True)

            if plot:
                x_new = np.arange(int(min(dict_scale.keys()) - 500), int(max(dict_scale.keys()) + 500), 5)
                plt.plot(dict_scale.keys(), dict_scale.values(), 'o', x_new, spline(x_new))
                plt.grid()
                plt.show()
                plt.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Accessing & Manipulating Text File Data
# ------------------------------------------------------------------------------------------------------------------- #

def get_file_dmnsn(file_name, title_rows=1):
    """
    Finds out the number of rows and columns in a text file.
    Args:
        file_name   : Text file whose dimensions are to be obtained
        title_rows  : No. of rows used as title description
    Returns:
        rows        : Number of rows in the text file
        columns     : Number of columns in the text file
    """
    with open(str(file_name), 'r') as f:
        columns = len(f.readline().rstrip().split())
    if columns == 0:
        print ("\nError : '{0}' Is An Empty File\n".format(str(file_name)))
        sys.exit(1)

    with open(str(file_name), 'r') as f:
        for _ in range(0, int(title_rows)):
            f.readline()
        length_data = len(f.read().split())
        rows = length_data / columns

    return rows, columns


def read_column(file_name, col_index, title_rows=1):
    """
    Extracts the specified column as a list from a text file.
    Args:
        file_name   : Text file from which the specified column has to be extracted
        col_index   : Index of the column to be extracted
        title_rows  : No. of rows used for title description
    Returns:
        list_col    : List of all the elements extracted from the column
    """
    rows, columns = get_file_dmnsn(str(file_name), title_rows=title_rows)

    with open(str(file_name), 'r') as f:
        for i in range(0, int(title_rows)):
            f.readline().rstrip()
        data_file = f.read().split()

    list_col = []
    for index in range(0, rows):
        list_col.append(data_file[col_index + index * columns])

    return list_col


def read_file(file_name, title_rows=1):
    """
    Extracts the file data as a list of columns from a text file.
    Args:
        file_name       : Text file from which file data has to be extracted
        title_rows      : No. of rows used for title description
    Returns:
        list_filedata   : List of all columns extracted from the text file
    """
    rows, columns = get_file_dmnsn(str(file_name), title_rows=title_rows)

    with open(str(file_name), 'r') as f:
        for i in range(0, int(title_rows)):
            f.readline().rstrip()
        data_file = f.read().split()

    list_filedata = []
    for col_index in range(0, columns):
        list_col = []
        for index in range(0, rows):
            list_col.append(data_file[col_index + index * columns])
        list_filedata.append(list_col)

    return list_filedata


def read_jd(file_name):
    """
    Reads JD of observation of the file "file_name".
    Args:
        file_name   : Name of the 1-D Spectra whose JD of observation is to be found out
    Returns:
        julian_day  : Julian day of the 1-D spectra
    """
    julian_day = fits.getval(str(file_name), str(JD_keyword))

    return julian_day


def read_specflux(file_name, file_specflux="output_sbands"):
    """
    Reads spectroscopic fluxes from a text file "file_specflux". The fluxes are determined
    for the spectra specified by the file "file_name"
    Args:
        file_name     : 1-D Spectra whose spectroscopic fluxes are to be extracted
        file_specflux : Text file from which spectroscopic fluxes are to be extracted
    Returns:
        dict_specflux : Dictionary of spectroscopic fluxes in different bands
    """
    dict_specflux = {}
    with open(str(file_specflux), 'r') as f:
        for line in f:
            line = line.rstrip()
            if re.search(str(file_name), line):
                dict_specflux[line.split()[1]] = line.split()[2]

    return dict_specflux


def read_photflux(list_photfiles, julian_day, flux=True):
    """
    Reads broadband photometric magnitudes from a text list containing names of files containing
    photometric magnitudes. The magnitudes are determined for the epoch specified by "julian_day".
    Args:
        list_photfiles  : Text list of files containing different broadband photometric magnitudes
        julian_day      : Julian day close to which the photometric magnitude has to be extracted
        flux            : Boolean describing whether flux has to be returned (or magnitude)
    Returns:
        dict_photmag    : Dictionary of photometric magnitudes in different bands
        dict_photflux   : Dictionary of photometric fluxes in different bands
    """
    list_photfiles = text_list_to_python_list(str(list_photfiles))

    dict_photmag = {}
    dict_photflux = {}
    for file_photmag in list_photfiles:
        file_data = read_file(str(file_photmag), title_rows=1)
        list_elmnts = file_data[0]
        list_mag = file_data[2]

        for index in range(0, len(list_elmnts)):
            if abs(float(list_elmnts[index]) - float(julian_day)) <= 0.25:
                dict_photmag[str.upper(file_photmag[-1:])] = list_mag[index]
                dict_photflux[str.upper(file_photmag[-1:])] = mag_to_flux(list_mag[index], file_band=file_photmag)
                break

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
    if file_band[-1].upper() == 'U':
        flux = "%7.5E" % 10 ** (-(float(mag) + zp_u + 21.100) / 2.5)
    elif file_band[-1].upper() == 'B':
        flux = "%7.5E" % 10 ** (-(float(mag) + zp_b + 21.100) / 2.5)
    elif file_band[-1].upper() == 'V':
        flux = "%7.5E" % 10 ** (-(float(mag) + zp_v + 21.100) / 2.5)
    elif file_band[-1].upper() == 'R':
        flux = "%7.5E" % 10 ** (-(float(mag) + zp_r + 21.100) / 2.5)
    elif file_band[-1].upper() == 'I':
        flux = "%7.5E" % 10 ** (-(float(mag) + zp_i + 21.100) / 2.5)
    else:
        print ("Error: Band Of Observation Not Recognised")
        sys.exit(1)

    return flux


def get_zflux(dict_phot, cntrl_wav=7000):
    """
    Obtains Z-band (narrow band, 7000 Angstroms) flux value.
    Args:
        dict_phot   : Dictionary containing broadband photometric flux values
        cntrl_wav   : Central wavelength of the Z-band
    Returns:
        dict_phot   : Modified dictionary with Z-band flux value included
    """
    if len(dict_phot.keys()) > 3:
        order = 3
    else:
        order = 2

    spline = UnivariateSpline(dict_phot.keys(), dict_phot.values(), k=int(order))
    dict_phot[int(cntrl_wav)] = "%7.5E" % spline(int(cntrl_wav))

    return dict_phot

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
copy_files(in_path=DIR_CURNT, out_path=DIR_SPECS, common_text="FILTER_*.asc")
copy_files(in_path=DIR_PHOT, out_path=DIR_SPECS, common_text="OUTPUT_InterpSNMag*")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files', choices=['Yes', 'No'])
ctext = eg.enterbox('Enter The Common Text Of Files To Be Flux Calibrated?', title='Flux Calibration',
                    default='cfwcbs_*.ms.fits')
bool_smooth = eg.boolbox('Perform Smoothening Of Spectra?', title='Smoothening 1-D Spectra', choices=['Yes', 'No'])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes Residual Files From Previous Run Of Flux Calibration
# ------------------------------------------------------------------------------------------------------------------- #
os.chdir(DIR_SPECS)
if rmv_files:
    for text in ['z_*.fits', 'fz_*.fits', 'tfz_*.fits', 'list_smspec']:
        remove_similar_files(common_text=str(text))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Smoothening, SBANDS Task & Finally Flux Calibration On 1-D Spectra
# ------------------------------------------------------------------------------------------------------------------- #
if bool_smooth:
    smooth_1dspec(ctext=str(ctext), sp=1.2, kernel="gaussian", prefix_str='z_')
    sbands(ctext="z_" + str(ctext))
else:
    sbands(ctext=str(ctext))

scale_spectra(ctext="OUTPUT_InterpSNMag*", list_photfiles="list_interpmag", plot=True)
# ------------------------------------------------------------------------------------------------------------------- #
