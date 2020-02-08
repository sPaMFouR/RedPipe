# !/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx---------------------------PLOT 1-D SPECTRA---------------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
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
from jdcal import gcal2jd
from astropy.io import fits
import matplotlib.pyplot as plt
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
JD_keyword = 'JD'
phase_plateauend = 110
phase_nebstart = 110
redshift = 0.008
date_explosion = 2457644.60
light_speed = 2.99792458e5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/Spec_Data/"
DIR_MODEL = "/home/avinash/Dropbox/IIP_Data/Model_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Converting Calendar Date To Julian Day
# ------------------------------------------------------------------------------------------------------------------- #

def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split('-')
    jd_tuple = gcal2jd(date_comp[0], date_comp[1], str(int(float(date_comp[2])) + 1))
    julian_day = jd_tuple[0] + jd_tuple[1]

    return julian_day

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

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Packages In IRAF
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.onedspec(_doprint=0)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Manipulating 1-D Spectra (Read, Write, Smoothen etc.)
# ------------------------------------------------------------------------------------------------------------------- #

def read_jd(file_name):
    """
    Reads JD of observation of the file 'file_name'.
    Args:
        file_name   : Name of the 1-D Spectra whose JD of observation is to be found out
    Returns:
        julian_day  : Julian day of the 1-D spectra
    """
    julian_day = float(fits.getval(filename=file_name, keyword=JD_keyword))

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
    list_spectra = group_similar_files('', common_text=common_text)
    usable_kernel = Gaussian1DKernel(int(sp))

    if kernel.lower() != 'gaussian':
        if kernel.lower() == 'box':
            usable_kernel = Box1DKernel(int(sp))
        else:
            print ("Error: Kernel '{0}' Not Recognised".format(kernel))
            sys.exit(1)

    for file_name in list_spectra:
        wav_data, flux_data = read_1dspec(file_name)
        smoothed_data = convolve(flux_data, usable_kernel)
        write_1dspec(ref_filename=file_name, flux_array=smoothed_data, prefix_str=prefix_str)

        if plot:
            plt.plot(wav_data, flux_data, 'g', label='Original Spectrum')
            plt.plot(wav_data, smoothed_data, 'r', label='Smooth Spectrum')
            plt.legend()
            plt.show()
            plt.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Conversion Of 1-D Spectra FITS File Into Text File
# ------------------------------------------------------------------------------------------------------------------- #

def wspectext(common_text, out_ext='.dat'):
    """
    Converts a list of 1-D spectra FITS files into a two-column text file.
    Args:
        common_text      : Common text of 1-D spectra FITS files
        out_ext          : Name of the extension to be used for the output file
    Returns:
        list_outspectra  : Python list of 1-D spectra DAT files
    """
    list_spectra = group_similar_files('', common_text=common_text)

    task = iraf.noao.onedspec.wspectext
    task.unlearn()

    list_outspectra = []
    for spectrum in list_spectra:
        output_spectrum = spectrum.split('.')[0] + out_ext
        task(input=spectrum, output=output_spectrum, header='no')
        list_outspectra.append(output_spectrum)

    return list_outspectra


def dopcor(common_text, prefix_str='r'):
    """
    Corrects a list of 1-D spectra FITS files for doppler redshift.
    Args:
        common_text      : Common text of 1-D spectra FITS files
        prefix_str       : Prefix to distinguish the doppler corrected 1-D spectra from the original
    Returns:
        list_outspectra  : Python list of 1-D spectra FITS files corrected for doppler redshift
    """
    list_spectra = group_similar_files('', common_text=common_text)

    task = iraf.noao.onedspec.dopcor
    task.unlearn()
    
    task.isvelocity = 'no'
    task.dispersion = 'yes'
    task.flux = 'no'

    list_outspectra = []
    for spectrum in list_spectra:
        output_spectrum = prefix_str + spectrum
        task(input=spectrum, output=output_spectrum, redshift=redshift)
        list_outspectra.append(output_spectrum)

    return list_outspectra


def convert_to_text(common_text, out_ext='.dat'):
    """
    Converts a list of 1-D spectra FITS files into a two-column text file.
    Args:
        common_text      : Common text of 1-D spectra FITS files
        out_ext          : Name of the extension to be used for the output file
    Returns:
        list_outspectra  : Python list of 1-D spectra DAT files
    """
    list_spectra = group_similar_files('', common_text=common_text)
    
    list_outspectra = []
    for spectrum in list_spectra:
        wave_array, flux_array = read_1dspec(spectrum)
        spec_df = pd.DataFrame(flux_array, index=wave_array)
        
        output_spectrum = spectrum.split('.')[0] + out_ext
        spec_df.to_csv(output_spectrum, sep=' ', index=True, header=None)
        list_outspectra.append(output_spectrum)

    return list_outspectra

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Subplots
# ------------------------------------------------------------------------------------------------------------------- #

def plot_singlespec(ax_obj, file_name):
    """
    Plots the spectrum with line identified and labelled accordingly.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
    Returns:
        None
    """
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]

    ax_obj.plot(data_df['Wavelength'], data_df['Flux'] * 1e16)
    ax_obj.grid()
#     ax_obj.tick_params(axis='y', which='both', labelleft='off')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.set_xlim(3400, 9300)
    ax_obj.xaxis.set_major_locator(MultipleLocator(500))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(100))
    ax_obj.tick_params(axis='x', which='both', direction='in', width=0.5, labelsize=12)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files', choices=['Yes', 'No'])
# ctext = eg.enterbox('Common Text Of Files To Be Flux Calibrated?', title='Flux Calibration', default='fz_*.fits')
# bool_smooth = eg.boolbox('Perform Smoothening Of Spectra?', title='Smoothening 1-D Spectra', choices=['Yes', 'No'])
# clip_str = eg.enterbox('Specify Clipping Section: ', title='Clipping Of 1-D Spectra', default='4000:8500')
rmv_files = True
ctext = 'fz_*.fits'
bool_smooth = False
clip_str = '3500:9100'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes Residual Files From Previous Run Of Plotting Spectra
# ------------------------------------------------------------------------------------------------------------------- #
if rmv_files:
    for text in ['rfz_*.fits', '*.dat']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Smoothening & Writes The 1-D Spectra Into A Text File (.dat)
# ------------------------------------------------------------------------------------------------------------------- #
list_spec = dopcor(common_text=ctext)

if bool_smooth:
    smooth_1dspec(common_text='r' + ctext, sp=1.5, kernel='gaussian', prefix_str='z_', plot=False)
    list_textspec = convert_to_text('z_r' + ctext)
else:
    list_textspec = convert_to_text('r' + ctext)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Creates A Pandas DataFrame With Information Of All 1-Dimensional Spectra
# ------------------------------------------------------------------------------------------------------------------- #
spec_df = pd.DataFrame(list_spec, columns=['FileName'])
spec_df['TextFile'] = list_textspec
spec_df['Phase'] = spec_df['FileName'].apply(lambda x: read_jd(x) - date_explosion)
spec_df['Label'] = spec_df['Phase'].apply(lambda x: '+{0:>.1f} d'.format(x) if x > 0 else '{0:>.1f} d'.format(x))
spec_df = spec_df.sort_values(by='Phase').reset_index(drop=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
early_df = spec_df[spec_df['Phase'] < 20].copy()
plateau_df = spec_df[spec_df['Phase'] < phase_plateauend].copy()
nebular_df = spec_df[spec_df['Phase'] > phase_nebstart].copy()
halpha_df = spec_df[spec_df['Phase'] < phase_nebstart].copy()
na1d_df = spec_df[1:].copy()

lower_lim, upper_lim = clip_str.split(':')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)

for index, file_name in plateau_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= 1.0e-15 * index
        
    ax.plot(data_df['Wavelength'], data_df['Flux'], label='_nolegend_')
    ax.text(x=data_df['Wavelength'].values[-1] + 100, y=data_df['Flux'].values[-1],
            s=plateau_df.loc[index, 'Label'], fontsize=10)
    
    if file_name == 'rfz_jul21_cfwcbs_ASASSN14dq-gr7.dat':
        jul21_df = data_df.copy()
    elif file_name == 'rfz_jul24_cfwcbs_ASASSN14dq-gr7.dat':
        jul24_df = data_df.copy()
    elif file_name == 'rfz_jul30_cfwcbs_ASASSN14dq.dat':
        jul30_df = data_df.copy()
    elif file_name == 'rfz_aug02_cfwcbs_ASASSN14dq.dat':
        aug02_df = data_df.copy()

dict_labels = {r'$\rm Ca\,II$ (H & K)': [3840, 60], r'$\rm H\gamma$ 4340': [4210, 70], 
               r'$\rm H\beta$ 4861': [4730, 70], r'$\rm Fe\,II$ Lines': [5000, 140], r'$\rm Na\,I$ D': [5810, 50],
               r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5460, 50], r'$\rm Sc\,II$ Multiplet': [5610, 50], 
               r'$\rm H\alpha$ 6563': [6380, 100], r'[$\rm O_2$ 6867 B band]': [6820, 50], 
               r'[$\rm O_2$ 7620 A band]': [7560, 50], r'$\rm Ca\,II$ Triplet': [8480, 180]}

ax.set_ylim(-20e-15, 6.3e-15)
ax.set_xlim(int(lower_lim) - 400, int(upper_lim) + 800)
ax.axvline(6563, color='k', linestyle='--', linewidth=0.7, alpha=0.2, label='_nolegend_')
ax.axvline(6712, color='k', linestyle='--', linewidth=0.7, alpha=0.2, label='Galactic Lines')

for (line, [wavelength, width]) in dict_labels.items():
    ax.text(wavelength - 30, 5.8e-15, line, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.1)

ax.legend(fontsize=10)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotSpecPlateau.eps', format='eps', dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Early Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(11, 7))
ax = fig.add_subplot(111)

for index, file_name in early_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= 1.2e-15 * index
        
    ax.plot(data_df['Wavelength'], data_df['Flux'], label='_nolegend_')
    ax.text(x=data_df['Wavelength'].values[-1] + 100, y=data_df['Flux'].values[-1],
            s=plateau_df.loc[index, 'Label'], fontsize=12)

dict_labels1 = {r'$\rm Ca\,II$ (H & K)': [3800, 60], r'$\rm H\delta$ 4102': [3960, 50], 
                r'$\rm H\gamma$ 4340': [4190, 50], r'$\rm H\beta$ 4861': [4700, 70], r'$\rm He\,I$ 5876': [5670, 80], 
                r'$\rm H\alpha$ 6563': [6340, 110], r'[$\rm O_2$ 7620 band]': [7560, 60], r'$\rm N\,I$': [6583, 50]}

ax.set_ylim(-4.8e-15, 4.5e-15)
ax.set_xlim(int(lower_lim) - 400, int(upper_lim) + 800)
ax.axvline(6563, color='k', linestyle='--', linewidth=0.7, alpha=0.1, label='_nolegend_')
ax.axvline(6712, color='k', linestyle='--', linewidth=0.7, alpha=0.1, label='Galactic Lines')

for (line, [wavelength, width]) in dict_labels1.items():
    ax.text(wavelength - 30, 4.2e-15, line, fontsize=9, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.1)

ax.legend(fontsize=11)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotSpecEarly.eps', format='eps', dpi=800, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Nebular Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)

for index, file_name in nebular_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= 2.3e-16 * index

    ax.plot(data_df['Wavelength'], data_df['Flux'], label='_nolegend_')
    ax.text(x=data_df['Wavelength'].values[-1] + 100, y=data_df['Flux'].values[-1],
            s=nebular_df.loc[index, 'Label'], fontsize=11)

dict_labels2 = {r'$\rm Na\,I$ D': [5820, 60], r'$\rm Fe\,II$ Lines': [4990, 160], r'$\rm H\alpha$ 6563': [6440, 50],
                r'$\rm [O\,I]$ 6300, 6364': [6290, 50], r'$\rm [Ca\,II]$ 7291, 7324': [7290, 50], 
                r'[$\rm O_2$ 7620 A band]': [7560, 50], r'$\rm Ca\,II$ Triplet': [8560, 130]}

ax.set_ylim(-7.5e-15, -3.2e-15)
ax.set_xlim(int(lower_lim) - 400, int(upper_lim) + 800)
ax.axvline(6563, color='k', linestyle='--', linewidth=0.7, alpha=0.2, label='_nolegend_')
ax.axvline(6722, color='k', linestyle='--', linewidth=0.7, alpha=0.2, label='Galactic Lines')

for (line, [wavelength, width]) in dict_labels2.items():
    ax.text(wavelength - 30, -3.3e-15, line, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.1)
    
ax.legend(fontsize=10)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotSpecNebular.eps', format='eps', dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The H-Alpha Feature In Spectra Of Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 14))
ax = fig.add_subplot(111)

lower_clip = 6200
upper_clip = 6800

for index, file_name in halpha_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_clip)) & (data_df['Wavelength'] <= int(upper_clip))]
    data_df['Velocity'] = data_df['Wavelength'].apply(lambda x: ((x - 6563) / 6563) * 3e2)
    data_df['Flux'] -= 5e-16 * index

    ax.plot(data_df['Wavelength'], data_df['Flux'])
    ax.text(x=data_df['Wavelength'].values[0] - 80, y=data_df['Flux'].values[0], 
            s=halpha_df.loc[index, 'Label'], fontsize=10)

ax.set_ylim(-9.5e-15, 1.3e-15)
ax.set_xlim(lower_clip - 110, upper_clip + 40)
ax.axvline(6563, linestyle='--', color='k', linewidth=0.7)
ax.axvline(6582, linestyle='--', color='k', linewidth=0.7)
ax.axvline(6717, linestyle='--', color='k', linewidth=0.7)
ax.axvline(6730, linestyle='--', color='k', linewidth=0.7)

ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(100))
ax.xaxis.set_minor_locator(MultipleLocator(25))
ax.tick_params(axis='y', which='both', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=0.5, labelsize=12)
# ax.set_xlabel(r'Velocity [$\rm 10^3 km\ s^{-1}$]', fontsize=12)
ax.set_xlabel(r'Wavelength ($\AA$)', fontsize=12)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=12)

fig.savefig('OUTPUT_PlotHalpha.eps', format='eps', dpi=1200, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The High Velocity Feature
# ------------------------------------------------------------------------------------------------------------------- #
jul21_df['Flux'] = jul21_df['Flux'].apply(lambda x: x / 1e-16)
jul24_df['Flux'] = jul24_df['Flux'].apply(lambda x: x / 1e-16 + 1)
jul30_df['Flux'] = jul30_df['Flux'].apply(lambda x: x / 1e-16 + 3)
aug02_df['Flux'] = aug02_df['Flux'].apply(lambda x: x / 1e-16 + 5)

jul24_df['HBetaVel'] = jul24_df['Wavelength'].apply(lambda x: (x - 4861) * 299.792458 / 4861)
jul21_df['HAlphaVel'] = jul21_df['Wavelength'].apply(lambda x: (x - 6563) * 299.792458 / 6563)
jul24_df['HAlphaVel'] = jul24_df['Wavelength'].apply(lambda x: (x - 6563) * 299.792458 / 6563)
jul30_df['HAlphaVel'] = jul30_df['Wavelength'].apply(lambda x: (x - 6563) * 299.792458 / 6563)
aug02_df['HAlphaVel'] = aug02_df['Wavelength'].apply(lambda x: (x - 6563) * 299.792458 / 6563)

fig_hv, (ax_hbeta, ax_halpha) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [2, 5]}, figsize=(6, 9))

ax_hbeta.plot(jul24_df['HBetaVel'], jul24_df['Flux'], color='k')
ax_hbeta.set_xlim(-24, 12)
ax_hbeta.set_ylim(2, 16)
ax_hbeta.set_yticklabels([])
ax_hbeta.set_xticklabels([])
ax_hbeta.yaxis.set_ticks_position('both')
ax_hbeta.xaxis.set_ticks_position('both')
ax_hbeta.xaxis.set_major_locator(MultipleLocator(5))
ax_hbeta.xaxis.set_minor_locator(MultipleLocator(2.5))
ax_hbeta.tick_params(which='both', direction='in', width=1, labelsize=14)

ax_halpha.plot(jul21_df['HAlphaVel'], jul21_df['Flux'], linewidth=1, color='k')
ax_halpha.plot(jul24_df['HAlphaVel'], jul24_df['Flux'], linewidth=1, color='k')
ax_halpha.plot(jul30_df['HAlphaVel'], jul30_df['Flux'], linewidth=1, color='k')
ax_halpha.plot(aug02_df['HAlphaVel'], aug02_df['Flux'], linewidth=1, color='k')

ax_halpha.set_xlim(-24, 12)
ax_halpha.set_ylim(-16, 12)
ax_halpha.set_yticklabels([])
ax_halpha.yaxis.set_ticks_position('both')
ax_halpha.xaxis.set_ticks_position('both')
ax_halpha.xaxis.set_major_locator(MultipleLocator(5))
ax_halpha.xaxis.set_minor_locator(MultipleLocator(2.5))
ax_halpha.tick_params(which='both', direction='in', width=0.8, labelsize=14)
ax_halpha.set_xlabel(r'Velocity ($\rm 10^3\ km\ s^{-1}$)', fontsize=14)

ax_hbeta.text(5, 10, '+21.9 d', fontsize=14)
ax_halpha.text(5, 5, '+18.9 d', fontsize=14)
ax_halpha.text(6, -1, '+21.9 d', fontsize=14)
ax_halpha.text(6, -6, '+27.9 d', fontsize=14)
ax_halpha.text(6, -14, '+30.9 d', fontsize=14)

ax_hbeta.text(7, 13, s=r'H$\beta$', fontsize=16)
ax_halpha.text(7, 10, s=r'H$\alpha$', fontsize=16)
ax_hbeta.text(-16.5, 13, s=r'B', fontsize=16)
ax_halpha.text(-16.5, 5, s=r'A', fontsize=16)
ax_halpha.text(-9, -15, s=r'NC', fontsize=16)
ax_halpha.text(-16.5, -15, s=r'HV', fontsize=16)

ax_hbeta.axvline(x=0, linestyle='-', color='k', linewidth=1)
ax_halpha.axvline(x=0, linestyle='-', color='k', linewidth=1)
ax_hbeta.axvline(x=-17, linestyle='--', color='k', linewidth=1)
ax_halpha.axvline(x=-17, linestyle='--', color='k', linewidth=1)
ax_hbeta.axvline(x=-9.3, linestyle='--', color='k', linewidth=1)
ax_halpha.axvline(x=-9.3, linestyle='--', color='k', linewidth=1)

fig_hv.text(0.04, 0.5, r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', rotation='vertical', fontsize=16)
fig_hv.subplots_adjust(hspace=0.000, top=0.9, right=0.95)
fig_hv.savefig('OUTPUT_PlotHVFeature.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig_hv)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Spectra In Comparison With Models From Dessart
# ------------------------------------------------------------------------------------------------------------------- #

def plot_modspec(ax_obj, file_name, phase='0d', offset=0, z=r'1 $\rm Z_{\odot}$'):
    data_df = pd.read_csv(file_name, sep='\s+', engine='python', header=None, comment='#')
    data_df = data_df[(data_df[0] > 3800)]
    data_df[1] = data_df[1] / data_df[1].mean()

    if file_name[0:3] != 'rfz':
        label = file_name.split('/')[-2]
        color = 'k'
        ax_obj.text(5250, offset + 0.4, s=label + ' +' + phase, color=color, fontsize=10)
        ax_obj.text(7000, offset + 0.3, s=z, color=color, fontsize=10)

    else:
        label = name_SN
        color = 'r'
        ax_obj.text(7250, 0.3, s=label + ' +' + phase, color=color, fontsize=10)
        data_df[1] = convolve(data_df[1].tolist(), Gaussian1DKernel(3))

        for index, row in data_df.iterrows():
            if 6500 < row[0] < 7500:
                row[1] = row[1] * 1.15
            if 8200 > row[0] >= 7500:
                row[1] = row[1] * 0.8
            elif row[0] >= 8200:
                row[1] = row[1] * 0.5

    ax_obj.plot(data_df[0], np.log10(data_df[1]) + offset, linewidth=1, c=color, alpha=0.8, label=label + ' ' + phase)


fig = plt.figure(figsize=(18, 9))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122, sharey=ax)

plot_modspec(ax, 'rfz_aug08_cfwcbs_ASASSN14dq.dat', '36.7 d', offset=0.05)
plot_modspec(ax, DIR_MODEL + 'm15z2m3/m15_du_sch_z2m3_FeC_mix0p4_15.fl', '38.1 d', -0.9, z=r'0.1 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z8m3/m15_du_sch_z8m3_FeC_mix0p4_15.fl', '37.4 d', -0.5, z=r'0.4 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z2m2/m15_du_sch_FeC_mix0p4_15.fl', '37.1 d', -1.3, z=r'1 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z4m2/m15_du_sch_z4m2_FeC_mix0p4_15.fl', '37.6 d', -1.7, z=r'2 $\rm Z_{\odot}$')

plot_modspec(ax2, 'rfz_aug31_cfwcbs_ASASSN14dq.dat', '59.7 d', offset=0.05)
plot_modspec(ax2, DIR_MODEL + 'm15z2m3/m15_du_sch_z2m3_FeC_mix0p4_20.fl', '61.3 d', -0.9, z=r'0.1 $\rm Z_{\odot}$')
plot_modspec(ax2, DIR_MODEL + 'm15z8m3/m15_du_sch_z8m3_FeC_mix0p4_20.fl', '60.2 d', -0.5, z=r'0.4 $\rm Z_{\odot}$')
plot_modspec(ax2, DIR_MODEL + 'm15z2m2/m15_du_sch_FeC_mix0p4_20.fl', '59.7 d', -1.3, z=r'1 $\rm Z_{\odot}$')
plot_modspec(ax2, DIR_MODEL + 'm15z4m2/m15_du_sch_z4m2_FeC_mix0p4_20.fl', '60.6 d', -1.7, z=r'2 $\rm Z_{\odot}$')

ax.set_xlim(3700, 9500)
ax.set_ylim(-1.9, 0.5)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=14)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
ax.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=15)

ax2.set_xlim(3700, 9500)
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_major_locator(MultipleLocator(1000))
ax2.xaxis.set_minor_locator(MultipleLocator(250))
ax2.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax2.tick_params(axis='x', which='both', direction='in', width=1, labelsize=14)
ax2.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)

fig.subplots_adjust(wspace=0.01)
fig.savefig('OUTPUT_PlotCompSpecMetal.eps', format='eps', dpi=1200, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
