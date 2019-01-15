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
name_SN = 'ASASSN-14dq'
JD_keyword = 'JD'
phase_plateauend = 90
phase_nebstart = 105
date_explosion = 2456841.50
light_speed = 2.99792458e5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SPEC = "/home/avinash/Dropbox/IIP_Data/Spec_Data/"
DIR_MODEL = "/home/avinash/Dropbox/IIP_Data/Model_Data/"
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
# Function For Conversion Of 1-D Spectra FITS File Into Text File
# ------------------------------------------------------------------------------------------------------------------- #

def wspectext(text_list_spectra, out_ext='.dat'):
    """
    Converts a list of 1-D spectra FITS files into a two-column text file.
    Args:
        text_list_spectra   : Text list of 1-D spectra FITS files
        out_ext             : Name of the extension to be used for the output file
    Returns:
        list_outspectra     : Python list of 1-D spectra DAT files
    """
    list_spectra = text_list_to_python_list(text_list_spectra)

    task = iraf.noao.onedspec.wspectext
    task.unlearn()

    list_outspectra = []
    for spectrum in list_spectra:
        output_spectrum = spectrum.split('.')[0] + out_ext
        task(input=spectrum, output=output_spectrum, header='no')
        list_outspectra.append(output_spectrum)

    return list_outspectra


def convert_to_text(text_list_spectra, out_ext='.dat'):
    """
    Converts a list of 1-D spectra FITS files into a two-column text file.
    Args:
        text_list_spectra   : Text list of 1-D spectra FITS files
        out_ext             : Name of the extension to be used for the output file
    Returns:
        None
    """
    
    list_spectra = text_list_to_python_list(text_list_spectra)
    
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


def smooth_1dspec(text_list_spectra, sp=1.2, kernel='gaussian', prefix_str='z_', plot=False):
    """
    Smoothens a 1-D spectra based on the smoothening parameter. Smoothening parameter
    is 'std.dev.' in case of isotropic Gaussian filter and is 'width' in the case of the
    non-isotropic box filter.
    Args:
        text_list_spectra   : Common text of 1-D spectra files which have to be smoothened
        sp                  : Smoothening parameter
        kernel              : Convolution Kernel used for smoothening (Gaussian or Box)
        prefix_str          : Prefix to distinguish the smoothened 1-D spectra from the original
        plot                : Boolean describing whether the smoothened spectra has to be plotted
    Returns:
        None
    """
    list_spectra = text_list_to_python_list(text_list_spectra)
    usable_kernel = Gaussian1DKernel(int(sp))

    if kernel.lower() != 'gaussian':
        if kernel.lower() == 'box':
            usable_kernel = Box1DKernel(int(sp))
        else:
            print ("Error: Kernel '{0}' Not Recognised".format(str(kernel)))
            sys.exit(1)

    for file_name in list_spectra:
        wav_data, flux_data = read_1dspec(file_name)
        smoothed_data = convolve(flux_data, usable_kernel)
        write_1dspec(ref_filename=file_name, flux_array=smoothed_data, prefix_str=prefix_str)

        if plot:
            plt.plot(wav_data, flux_data, 'g', label="Original Spectrum")
            plt.plot(wav_data, smoothed_data, 'r', label="Smooth Spectrum")
            plt.legend()
            plt.show()
            plt.close()

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
ctext = '*.fits'
bool_smooth = False
clip_str = '3700:9100'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes Residual Files From Previous Run Of Plotting Spectra
# ------------------------------------------------------------------------------------------------------------------- #
if rmv_files:
    for text in ['*.dat']:
        remove_similar_files(common_text=text)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Smoothening & Writes The 1-D Spectra Into A Text File (.dat)
# ------------------------------------------------------------------------------------------------------------------- #
text_list_spectra = 'list_spectra'
lower_lim, upper_lim = clip_str.split(':')
list_fits = group_similar_files(text_list_spectra, str(ctext))

if bool_smooth:
    smooth_1dspec(text_list_spectra, sp=1.5, kernel='gaussian', prefix_str='z_', plot=False)
    group_similar_files(text_list_spectra, common_text='z_' + str(ctext))

list_spectra = wspectext(text_list_spectra)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Creates A Pandas DataFrame With Information Of All 1-Dimensional Spectra
# ------------------------------------------------------------------------------------------------------------------- #
spec_df = pd.DataFrame(list_fits, columns=['FileName'])
spec_df['TextFile'] = list_spectra
spec_df['Phase'] = spec_df['FileName'].apply(lambda x: read_jd(x) - date_explosion)
spec_df['Label'] = spec_df['Phase'].apply(lambda x: '+{0:>.1f}d'.format(x) if x > 0 else '{0:>.1f}d'.format(x))
spec_df = spec_df.sort_values(by='Phase')
spec_df = spec_df.reset_index(drop=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
plateau_df = spec_df[spec_df['Phase'] < phase_plateauend].copy()
nebular_df = spec_df[spec_df['Phase'] > phase_nebstart].copy()
halpha_df = spec_df[spec_df['Phase'] > phase_nebstart][1:].copy()
na1d_df = spec_df[1:].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

for index, file_name in plateau_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= 0.6e-15 * index
        
    ax.plot(data_df['Wavelength'], data_df['Flux'], c='k', linewidth=1.0, alpha=0.8)
    ax.text(x=data_df['Wavelength'].values[-1] + 100, y=data_df['Flux'].values[-1],
            s=plateau_df.loc[index, 'Label'], fontsize=11)
    
    if file_name == 'rfz_jul21_cfwcbs_ASASSN14dq-gr7.dat':
        jul21_df = data_df.copy()
    elif file_name == 'rfz_jul24_cfwcbs_ASASSN14dq-gr7.dat':
        jul24_df = data_df.copy()
    elif file_name == 'rfz_jul30_cfwcbs_ASASSN14dq.dat':
        jul30_df = data_df.copy()
    elif file_name == 'rfz_aug02_cfwcbs_ASASSN14dq.dat':
        aug02_df = data_df.copy()

dict_labels = {r'H$\delta$ 4102': [3970, 1.40], r'H$\gamma$ 4340': [4210, 1.20], r'H$\beta$ 4861': [4715, 1.6],
               r'$\rm Fe\,II 5169$': [5060, 1.4], r'$\rm He\,I 5876$': [5720, 0.60], r'H$\alpha$ 6563': [6380, 1.6]}

ax.annotate(r'A', xy=(6220, 3.5e-16), xytext=(6170, 6e-16), arrowprops=dict(arrowstyle='-'), rotation='vertical')
ax.annotate(r'A', xy=(6220, -3e-16), xytext=(6170, 0), arrowprops=dict(arrowstyle='-'), rotation='vertical')
ax.annotate(r'A', xy=(6220, -9e-16), xytext=(6170, -6e-16), arrowprops=dict(arrowstyle='-'), rotation='vertical')
ax.annotate(r'A', xy=(6220, -15e-16), xytext=(6170, -12e-16), arrowprops=dict(arrowstyle='-'), rotation='vertical')
ax.annotate(r'B', xy=(4590, 1.1e-15), xytext=(4540, 1.4e-15), arrowprops=dict(arrowstyle='-'), rotation='vertical')

for (line, [wavelength, flux]) in dict_labels.items():
    ax.annotate(line, xy=(wavelength, flux * 1e-15), xytext=(wavelength - 50, 3.25e-15),
                arrowprops=dict(arrowstyle='-'), rotation='vertical')

ax.set_ylim(-7.5e-15, 3.5e-15)
ax.set_xlim(int(lower_lim) - 200, int(upper_lim) + 600)
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotSpecPlateau.eps', format='eps', dpi=1000, bbox_inches='tight')
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

fig_hv.text(0.04, 0.5, r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', va='center',
            rotation='vertical', fontsize=16)
fig_hv.subplots_adjust(hspace=0.000, top=0.9, right=0.95)
fig_hv.savefig('OUTPUT_PlotHVFeature.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig_hv)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Nebular Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111)

for index, file_name in nebular_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= 2e-16 * index

    ax.plot(data_df['Wavelength'], data_df['Flux'], c='k', linewidth=0.9, alpha=0.8)
    ax.text(x=data_df['Wavelength'].values[-1] + 100, y=data_df['Flux'].values[-1],
            s=nebular_df.loc[index, 'Label'], fontsize=11)

dict_labels2 = {r'$\rm He\,I 5876$': [5820, 2.80], r'H$\alpha$ 6563': [6380, 2.8], 
                r'$\rm Ca\,II 8498 8542$': [8430, 3.0], r'$\rm Ca\,II 8662$': [8600, 3.00]}

ax.text(4700, -2.7e-15, 'FeII lines')
for (line, [wavelength, flux]) in dict_labels2.items():
    ax.annotate(line, xy=(wavelength, flux * -1e-15), xytext=(wavelength - 50, -2.2e-15),
                arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')

ax.set_ylim(-4.5e-15, -2.1e-15)
ax.set_xlim(int(lower_lim) - 200, int(upper_lim) + 600)
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=0.5, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotSpecNebular.eps', format='eps', dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The H-Alpha Feature In Spectra Of Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(7, 9))
ax = fig.add_subplot(111)
halpha_df = spec_df[spec_df['Phase'] > 100][1:].copy()
lower_clip = 6350
upper_clip = 6800

for index, file_name in halpha_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_clip)) & (data_df['Wavelength'] <= int(upper_clip))]
    data_df['Velocity'] = data_df['Wavelength'].apply(lambda x: ((x - 6563) / 6563) * 3e2)
    data_df['Flux'] -= 1e-16 * index

    ax.plot(data_df['Velocity'], data_df['Flux'], c='k', linewidth=1.2, alpha=0.8)
    ax.text(x=data_df['Velocity'].values[0] - 1000, y=data_df['Flux'].values[0] * 0.995, 
            s=halpha_df.loc[index, 'Label'], color='k', fontsize=14)

ax.set_xlim(-11.5, 11)
ax.axvline(0, linestyle='--', color='k', linewidth=1)
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(axis='y', which='both', direction='in', width=1, labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=14)
ax.set_xlabel(r'Velocity [$\rm 10^3 km\ s^{-1}$]', fontsize=14)
ax.set_ylabel(r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('OUTPUT_PlotHalpha.eps', format='eps', dpi=1200, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The NaI D Feature In Spectra Of Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 12))
ax = fig.add_subplot(111)

lower_clip = 5850
upper_clip = 5980

for index, file_name in na1d_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[(data_df['Wavelength'] >= int(lower_clip)) & (data_df['Wavelength'] <= int(upper_clip))]
    data_df['Flux'] = data_df['Flux'] / data_df['Flux'].mean()
    data_df['Flux'] -= 0.3 * index

    ax.plot(data_df['Wavelength'], data_df['Flux'])
    ax.text(x=data_df['Wavelength'].values[-1] + 5, y=data_df['Flux'].values[-1], s=na1d_df.loc[index, 'Label'])

ax.set_xlim(lower_clip - 10, upper_clip + 30)
ax.tick_params(axis='y', which='both', labelleft='off')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.tick_params(axis='x', which='both', direction='in', width=0.2, labelsize=10)
ax.set_xlabel('Rest Wavelength ($\AA$)', fontsize=12)
ax.set_ylabel('Scaled Flux ($erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$) + Const.', fontsize=12)

fig.savefig('OUTPUT_PlotNaID.eps', format='eps', dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra With The Highest SNR
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(15, 6))
ax = fig.add_subplot(111)

group_similar_files('single_spec', common_text='z_jul30_cfwcbs_ASASSN14dq.ms.fits')
file_name = wspectext('single_spec')[0]
data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]

ax.plot(data_df['Wavelength'], data_df['Flux'], label=r'+$28$ days')
# ax.text(9200, 2e-15, r'+28 d')

dict_arrow = dict(facecolor='blue', arrowstyle='-')
dict_lines = {r'$\rm Ca\,II\ 3934,3968$': [3890, 0.70], r'H$\delta$ 4102': [4020, 1.30],
              r'H$\gamma$ 4340': [4260, 0.80], r'$\rm Ba\,II\ 4554 + Fe\,II$': [4500, 1.10],
              r'H$\beta$ 4861': [4775, 1.00],r'$\rm Fe\,II\ 4924$': [4860, 1.48],
              r'$\rm Fe\,II\ 5018 + Sc\,II$': [4970, 1.30], r'$\rm Fe\,II\ 5169$': [5110, 1.05],
              r'$\rm Fe\,II\ 5535 + Sc\,II\ 5527$': [5470, 1.10], r'$\rm Sc\,II\ multiplet$': [5620, 1.05],
              r'$\rm Na\,ID\ 5890,5896$': [5890, 1.00], '(HV)': [6300, 0.80],r'H$\alpha$ 6563': [6430, 0.75],
              r'[O$_2$] 6867 (Atm.)': [6870, 0.75],r'[H$_2$O] 7165 (Atm.)': [7175, 0.70],
              r'[O$_2$] 7620 (Atm.)': [7620, 0.60], 'NI': [8150, 0.55], r'$\rm Ca\,II\ 8498,8542$': [8410, 0.45],
              r'$\rm Ca\,II\ 8662$': [8550, 0.50]}

for (line, [wavelength, flux]) in dict_lines.items():
    ax.annotate(line, xy=(wavelength, flux * 1e-15), xytext=(wavelength - 30, 2.0e-15), arrowprops=dict_arrow,
                rotation='vertical')

ax.legend(fontsize=14)
ax.set_ylim(0.2e-15, 2.2e-15)
ax.set_xlim(int(lower_lim) - 400, int(upper_lim) + 500)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.20e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.05e-15))
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(axis='y', which='both', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=0.5, labelsize=16)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=16)

fig.savefig('OUTPUT_PlotSingleSpec.eps', format='eps', dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra Across Different Epochs For Line Identification
# ------------------------------------------------------------------------------------------------------------------- #
dict_arrow = dict(facecolor='blue', arrowstyle='-')

dict_lines2 = {'CaII 3934,3968': [3890, 0.70], r'H$\gamma$ 4340': [4260, 0.80], 'BaII 4554 + FeII': [4500, 1.10],
               r'H$\beta$ 4861': [4775, 1.00], 'FeII 4924': [4920, 1.45], 'FeII 5018 + ScII': [5020, 1.30],
               'FeII 5169': [5160, 1.05], 'HeI 5876': [5835, 1.00],
               'NaID 5890,5896': [5896, 1.00], r'H$\alpha$ 6563': [6485, 0.75], r'[O$_2$ 7620 band]': [7620, 0.60],
               'CaII 8498 8542': [8500, 0.80], 'CaII 8662': [8660, 1.2]}

dict_lines3 = {'CaII 3934,3968': [3890, 0.70], r'H$\gamma$ 4340': [4260, 0.80], 'BaII 4554 + FeII': [4500, 1.10],
               r'H$\beta$ 4861': [4775, 1.00], 'FeII 4924': [4860, 1.45], 'FeII 5018 + ScII': [4970, 1.30],
               'FeII 5169': [5110, 1.05], 'FeII 5535 + ScII 5527': [5450, 1.10], 'HeI 5876': [5820, 1.00],
               'NaID 5890,5896': [5890, 1.00], r'H$\alpha$ 6563': [6430, 0.75], r'[O$_2$ 6867 band]': [6870, 0.75],
               '[CaII 7291, 7324]': [7300, 0.70], r'[O$_2$ 7620 band]': [7620, 0.60], 'NI': [8150, 0.48],
               'CaII 8498 8542': [8410, 0.45], 'CaII 8662': [8550, 0.50]}

list_files = ['z_jul30_cfwcbs_ASASSN14dq.ms.fits', 'nov01_cfwcbs_ASASSN14dq.ms.fits.2.fits',
              'dec08_cfwcbs_ASASSN14dq.ms.fits.1']

python_list_to_text_list(list_files, 'list_singlespec')
list_specfiles = wspectext('list_singlespec')

fig = plt.figure(figsize=(14, 18))
fig.subplots_adjust(hspace=0.001, top=0.9, right=0.95)

ax1 = fig.add_subplot(311)
plot_singlespec(ax1, list_specfiles[0])
ax1.set_ylim(2, 22)
ax1.yaxis.set_major_locator(MultipleLocator(2))
ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
for (line, [wave, flux]) in dict_lines.items():
    ax1.annotate(line, xy=(wave, flux * 10), xytext=(wave - 30, 20), arrowprops=dict_arrow, rotation='vertical')

ax2 = fig.add_subplot(312, sharex=ax1)
plot_singlespec(ax2, list_specfiles[1])
ax2.set_ylim(-0.8, 5.4)
ax2.yaxis.set_major_locator(MultipleLocator(0.6))
ax2.yaxis.set_minor_locator(MultipleLocator(0.15))
ax2.set_ylabel('Flux ($erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$) + Const.', fontsize=16)
for (line, [wave, flux]) in dict_lines2.items():
    ax2.annotate(line, xy=(wave, flux * 1), xytext=(wave - 30, 5.0), arrowprops=dict_arrow, rotation='vertical')

ax3 = fig.add_subplot(313, sharex=ax1)
plot_singlespec(ax3, list_specfiles[2])
ax3.set_ylim(-0.4, 3.2)
ax3.yaxis.set_major_locator(MultipleLocator(0.40))
ax3.yaxis.set_minor_locator(MultipleLocator(0.10))
for (line, [wave, flux]) in dict_lines3.items():
    ax3.annotate(line, xy=(wave, flux * 1), xytext=(wave - 30, 3.0), arrowprops=dict_arrow, rotation='vertical')

plt.xlabel('Rest Wavelength ($\AA$)', fontsize=16)
fig.savefig('OUTPUT_PlotMultipleSpec.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
#  Plots The Spectra In Comparison With Spectra Of Other Well-Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def plot_spec(ax_obj, file_name, phase='0d', offset=0):
    if file_name[0:3] != 'rfz':
        name = file_name.split('/')[-1].split('_')[0]
        phase = file_name.split('/')[-1].split('_')[-1].split('.')[0]
    else:
        name = file_name.split('.')[0].split('_')[-1].split('-')[0]

    data_df = pd.read_csv(file_name, sep='\s+', engine='python', header=None, comment='#')
    data_df = data_df[(data_df[0] > 3750) & (data_df[0] < 9150)]
    data_df[1] = data_df[1] / data_df[1].mean()

    ax_obj.plot(data_df[0], data_df[1] - offset, label=phase + ' ' + name)
    ax_obj.text(data_df[0].tolist()[-1] + 150, data_df[1].tolist()[-1] - offset, s=phase + ' ' + name, fontsize=7)

    ax_obj.set_xlim(3600, 10800)
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(1000))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(500))
    ax_obj.tick_params(axis='y', which='both', direction='in', labelleft='off')
    ax_obj.tick_params(axis='x', which='both', direction='in', width=0.5, labelsize=12)


list_spec1 = ['1999em_24d.dat', '2004et_25d.dat', '2012aw_24d.dat', '2013ej_24d.dat']
list_spec2 = ['1999em_76d.dat', '2004et_74d.dat', '2012aw_76d.dat', '2013ej_70d.dat']
list_spec3 = ['1999em_163d.dat', '2004et_185d.dat', '2012aw_250d.dat', '2013ej_183d.dat']

fig = plt.figure(figsize=(6, 9))

ax1 = fig.add_subplot(311)
plot_spec(ax1, 'rfz_jul30_cfwcbs_ASASSN14dq.dat', phase='25d')
for index, file_name in enumerate(list_spec1):
    plot_spec(ax1, DIR_SPEC + file_name, offset=index + 0.8)

ax2 = fig.add_subplot(312, sharex=ax1)
plot_spec(ax2, 'rfz_sep12_cfwcbs_ASASSN14dq.dat', phase='74d')
for index, file_name in enumerate(list_spec2):
    plot_spec(ax2, DIR_SPEC + file_name, offset=(index + 1) * 1.5)
ax2.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

ax3 = fig.add_subplot(313, sharex=ax1)
plot_spec(ax3, 'rfz_dec08_cfwcbs_ASASSN14dq.dat', phase='158d')
for index, file_name in enumerate(list_spec3):
    plot_spec(ax3, DIR_SPEC + file_name, offset=(index + 1) * 4.5)

ax3.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
fig.subplots_adjust(hspace=0.0, top=0.9, right=0.95)
fig.savefig('OUTPUT_PlotSpecComp.eps', format='eps', dpi=500, bbox_inches='tight')
plt.show()
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
            if 8200 > row[0] >= 7500:
                row[1] = row[1] * 0.75
            elif row[0] >= 8200:
                row[1] = row[1] * 0.45

    ax_obj.plot(data_df[0], np.log10(data_df[1]) + offset, linewidth=1, c=color, alpha=0.8, label=label + ' ' + phase)


fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

# plot_modspec(ax, 'rfz_aug31_cfwcbs_ASASSN14dq.dat', phase='59.7 d', offset=0.05)
# plot_modspec(ax, DIR_MODEL + 'm15z2m3/m15_du_sch_z2m3_FeC_mix0p4_20.fl', phase='61.3 d', offset=-0.5)
# plot_modspec(ax, DIR_MODEL + 'm15z8m3/m15_du_sch_z8m3_FeC_mix0p4_20.fl', phase='60.2 d', offset=-1.0)
# plot_modspec(ax, DIR_MODEL + 'm15z2m2/m15_du_sch_FeC_mix0p4_20.fl', phase='59.7 d', offset=-1.5)
# plot_modspec(ax, DIR_MODEL + 'm15z4m2/m15_du_sch_z4m2_FeC_mix0p4_20.fl', phase='60.6 d', offset=-2.0)

plot_modspec(ax, 'rfz_aug08_cfwcbs_ASASSN14dq.dat', '36.7 d', offset=0.3)
plot_modspec(ax, DIR_MODEL + 'm15z2m3/m15_du_sch_z2m3_FeC_mix0p4_15.fl', '38.1 d', -0.5, z=r'0.1 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z8m3/m15_du_sch_z8m3_FeC_mix0p4_15.fl', '37.4 d', -0.9, z=r'0.4 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z2m2/m15_du_sch_FeC_mix0p4_15.fl', '37.1 d', -1.3, z=r'1 $\rm Z_{\odot}$')
plot_modspec(ax, DIR_MODEL + 'm15z4m2/m15_du_sch_z4m2_FeC_mix0p4_15.fl', '37.6 d', -1.7, z=r'2 $\rm Z_{\odot}$')

ax.set_xlim(3700, 9500)
ax.set_ylim(-1.9, 0.7)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(axis='y', which='both', direction='in', labelleft='off')
ax.tick_params(axis='x', which='both', direction='in', width=1, labelsize=14)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
ax.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=15)

fig.savefig('OUTPUT_PlotCompSpecModel.eps', format='eps', dpi=1200, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
