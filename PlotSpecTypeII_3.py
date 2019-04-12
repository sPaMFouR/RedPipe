# !/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx-----------------------PLOT SPECTRA OF SUPERNOVA----------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
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
import seaborn as sns
from pyraf import iraf
from jdcal import gcal2jd
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

sns.set_style('ticks')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
JD_keyword = 'JD'
lower_lim = 3600
upper_lim = 9050
light_speed = 2.99792458e5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Details Of SN In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
name_hostgal = 'NGC2276'
redshift = 0.008059
phase_early = 30
phase_nebstart = 115
date_explosion = 2457641.4
dist_val = 29.64
dist_err = 2.65
mNi = 0.033
tc = 486
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
DIR_GAL = "/home/avinash/Supernovae_Data/2016gfy/NGC2276/HCTSpec/"
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
            flux_array = hdulist[0].data * u.Angstrom
            wave_array = spec.read_fits_spectrum1d(file_name).dispersion * u.Unit('erg cm-2 s-1 angstrom-1')
        else:
            flux_array = hdulist[0].data[0][0] * u.Angstrom
            wave_array = spec.read_fits_spectrum1d(file_name)[0].dispersion * u.Unit('erg cm-2 s-1 angstrom-1')

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


def smooth_1dspec(common_text='', textlist_spectra='', sp=1.5, kernel='gaussian', prefix_str='z_', plot=False):
    """
    Smoothens a 1-D spectra based on the smoothening parameter. Smoothening parameter
    is 'std.dev.' in case of isotropic Gaussian filter and is 'width' in the case of the
    non-isotropic box filter.
    Args:
        common_text         : Common text of 1-D spectra files which have to be smoothened
        textlist_spectra    : Text list of 1-D spectra FITS files
        sp                  : Smoothening parameter
        kernel              : Convolution Kernel used for smoothening (Gaussian or Box)
        prefix_str          : Prefix to distinguish the smoothened 1-D spectra from the original
        plot                : Boolean describing whether the smoothened spectra has to be plotted
    Returns:
        None
    """
    if textlist_spectra == '':
        list_spectra = group_similar_files('', common_text=common_text)
    elif common_text == '':
        list_spectra = text_list_to_python_list(textlist_spectra)
    else:
        print ("Error: Specify Only 1 Option For Selection Of Files")
        sys.exit(1)

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


def convert_to_text(common_text='', textlist_spectra='', out_ext='.dat'):
    """
    Converts a list of 1-D spectra FITS files into a two-column text file.
    Args:
        common_text      : Common text of 1-D spectra FITS files
        textlist_spectra : Text list of 1-D spectra FITS files
        out_ext          : Name of the extension to be used for the output file
    Returns:
        list_outspectra  : Python list of 1-D spectra DAT files
    """
    if textlist_spectra == '':
        list_spectra = group_similar_files('', common_text=common_text)
    elif common_text == '':
        list_spectra = text_list_to_python_list(textlist_spectra)
    else:
        print ("Error: Specify Only 1 Option For Selection Of Files")
        sys.exit(1)

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
# Function To Plot Individual Spectrum
# ------------------------------------------------------------------------------------------------------------------- #

def plot_velevol(ax_obj, file_df, index, wavelength=6563, offset=4.25e-16, smooth=False, sp=1, plot_label=True):
    """
    Plots the balmer feature evolution with labelled epochs and sets plot parameters.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_df     : Pandas DataFrame containing 1D Spectrum data to be plotted
        index       : Index of the file 'file_name' which is to be plotted
        offset      : Offset to be applied to the 1D Spectrum to be plotted
        wavelength  : Central wavelength of the balmer feature to be plotted
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        plot_label  : Boolean depicting whether the label associated with the 1D Spectrum is to be plotted
    Returns:
        None
    """
    if wavelength == 6563:
        file_df.loc[:, 'Flux'] += 0.81e-15

    if smooth:
        file_df.loc[:, 'Flux'] = convolve(file_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    file_df['Velocity'] = file_df['Wavelength'].apply(lambda x: ((x - wavelength) / wavelength) * light_speed / 1e3)
    file_df = file_df[(file_df['Velocity'] > -18) & (file_df['Velocity'] < 6)]
    file_df.loc[:, 'Flux'] -= offset * index

    if plot_label:
        ax_obj.text(x=file_df['Velocity'].values[-1] + 1.2, y=file_df['Flux'].values[-1],
                    s=evolution_df.loc[index, 'Label'], fontsize=14)

    ax_obj.plot(file_df['Velocity'], file_df['Flux'])


def plot_epoch(ax_obj, file_name, index, master_df, offset=0.8e-15, smooth=False, sp=2):
    """
    Plots the spectrum with line identified and labelled accordingly.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        index       : Index of the file 'file_name' which is to be plotted
        master_df   : Master Pandas DataFrame from which the label is to be read
        offset      : Offset in Flux units to be applied to the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
    Returns:
        None
    """
    data_df = pd.read_csv(file_name, names=['Wave', 'Flux'], sep='\s+', dtype='float64')
    data_df['Flux'] -= offset * index

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    if data_df['Wave'].max() >= int(upper_lim):
        data_df = data_df[(data_df['Wave'] >= int(lower_lim)) & (data_df['Wave'] <= int(upper_lim))]
    else:
        data_df = data_df[(data_df['Wave'] >= int(lower_lim)) & (data_df['Wave'] <= data_df['Wave'].max() - 30)]

    ax_obj.plot(data_df['Wave'], data_df['Flux'], linewidth=1, label=None)
    ax_obj.text(x=data_df['Wave'].values[-1] + 50, y=data_df['Flux'].values[-1],
                s=master_df.loc[index, 'Label'], fontsize=10)


def plot_galacticlines(ax_obj, lines=[4861, 4959, 5007, 6563, 6717, 6731]):
    """
    Plots host galactic lines in the plot denoted by axes 'ax_obj' with dashed lines.
    Args:
        ax_obj  : Axes object to be used for plotting
        lines   : Python list of all the galactic lines to be marked
    Returns:
        None
    """
    for gline in lines:
        ax.axvline(gline, color='k', linestyle='--', linewidth=0.8, alpha=0.6, label='_nolegend_')


def set_plotparams(ax_obj):
    """
    Sets plot parameters to the axes object 'ax_obj'.
    Args:
        ax_obj  : Axes object to be used for plotting and setting plot parameters
    Returns:
        None
    """
    ax_obj.set_yticklabels([])
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
    ax_obj.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
    ax_obj.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
    ax_obj.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files',
#                         choices=['Yes', 'No'])
# ctext = eg.enterbox('Common Text Of Files To Be Flux Calibrated?', title='Flux Calibration', default='fz_*.fits')
# bool_smooth = eg.boolbox('Perform Smoothening Of Spectra?', title='Smoothening 1-D Spectra', choices=['Yes', 'No'])
rmv_files = True
ctext = 'fz_*.fits'
bool_smooth = False
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
    smooth_1dspec(common_text='r' + ctext, sp=3, kernel='gaussian', prefix_str='z_', plot=False)
    list_textspec = convert_to_text(common_text='z_r' + ctext)
else:
    list_textspec = convert_to_text(common_text='r' + ctext)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Creates A Pandas DataFrame With Information Of All 1-Dimensional Spectra
# ------------------------------------------------------------------------------------------------------------------- #
spec_df = pd.DataFrame(list_spec, columns=['FileName'])
spec_df['TextFile'] = list_textspec
spec_df['Phase'] = spec_df['FileName'].apply(lambda x: read_jd(x) - date_explosion)
spec_df['Label'] = spec_df['Phase'].apply(lambda x: '+{0:>.1f} d'.format(x) if x > 0 else '{0:>.1f} d'.format(x))
spec_df = spec_df.sort_values(by='Phase')
spec_df = spec_df.reset_index(drop=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Copy Pandas DataFrame Into Subset DataFrames Required For Various Plots
# ------------------------------------------------------------------------------------------------------------------- #
early_df = spec_df[spec_df['Phase'] < phase_early].copy()
late_df = spec_df[(spec_df['Phase'] > phase_early) & (spec_df['Phase'] < phase_nebstart)].copy()
plateau_df = spec_df[spec_df['Phase'] < phase_nebstart].copy()
nebular_df = spec_df[spec_df['Phase'] > phase_nebstart].copy()
evolution_df = spec_df[(spec_df['Phase'] > 22) & (spec_df['Phase'] < phase_nebstart)].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plots The 1-Dimensional Spectra During The Plateau Phase
# # ------------------------------------------------------------------------------------------------------------------- #
# sns.set_palette(sns.color_palette('rocket', 4))

# fig = plt.figure(figsize=(12, 12))
# ax = fig.add_subplot(111)

# for index, file_name in plateau_df['TextFile'].items():
#     plot_epoch(ax, file_name, index, plateau_df, offset=0.8e-15)

# dict_labelsplat = {r'$\rm Ca\,II$ (H & K)': [3830, 60], r'$\rm H\delta$ 4102': [3990, 40],
#                    r'$\rm H\gamma$ 4340': [4210, 70], r'$\rm Ba\,II\ 4554\ +\ Fe\,II$': [4490, 50],
#                    r'$\rm H\beta$ 4861': [4730, 70], r'$\rm Fe\,II$ Lines': [5000, 140], r'$\rm Na\,I$ D': [5810, 50],
#                    r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5460, 50], r'$\rm Sc\,II$ 5663 Multiplet': [5610, 50],
#                    r'$\rm H\alpha$ 6563': [6380, 100], r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40],
#                    r'$\rm \oplus[O_2$ 7620 A-band]': [7560, 40], r'$\rm O\,I$ 7774': [7700, 50],
#                    r'$\rm Ca\,II$ NIR Triplet': [8480, 180]}

# ax.set_ylim(-15.5e-15, 6.3e-15)
# ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

# for (line, [wavelength, width]) in dict_labelsplat.items():
#     ax.text(wavelength - 30, 5.6e-15, line, rotation='vertical', fontsize=10)
#     ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

# set_plotparams(ax)
# plot_galacticlines(ax)
# ax.xaxis.set_major_locator(MultipleLocator(500))
# ax.xaxis.set_minor_locator(MultipleLocator(50))
# ax.yaxis.set_major_locator(MultipleLocator(4e-15))
# ax.yaxis.set_minor_locator(MultipleLocator(0.4e-15))

# fig.savefig('PLOT_SpecPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plots The 1-Dimensional Spectra During The Early Plateau Phase
# # ------------------------------------------------------------------------------------------------------------------- #
# sns.set_palette(sns.color_palette('rocket_d', 4))

# fig = plt.figure(figsize=(12, 6))
# ax = fig.add_subplot(111)

# for index, file_name in early_df['TextFile'].items():
#     plot_epoch(ax, file_name, index, early_df, offset=0.7e-15)

# dict_labelsearly = {r'$\rm Ca\,II$ (H & K)': [3800, 60], r'$\rm H\delta$ 4102': [3960, 50],
#                     r'$\rm H\gamma$ 4340': [4190, 50], r'$\rm H\beta$ 4861': [4700, 70],
#                     r'$\rm Fe\,I$ 5169': [5070, 50],
#                     r'$\rm He\,I$ 5876': [5660, 70], r'$\rm H\alpha$ 6563': [6360, 60],
#                     r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40], r'$\rm \oplus[H_{2}O$ 7165]': [7115, 40],
#                     r'$\rm \oplus[O_2$ 7620 A-band]': [7560, 50], r'$\rm \oplus[H_{2}O$ band]': [8100, 40], }

# ax.set_ylim(-3.4e-15, 4.1e-15)
# ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

# for (line, [wavelength, width]) in dict_labelsearly.items():
#     ax.text(wavelength - 30, 3.7e-15, line, fontsize=9, rotation='vertical')
#     ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

# set_plotparams(ax)
# plot_galacticlines(ax)
# ax.xaxis.set_major_locator(MultipleLocator(500))
# ax.xaxis.set_minor_locator(MultipleLocator(50))
# ax.yaxis.set_major_locator(MultipleLocator(2e-15))
# ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))

# fig.savefig('PLOT_SpecEarlyPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plots The 1-Dimensional Spectra During The Late Plateau Phase
# # ------------------------------------------------------------------------------------------------------------------- #
# sns.set_palette(sns.color_palette('rocket', 4))

# fig = plt.figure(figsize=(12, 10))
# ax = fig.add_subplot(111)

# for index, file_name in late_df['TextFile'].items():
#     plot_epoch(ax, file_name, index, late_df, offset=0.9e-15)

# dict_labelslate = {r'$\rm Ca\,II$ (H & K)': [3840, 60], r'$\rm H\gamma$ 4340': [4190, 50],
#                    r'$\rm Ba\,II\ 4554\ +\ Fe\,II$': [4490, 50], r'$\rm H\beta$ 4861': [4750, 60],
#                    r'$\rm Fe\,II$ Triplet': [5010, 130], r'$\rm Fe\,II$ 5267, 5363 Blend': [5240, 50],
#                    r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5480, 40], r'$\rm Sc\,II$ 5663 Multiplet': [5610, 40],
#                    r'$\rm Na\,ID$ 5890, 5896': [5800, 50], r'$\rm Ba\,I$ 6142': [6090, 30],
#                    r'$\rm Sc\,II$ 6246': [6200, 30], r'$\rm H\alpha$ 6563': [6410, 70],
#                    r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40], r'$\rm \oplus[O_2$ 7620 A-band]': [7550, 40],
#                    r'$\rm O\,I$ 7774': [7690, 40], r'$\rm Ca\,II$ 8498, 8542': [8370, 60],
#                    r'$\rm Ca\,II$ 8662': [8530, 60]}

# ax.set_ylim(-17.4e-15, -0.5e-15)
# ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

# for (line, [wavelength, width]) in dict_labelslate.items():
#     ax.text(wavelength - 30, -1.0e-15, line, fontsize=9, rotation='vertical')
#     ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

# set_plotparams(ax)
# plot_galacticlines(ax)
# ax.xaxis.set_major_locator(MultipleLocator(500))
# ax.xaxis.set_minor_locator(MultipleLocator(50))
# ax.yaxis.set_major_locator(MultipleLocator(2e-15))
# ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))

# fig.savefig('PLOT_SpecLatePlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Plots The 1-Dimensional Spectra During The Nebular Phase
# # ------------------------------------------------------------------------------------------------------------------- #
# sns.set_palette(sns.color_palette('rocket_d', 4))

# fig = plt.figure(figsize=(12, 8))
# ax = fig.add_subplot(111)

# for index, file_name in nebular_df['TextFile'].items():
#     plot_epoch(ax, file_name, index, nebular_df, offset=2.3e-16, smooth=True)

# dict_labelsneb = {r'$\rm Na\,ID$ 5890, 5896': [5870, 60], r'$\rm H\alpha$ 6563': [6560, 80],
#                   r'$\rm [O\,I]$ 6300, 6364': [6300, 60], r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40],
#                   r'$\rm [Ca\,II]$ 7291, 7324': [7290, 50], r'$\rm \oplus[O_2$ 7620 A-band]': [7550, 40],
#                   r'$\rm OI$ 7774': [7720, 50], r'$\rm Ca\,II$ 8498, 8542': [8510, 50],
#                   r'$\rm Ca\,II$ 8662': [8670, 50]}

# ax.set_ylim(-7.5e-15, -3.3e-15)
# ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

# for (line, [wavelength, width]) in dict_labelsneb.items():
#     ax.text(wavelength - 30, -3.4e-15, line, rotation='vertical', fontsize=9)
#     ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

# set_plotparams(ax)
# plot_galacticlines(ax)
# ax.xaxis.set_major_locator(MultipleLocator(500))
# ax.xaxis.set_minor_locator(MultipleLocator(50))
# ax.yaxis.set_major_locator(MultipleLocator(0.5e-15))
# ax.yaxis.set_minor_locator(MultipleLocator(0.05e-15))

# fig.savefig('PLOT_SpecNebular.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Balmer Features Evolution Feature Across Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
# sns.set_palette(sns.color_palette('Paired', 10))
# fig, (ax_hbeta, ax_halpha) = plt.subplots(1, 2, figsize=(9, 16), sharey=True)

# for index, file_name in evolution_df['TextFile'].items():
#     data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
#     data_df = data_df[data_df['Wavelength'] <= 6800]

#     plot_velevol(ax_hbeta, data_df, index, wavelength=4861, plot_label=False, smooth=True)
#     plot_velevol(ax_halpha, data_df, index, offset=4.73e-16, smooth=True)

# dict_halpha = {(-9.9, -9.0): [-1.45e-15, -4.4e-15, '--'], (-9.0, -8.9): [-4.4e-15, -4.96e-15, '--'],
#                (-8.9, -8.7): [-4.96e-15, -7.91e-15, '--'], (-8.7, -7.8): [-1.5e-15, -3.15e-15, '-'],
#                (-7.8, -7.0): [-3.15e-15, -3.95e-15, '-'], (-7.0, -6): [-3.95e-15, -7.98e-15, '-']}

# dict_hbeta = {(-9.8, -9.0): [-1.8e-15, -7.78e-15, '--'], (-7.4, -6.9): [-1.95e-15, -2.3e-15, '-'],
#               (-6.9, -6.6): [-2.30e-15, -3.0e-15, '-'], (-6.6, -6.0): [-3.0e-15, -3.57e-15, '-'],
#               (-6.0, -4.8): [-3.57e-15, -7.90e-15, '-']}

# for ((xstart, xend), [ystart, yend, ls]) in dict_halpha.items():
#     ax_halpha.plot([xstart, xend], [ystart, yend], linewidth=3, alpha=0.4, linestyle=ls, color='r')
# for ((xstart, xend), [ystart, yend, ls]) in dict_hbeta.items():
#     ax_hbeta.plot([xstart, xend], [ystart, yend], linewidth=3, alpha=0.4, linestyle=ls, color='b')

# ax_halpha.text(x=-15, y=0.2e-15, s=r'$\rm H_{\alpha}\ 6563\ \AA$', fontsize=18)
# ax_hbeta.text(x=-15, y=0.2e-15, s=r'$\rm H_{\beta}\ 4861\ \AA$', fontsize=18)
# ax_halpha.axvline(0, linestyle='--', color='k', linewidth=0.8)
# ax_hbeta.axvline(0, linestyle='--', color='k', linewidth=0.8)

# ax_halpha.set_xlim(-16, 7)
# ax_halpha.set_ylim(-8.1e-15, 1.0e-15)
# ax_halpha.set_yticklabels([])
# ax_halpha.xaxis.set_ticks_position('both')
# ax_halpha.yaxis.set_ticks_position('both')
# ax_halpha.xaxis.set_major_locator(MultipleLocator(5))
# ax_halpha.xaxis.set_minor_locator(MultipleLocator(1))
# ax_halpha.yaxis.set_major_locator(MultipleLocator(1e-15))
# ax_halpha.yaxis.set_minor_locator(MultipleLocator(1e-16))
# ax_halpha.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=15)
# ax_halpha.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=15)

# ax_hbeta.set_xlim(-16, 7)
# ax_hbeta.set_yticklabels([])
# ax_hbeta.xaxis.set_ticks_position('both')
# ax_hbeta.yaxis.set_ticks_position('both')
# ax_hbeta.xaxis.set_major_locator(MultipleLocator(5))
# ax_hbeta.xaxis.set_minor_locator(MultipleLocator(1))
# ax_hbeta.yaxis.set_major_locator(MultipleLocator(1e-15))
# ax_hbeta.yaxis.set_minor_locator(MultipleLocator(1e-16))
# ax_hbeta.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=15)
# ax_hbeta.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=15)

# ax_halpha.set_xlabel(r'Velocity [$\rm \times 10^3\ km\ s^{-1}$]', fontsize=18)
# ax_hbeta.set_xlabel(r'Velocity [$\rm \times 10^3\ km\ s^{-1}$]', fontsize=18)
# ax_hbeta.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=18)

# fig.subplots_adjust(wspace=0.02)
# fig.savefig('PLOT_BalmerEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
# plt.show()
# plt.close(fig)
# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Halpha Early Evolution
# ------------------------------------------------------------------------------------------------------------------- #

def plot_hacomp(ax_obj, path_name, phase='0d', offset=0, smooth=False, sp=1, legend=True, clip='6000:7000', ls='-'):
    """
    Plots the line evolution of line the specified by the wavelength region in 'clip_str'.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        phase       : Phase of the spectrum to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        legend      : Whether the legend will be used to label the plots
        clip        : Wavelength region to be used for plotting
        ls          : Linestyle to be used for plotting
    Returns:
        None
    """
    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'] / np.median(data_df['Flux'])

    lower_clip, upper_clip = clip.split(':')
    file_name = path_name.split('/')[-1]

    if file_name[0:3] != 'rfz':
        name = file_name.split('_')[0]
        phase = file_name.split('_')[-1].rstrip('d') + ' d'
        if name == '2016esw':
            offset = -0.1
    else:
        name = name_SN

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]
    data_df['Vel'] = data_df['Wave'].apply(lambda x: (x / 6562.63 - 1) * light_speed * 1e-3)
    data_df['Flux'] = data_df['Flux'] - offset

    if offset == 0:
        ax_obj.plot(data_df['Vel'], np.log10(data_df['Flux']), lw=1.2, ls=ls, label=' +' + phase + '')
    else:
        ax_obj.plot(data_df['Vel'], np.log10(data_df['Flux']), lw=1.2, ls=ls, label=name + ' [+' + phase + ']')

    if not legend:
        ax_obj.text(-27, np.log10(data_df['Flux']).tolist()[0] + 0.03, fontsize=12, s=name + ' [+' + phase + ']')


list_compHalpha = ['BoxyHalpha/2007od_5.5d', 'BoxyHalpha/2007od_9.2d', 'BoxyHalpha/2016esw_19.5d']
list_compCaII = ['BoxyHalpha/2007od_5.5d' + 'BoxyHalpha/2007od_5.5d', 'BoxyHalpha/2007od_5.5d']

# sns.set_palette(sns.color_palette('colorblind', 10))
sns.set_palette(sns.color_palette('rocket', 5))
fig_ha, (ax_ha, ax_comp) = plt.subplots(2, 1, figsize=(8, 14), gridspec_kw={'height_ratios': [2, 3]}, sharex=True)

plot_hacomp(ax_ha, DIR_SPEC + 'rfz_2016-09-20_2016gfy.dat', phase='11 d', smooth=True, sp=1)
plot_hacomp(ax_ha, DIR_SPEC + 'rfz_2016-09-27_2016gfy.dat', phase='18 d', smooth=True, sp=1, ls='--')
plot_hacomp(ax_ha, DIR_SPEC + 'rfz_2016-09-30_2016gfy.dat', phase='21 d', smooth=True, sp=1, ls=':')

ax_ha.set_yticklabels([])
ax_ha.set_ylim(-0.18, 0.44)
ax_ha.text(-25, 0.38, s=r'SN 2016gfy', fontsize=14)
ax_ha.axvline(0, lw=1, ls='--', color='k')

ax_ha.legend(markerscale=3, frameon=False, fontsize=16)
ax_ha.yaxis.set_ticks_position('both')
ax_ha.xaxis.set_ticks_position('both')
ax_ha.xaxis.set_major_locator(MultipleLocator(10))
ax_ha.xaxis.set_minor_locator(MultipleLocator(1))
ax_ha.yaxis.set_major_locator(MultipleLocator(0.1))
ax_ha.yaxis.set_minor_locator(MultipleLocator(0.01))
ax_ha.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=16)
ax_ha.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=16)
ax_ha.set_ylabel(r'Log $\rm\ F_{\lambda}$', fontsize=18)

plot_hacomp(ax_comp, DIR_SPEC + 'rfz_2016-09-20_2016gfy.dat', phase='11 d', offset=0.1, legend=False, smooth=False)
for index, path_name in enumerate(list_compHalpha):
    plot_hacomp(ax_comp, DIR_SNe + 'Spec_Data/' + path_name, offset=index * 0.2 + 0.2, legend=False, smooth=False)

ax_comp.set_ylim(-0.66, 0.4)
ax_comp.set_yticklabels([])
ax_comp.axvline(0, lw=1, ls='--', color='k')

ax_comp.yaxis.set_ticks_position('both')
ax_comp.xaxis.set_ticks_position('both')
ax_comp.xaxis.set_major_locator(MultipleLocator(10))
ax_comp.xaxis.set_minor_locator(MultipleLocator(1))
ax_comp.yaxis.set_major_locator(MultipleLocator(0.2))
ax_comp.yaxis.set_minor_locator(MultipleLocator(0.02))
ax_comp.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=16)
ax_comp.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=16)
ax_comp.set_xlabel(r'Velocity $\rm [\times 10^3\ km\ s^{-1}$]', fontsize=18)
ax_comp.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=18)

fig_ha.subplots_adjust(hspace=0.01)
fig_ha.savefig('PLOT_SpecBoxyHalpha.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_ha)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Individual Spectra Of Galaxy
# Plots The Host Galaxy Spectrum [Nucleus  + Host HII Region]
# ------------------------------------------------------------------------------------------------------------------- #

def plot_galspec(ax_obj, path_name, offset=0, smooth=False, sp=1, clip_str='3900:8700'):
    """
    Plots the spectrum of the nucleus of the host galaxy and the host HII region.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        path_name   : Name of the 1-D Spectrum FITS file to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        clip_str    : Wavelength region to be used for plotting
    Returns:
        None
    """
    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    ax_obj.plot(data_df['Wave'], np.log10(data_df['Flux']) + offset, lw=1.5, c='dimgrey',
                alpha=0.5, label='_nolegend_')

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    lower_clip, upper_clip = clip_str.split(':')
    if int(upper_clip) > data_df['Wave'].max():
        data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] <= data_df['Wave'].max() - 50)]
    else:
        data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]

    ax_obj.plot(data_df['Wave'], np.log10(data_df['Flux']) + offset, lw=1.2, c='k', alpha=0.7, label='_nolegend_')


list_specgal = group_similar_files('', common_text=DIR_GAL + 'Iter6/z_*' + name_hostgal + '*.dat')

fig_host = plt.figure(figsize=(12, 8))
ax_host = fig_host.add_subplot(111)

for index, path_name in enumerate(list_specgal):
    plot_galspec(ax_host, path_name, offset=index * -1.0, smooth=True, sp=2)

# ax_host.axvline(4861, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(4959, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(5007, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6548, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6563, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6584, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6717, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6731, color='k', ls='-.', lw=0.8, alpha=0.4, label='_nolegend_')

ax_host.text(4370, -15.25, color='k', s=r'$\rm H\gamma$', fontsize=13)
ax_host.text(4740, -15.2, color='k', s=r'$\rm H\beta$', fontsize=13)
ax_host.text(5040, -15.0, color='k', s=r'$\rm [O\,III]$', fontsize=13)
ax_host.text(6200, -14.9, color='k', s=r'$\rm H\alpha\,+\,[N\,II]$', fontsize=13)
ax_host.text(6750, -15.25, color='k', s=r'$\rm [S\,II]$', fontsize=13)
ax_host.text(3950, -15.15, color='r', s=r'HII Region', fontsize=15)
ax_host.text(3950, -16.05, color='b', s=r'Nucleus', fontsize=15)

ax_host.set_ylim(-16.5, -14.55)
ax_host.set_xlim(3820, 7150)
ax_host.set_yticklabels([])
ax_host.yaxis.set_ticks_position('both')
ax_host.xaxis.set_ticks_position('both')
ax_host.xaxis.set_major_locator(MultipleLocator(500))
ax_host.xaxis.set_minor_locator(MultipleLocator(50))
ax_host.yaxis.set_major_locator(MultipleLocator(0.2))
ax_host.yaxis.set_minor_locator(MultipleLocator(0.04))
ax_host.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=18)
ax_host.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=18)

ax_host.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=18)
ax_host.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=18)

fig_host.savefig('PLOT_SpecHostGal.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_host)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
#  Plots The Spectra In Comparison With Spectra Of Other Well-Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def plot_compspec(ax_obj, path_name, phase='0d', offset=0, smooth=False, phaseoffset=0.05, sp=2, clip_str='3600:9000'):
    """
    Plots the spectrum for comparison with other well-studied SNe.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        phase       : Phase of the spectrum to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        phaseoffset : Offset in Flux units to be applied to the label for the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        clip_str    : Wavelength region to be used for plotting
    Returns:
        None
    """
    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'].apply(lambda x: x * 1e15 if x > 0 else np.nan)
    data_df = data_df.dropna(how='any')

    lower_clip, upper_clip = clip_str.split(':')
    file_name = path_name.split('/')[-1]

    if file_name[0:3] != 'rfz':
        name = file_name.split('_')[0][2:]
        phase = file_name.split('_')[-1].rstrip('d') + ' d'
    else:
        name = name_SN[2:]
#         copy_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]
#         ax_obj.plot(copy_df['Wave'], np.log10(copy_df['Flux']), lw=1, alpha=0.5, color='dimgrey', label='_nolegend')

    if float(phase.rstrip('d')) > 80:
        wave_offset = 4870
    else:
        wave_offset = 5540

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    if name in dict_redshift.keys():
        data_df['Wave'] = data_df['Wave'].apply(lambda x: (1 - dict_redshift[name]) * x)
        upper_clip = 8000

    data_df['Flux'] = np.log10(data_df['Flux']) - offset
    data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]

    ax_obj.plot(data_df['Wave'], data_df['Flux'], lw=1.2, label=phase + ' ' + name)
    ax_obj.text(s=name + ' [+' + phase + ']', color='k', fontsize=14, x=wave_offset,
                y=np.mean(data_df[(data_df['Wave'] > 5000) & (data_df['Wave'] < 6300)]['Flux']) + phaseoffset)

    ax_obj.set_yticklabels([])
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(1000))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(100))
    ax_obj.tick_params(which='major', direction='in', width=1.6, length=8, labelsize=18)
    ax_obj.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=18)


sns.set_palette(sns.cubehelix_palette(6, start=1.2, rot=8, dark=0.2, light=.9, reverse=True))

dict_redshift = {'16esw': 0.02831}
list_specmax = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/NearExp/*')
list_specearlyplat = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/EarlyPlat/*')
list_speclateplat = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/LatePlat/*')
list_specnebular = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/Nebular/*')

fig1 = plt.figure(figsize=(12, 16))

ax1 = fig1.add_subplot(211)
dict_offset1 = {0: 0.2, 1: 0.4, 2: 1.3, 3: 13.2}
plot_compspec(ax1, DIR_SPEC + 'rfz_2016-09-13_2016gfy.dat', phase='04 d', smooth=True, sp=4)
for index, path_name in enumerate(list_specmax):
    plot_compspec(ax1, path_name, offset=dict_offset1[index], phaseoffset=0.06, smooth=False)

ax1.set_ylim(-0.65, 1.7)
ax1.yaxis.set_major_locator(MultipleLocator(0.4))
ax1.yaxis.set_minor_locator(MultipleLocator(0.04))
ax1.text(6785, 1.05, s=r'$\rm \oplus$', fontsize=14)
ax1.text(7520, 1.05, s=r'$\rm \oplus$', fontsize=14)
ax1.text(7700, 1.4, s='Near Maximum', color='r', fontsize=20)

ax2 = fig1.add_subplot(212, sharex=ax1)
dict_offset2 = {0: 0.25, 1: 0.3, 2: 0.2}
plot_compspec(ax2, DIR_SPEC + 'rfz_2016-10-04_2016gfy.dat', phase='25 d', phaseoffset=0.05, smooth=True, sp=4)
for index, path_name in enumerate(list_specearlyplat):
    plot_compspec(ax2, path_name, offset=dict_offset2[index], phaseoffset=.05, smooth=False)

ax2.set_ylim(-0.45, 1.5)
ax2.set_xlim(3400, 9200)
ax2.yaxis.set_major_locator(MultipleLocator(0.4))
ax2.yaxis.set_minor_locator(MultipleLocator(0.04))
ax2.text(6795, 1.2, s=r'$\rm \oplus$', fontsize=14)
ax2.text(7520, 1.1, s=r'$\rm \oplus$', fontsize=14)
ax2.text(7400, 1.3, s='Early Plateau Phase', color='r', fontsize=20)

ax1.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=20)
ax2.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=20)
ax2.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=20)

fig1.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig1.savefig('PLOT_SpecComp1.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig1)

fig2, (ax3, ax4) = plt.subplots(2, 1, figsize=(12, 16), gridspec_kw={'height_ratios': [4, 5]}, sharex=True)

dict_offset3 = {0: 12.85, 1: -0.5, 2: 0.0}
plot_compspec(ax3, DIR_SPEC + 'rfz_2016-11-24_2016gfy.dat', phase='76 d', phaseoffset=0.15,
              smooth=True, sp=2.2, clip_str='3700:9000')
for index, path_name in enumerate(list_speclateplat):
    plot_compspec(ax3, path_name, offset=dict_offset3[index], phaseoffset=0.3, smooth=True, clip_str='3700:9000')

dict_offset4 = {0: -2, 1: -1.4, 2: -0.8}
plot_compspec(ax4, DIR_SPEC + 'rfz_2017-04-13_2016gfy.dat', phase='216 d', phaseoffset=-0.45,
              smooth=True, sp=3, clip_str='3800:9000')
for index, path_name in enumerate(list_specnebular):
    plot_compspec(ax4, path_name, offset=dict_offset4[index], phaseoffset=-0.45, smooth=True, clip_str='3800:9000')

ax3.set_ylim(-0.9, 3.1)
ax3.yaxis.set_major_locator(MultipleLocator(0.8))
ax3.yaxis.set_minor_locator(MultipleLocator(0.08))
ax3.text(6795, 1.6, s=r'$\rm \oplus$', fontsize=14)
ax3.text(7520, 1.6, s=r'$\rm \oplus$', fontsize=14)
ax3.text(7400, 2.7, s='Late Plateau Phase', color='r', fontsize=20)

ax4.set_ylim(-2.5, 3.8)
ax4.set_xlim(3500, 9200)
ax4.yaxis.set_major_locator(MultipleLocator(1))
ax4.yaxis.set_minor_locator(MultipleLocator(0.1))
ax4.text(6795, 2.6, s=r'$\rm \oplus$', fontsize=14)
ax4.text(7520, 2.4, s=r'$\rm \oplus$', fontsize=14)
ax4.text(7700, 3.3, s='Nebular Phase', color='r', fontsize=20)

ax3.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=20)
ax4.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=20)
ax4.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=20)

fig2.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig2.savefig('PLOT_SpecComp2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Spectra In Comparison With Models From Jerkstrand
# ------------------------------------------------------------------------------------------------------------------- #

def plot_jerkspecha(ax_obj, file_name, phase='0d', offset=0, color='k', smooth=False, sp=2, mass=12, log=True):
    """
    Plots the spectrum of the object along with model spectra from Dessart et al.(2013).
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        phase       : Phase of the spectrum to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        color       : Color of the spectrum to be plotted
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        mass        : String to specify the progenitor mass of the model spectrum
    Returns:
        None
    """
    data_df = pd.read_csv(file_name, sep='\s+', names=['Wavelength', 'Flux'], header=None, comment='#')
    masslabel = str(mass) + r' $\rm M_{\odot}$'

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    if file_name[0:3] != 'rfz':
        label = masslabel
        ls = '--'
        lw = 1.2
        redshift = 0.00133
        data_df['Wavelength'] = data_df['Wavelength'].apply(lambda x: (1 + redshift) * x)
        data_df['Flux'] = data_df['Flux'].apply(lambda x: x * (mNi / 0.062) * ((5.5 / dist_val) ** 2))
        data_df['Flux'] = data_df['Flux'].apply(lambda x: x / (1 - np.exp(-dict_tc[mass] / float(phase.split()[0]))))
    else:
        label = name_SN
        ls = '-'
        lw = 1.5
        data_df['Flux'] = data_df['Flux'].apply(lambda x: x / (1 - np.exp(-tc / float(phase.split()[0]))))

    data_df = data_df[(data_df['Wavelength'] > 6150) & (data_df['Wavelength'] < 6800)]

    if log:
        ax_obj.plot(data_df['Wavelength'], 17 + np.log10(data_df['Flux']) + offset, lw=lw, ls=ls,
                    c=color, alpha=0.9, label=label + ' [+' + phase + ']')
    else:
        ax_obj.plot(data_df['Wavelength'], data_df['Flux'] * 1e15 + offset, lw=lw, ls=ls,
                    c=color, alpha=0.9, label=label + ' [+' + phase + ']')


dict_tc = {12: 502, 15: 560, 19: 660, 25: 706}

fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(111)

plot_jerkspecha(ax, 'rfz_2017-04-13_2016gfy.dat', '216 d', log=False)
plot_jerkspecha(ax, DIR_SNe + 'ModJerk_Data/mzams12_212d.dat', '212 d', color='darkorange', mass=12, log=False)
plot_jerkspecha(ax, DIR_SNe + 'ModJerk_Data/mzams15_212d.dat', '212 d', color='r', mass=15, log=False)
plot_jerkspecha(ax, DIR_SNe + 'ModJerk_Data/mzams19_212d.dat', '212 d', color='b', mass=19, log=False)
plot_jerkspecha(ax, DIR_SNe + 'ModJerk_Data/mzams25_212d.dat', '212 d', color='g', mass=25, log=False)

ax.text(6270, 0.35, s=r'$\rm [O\,I]\ 6300,\ 6364\ \AA$', fontsize=16)
ax.text(6560, 0.4, s=r'$\rm H\alpha$', fontsize=18, rotation='vertical')

ax.set_ylim(0, 1.45)
ax.set_xlim(6150, 6750)
ax.legend(frameon=False, markerscale=3, fontsize=14)
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(100))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.3))
ax.yaxis.set_minor_locator(MultipleLocator(0.03))
ax.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=16)
ax.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=16)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=18)
ax.set_ylabel(r'Flux [$\rm \times\ 10^{-15}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$]', fontsize=18)

fig.subplots_adjust(wspace=0.01)
fig.savefig('PLOT_CompJerkSpecHalpha.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
