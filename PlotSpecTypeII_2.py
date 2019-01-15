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
redshift = 0.00806
phase_early = 30
phase_nebstart = 115
date_explosion = 2457641.4
dist_val = 31.30
dist_err = 2.36
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
DIR_GAL = "/home/avinash/Supernovae_Data/2016gfy/NGC2276/HCTSpec/Iter3/"
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

def plot_velevol(ax_obj, file_df, index, wavelength=6563, offset=4.25e-16, smooth=False, sp=2, plot_label=True):
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
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df['Flux'] -= offset * index

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    if data_df['Wavelength'].max() >= int(upper_lim):
        data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    else:
        data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) &
                          (data_df['Wavelength'] <= data_df['Wavelength'].max() - 30)]

    ax_obj.plot(data_df['Wavelength'], data_df['Flux'], linewidth=1, label=None)
    ax_obj.text(x=data_df['Wavelength'].values[-1] + 50, y=data_df['Flux'].values[-1],
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
evolution_df = spec_df[(spec_df['Phase'] > phase_early) & (spec_df['Phase'] < phase_nebstart)].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket', 4))

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)

for index, file_name in plateau_df['TextFile'].items():
    plot_epoch(ax, file_name, index, plateau_df, offset=0.8e-15)

dict_labelsplat = {r'$\rm Ca\,II$ (H & K)': [3830, 60], r'$\rm H\delta$ 4102': [3990, 40],
                   r'$\rm H\gamma$ 4340': [4210, 70], r'$\rm Ba\,II\ 4554\ +\ Fe\,II$': [4490, 50],
                   r'$\rm H\beta$ 4861': [4730, 70], r'$\rm Fe\,II$ Lines': [5000, 140], r'$\rm Na\,I$ D': [5810, 50],
                   r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5460, 50], r'$\rm Sc\,II$ 5663 Multiplet': [5610, 50],
                   r'$\rm H\alpha$ 6563': [6380, 100], r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40],
                   r'$\rm \oplus[O_2$ 7620 A-band]': [7560, 40], r'$\rm O\,I$ 7774': [7700, 50],
                   r'$\rm Ca\,II$ NIR Triplet': [8480, 180]}

ax.set_ylim(-15.5e-15, 6.3e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

for (line, [wavelength, width]) in dict_labelsplat.items():
    ax.text(wavelength - 30, 5.6e-15, line, rotation='vertical', fontsize=10)
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

set_plotparams(ax)
plot_galacticlines(ax)
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(4e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.4e-15))

fig.savefig('PLOT_SpecPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Early Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket_d', 4))

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(111)

for index, file_name in early_df['TextFile'].items():
    plot_epoch(ax, file_name, index, early_df, offset=0.7e-15)

dict_labelsearly = {r'$\rm Ca\,II$ (H & K)': [3800, 60], r'$\rm H\delta$ 4102': [3960, 50],
                    r'$\rm H\gamma$ 4340': [4190, 50], r'$\rm H\beta$ 4861': [4700, 70],
                    r'$\rm Fe\,I$ 5169': [5070, 50],
                    r'$\rm He\,I$ 5876': [5660, 70], r'$\rm H\alpha$ 6563': [6360, 60],
                    r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40], r'$\rm \oplus[H_{2}O$ 7165]': [7115, 40],
                    r'$\rm \oplus[O_2$ 7620 A-band]': [7560, 50], r'$\rm \oplus[H_{2}O$ band]': [8100, 40], }

ax.set_ylim(-3.4e-15, 4.1e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

for (line, [wavelength, width]) in dict_labelsearly.items():
    ax.text(wavelength - 30, 3.7e-15, line, fontsize=9, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

set_plotparams(ax)
plot_galacticlines(ax)
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(2e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))

fig.savefig('PLOT_SpecEarlyPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Late Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket', 4))

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

for index, file_name in late_df['TextFile'].items():
    plot_epoch(ax, file_name, index, late_df, offset=0.9e-15)

dict_labelslate = {r'$\rm Ca\,II$ (H & K)': [3840, 60], r'$\rm H\gamma$ 4340': [4190, 50],
                   r'$\rm Ba\,II\ 4554\ +\ Fe\,II$': [4490, 50], r'$\rm H\beta$ 4861': [4750, 60],
                   r'$\rm Fe\,II$ Triplet': [5010, 130], r'$\rm Fe\,II$ 5267, 5363 Blend': [5240, 50],
                   r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5480, 40], r'$\rm Sc\,II$ 5663 Multiplet': [5610, 40],
                   r'$\rm Na\,ID$ 5890, 5896': [5800, 50], r'$\rm Ba\,I$ 6142': [6090, 30],
                   r'$\rm Sc\,II$ 6246': [6200, 30], r'$\rm H\alpha$ 6563': [6410, 70],
                   r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40], r'$\rm \oplus[O_2$ 7620 A-band]': [7550, 40],
                   r'$\rm O\,I$ 7774': [7690, 40], r'$\rm Ca\,II$ 8498, 8542': [8370, 60],
                   r'$\rm Ca\,II$ 8662': [8530, 60]}

ax.set_ylim(-17.4e-15, -0.5e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

for (line, [wavelength, width]) in dict_labelslate.items():
    ax.text(wavelength - 30, -1.0e-15, line, fontsize=9, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

set_plotparams(ax)
plot_galacticlines(ax)
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(2e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))

fig.savefig('PLOT_SpecLatePlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Nebular Phase
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket_d', 4))

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111)

for index, file_name in nebular_df['TextFile'].items():
    plot_epoch(ax, file_name, index, nebular_df, offset=2.3e-16, smooth=True)

dict_labelsneb = {r'$\rm Na\,ID$ 5890, 5896': [5870, 60], r'$\rm H\alpha$ 6563': [6560, 80],
                  r'$\rm [O\,I]$ 6300, 6364': [6300, 60], r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40],
                  r'$\rm [Ca\,II]$ 7291, 7324': [7290, 50], r'$\rm \oplus[O_2$ 7620 A-band]': [7550, 40],
                  r'$\rm OI$ 7774': [7720, 50], r'$\rm Ca\,II$ 8498, 8542': [8510, 50],
                  r'$\rm Ca\,II$ 8662': [8670, 50]}

ax.set_ylim(-7.5e-15, -3.3e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)

for (line, [wavelength, width]) in dict_labelsneb.items():
    ax.text(wavelength - 30, -3.4e-15, line, rotation='vertical', fontsize=9)
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

set_plotparams(ax)
plot_galacticlines(ax)
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(0.5e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.05e-15))

fig.savefig('PLOT_SpecNebular.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Balmer Features Evolution Feature Across Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('Paired', 10))
fig, (ax_hbeta, ax_halpha) = plt.subplots(1, 2, figsize=(9, 16), sharey=True)

for index, file_name in evolution_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[data_df['Wavelength'] <= 6800]

    plot_velevol(ax_hbeta, data_df, index, wavelength=4861, plot_label=False, smooth=True)
    plot_velevol(ax_halpha, data_df, index, offset=4.7e-16, smooth=True)

dict_halpha = {(-9.6, -9.0): [-1.45e-15, -4.4e-15, '-.'], (-9.0, -8.9): [-4.4e-15, -4.96e-15, '-.'],
               (-8.9, -8.7): [-4.96e-15, -7.81e-15, '-.'], (-8.5, -7.8): [-1.45e-15, -3.15e-15, '-'],
               (-7.8, -7.0): [-3.15e-15, -3.95e-15, '-'], (-7.0, -6): [-3.95e-15, -7.91e-15, '-']}

dict_hbeta = {(-9.6, -9.0): [-1.8e-15, -7.78e-15, '-.'], (-6.9, -6.6): [-1.90e-15, -3.0e-15, '-'],
              (-6.6, -6.0): [-3.0e-15, -3.57e-15, '-'], (-6.0, -4.8): [-3.57e-15, -7.90e-15, '-']}

for ((xstart, xend), [ystart, yend, ls]) in dict_halpha.items():
    ax_halpha.plot([xstart, xend], [ystart, yend], linewidth=4, alpha=0.4, linestyle=ls, color='r')
for ((xstart, xend), [ystart, yend, ls]) in dict_hbeta.items():
    ax_hbeta.plot([xstart, xend], [ystart, yend], linewidth=4, alpha=0.4, linestyle=ls, color='b')

ax_halpha.text(x=-15, y=-0.6e-15, s=r'$\rm H_{\alpha}\ 6563\ \AA$', fontsize=18)
ax_hbeta.text(x=-15, y=-0.6e-15, s=r'$\rm H_{\beta}\ 4861\ \AA$', fontsize=18)
ax_halpha.axvline(0, linestyle='--', color='k', linewidth=0.8)
ax_hbeta.axvline(0, linestyle='--', color='k', linewidth=0.8)

ax_halpha.set_xlim(-16, 7)
ax_halpha.set_ylim(-8.1e-15, 0.8e-15)
ax_halpha.set_yticklabels([])
ax_halpha.xaxis.set_ticks_position('both')
ax_halpha.yaxis.set_ticks_position('both')
ax_halpha.xaxis.set_major_locator(MultipleLocator(5))
ax_halpha.xaxis.set_minor_locator(MultipleLocator(1))
ax_halpha.yaxis.set_major_locator(MultipleLocator(1e-15))
ax_halpha.yaxis.set_minor_locator(MultipleLocator(1e-16))
ax_halpha.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=15)
ax_halpha.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=15)

ax_hbeta.set_xlim(-16, 7)
ax_hbeta.set_yticklabels([])
ax_hbeta.xaxis.set_ticks_position('both')
ax_hbeta.yaxis.set_ticks_position('both')
ax_hbeta.xaxis.set_major_locator(MultipleLocator(5))
ax_hbeta.xaxis.set_minor_locator(MultipleLocator(1))
ax_hbeta.yaxis.set_major_locator(MultipleLocator(1e-15))
ax_hbeta.yaxis.set_minor_locator(MultipleLocator(1e-16))
ax_hbeta.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=15)
ax_hbeta.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=15)

ax_halpha.set_xlabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=18)
ax_hbeta.set_xlabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=18)
ax_hbeta.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=18)

fig.subplots_adjust(wspace=0.02)
fig.savefig('PLOT_BalmerEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Spectra In Comparison With Models From Dessart
# ------------------------------------------------------------------------------------------------------------------- #

def plot_desspec(ax_obj, file_name, phase='0d', offset=0, smooth=True, sp=3, z=1):
    """
    Plots the spectrum of the object along with model spectra from Dessart et al.(2013).
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        phase       : Phase of the spectrum to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        z           : String to specify the metallicity of the model spectrum
    Returns:
        None
    """
    data_df = pd.read_csv(file_name, sep='\s+', names=['Wavelength', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'] / data_df['Flux'].mean()

    z = str(z) + r' $\rm Z_{\odot}$'

    if file_name[0:3] != 'rfz':
        label = file_name.split('/')[-2]
        color = 'k'
        ax_obj.text(5000, offset + 0.6, s=label + ' +' + phase, color=color, fontsize=10)
        ax_obj.text(7000, offset + 0.5, s=z, color=color, fontsize=10)
    else:
        label = name_SN
        color = 'r'
        ax_obj.text(7250, 0.3, s=label + ' +' + phase, color=color, fontsize=10)

        for index, row in data_df.iterrows():
            if 3500 < row[0] < 5100:
                row[1] = row[1] * 1.35
            if 5200 < row[0] < 6100:
                row[1] = row[1] * 1.1
            if 6500 < row[0] < 7500:
                row[1] = row[1] * 0.9
            if 8200 > row[0] >= 7500:
                row[1] = row[1] * 0.8
            elif row[0] >= 8200:
                row[1] = row[1] * 0.6

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    data_df = data_df[(data_df['Wavelength'] > 3650) & (data_df['Wavelength'] < 9000)]
    ax_obj.plot(data_df['Wavelength'], np.log10(data_df['Flux']) + offset, lw=1, c=color, alpha=0.8,
                label=label + ' ' + phase)


fig = plt.figure(figsize=(16, 8))
ax = fig.add_subplot(121)
ax2 = fig.add_subplot(122, sharey=ax)

plot_desspec(ax, 'rfz_2016-10-14_2016gfy.dat', '31.8 d', offset=0.2)
plot_desspec(ax, DIR_SNe + 'ModDes_Data/m15z2m3/m15_du_sch_z2m3_FeC_mix0p4_13.fl', '31.5 d', -0.9, z='0.1')
plot_desspec(ax, DIR_SNe + 'ModDes_Data/m15z8m3/m15_du_sch_z8m3_FeC_mix0p4_13.fl', '30.9 d', -1.3, z='0.4')
plot_desspec(ax, DIR_SNe + 'ModDes_Data/m15z2m2/m15_du_sch_FeC_mix0p4_13.fl', '30.7 d', -0.5, z='1')
plot_desspec(ax, DIR_SNe + 'ModDes_Data/m15z4m2/m15_du_sch_z4m2_FeC_mix0p4_13.fl', '31.2 d', -1.7, z='2')

plot_desspec(ax2, 'rfz_2016-11-16_2016gfy.dat', '64.8 d', offset=0.25)
plot_desspec(ax2, DIR_SNe + 'ModDes_Data/m15z2m3/m15_du_sch_z2m3_FeC_mix0p4_20.fl', '61.3 d', -0.9, z='0.1')
plot_desspec(ax2, DIR_SNe + 'ModDes_Data/m15z8m3/m15_du_sch_z8m3_FeC_mix0p4_20.fl', '60.2 d', -1.3, z='0.4')
plot_desspec(ax2, DIR_SNe + 'ModDes_Data/m15z2m2/m15_du_sch_FeC_mix0p4_20.fl', '59.7 d', -0.5, z='1')
plot_desspec(ax2, DIR_SNe + 'ModDes_Data/m15z4m2/m15_du_sch_z4m2_FeC_mix0p4_20.fl', '60.6 d', -1.7, z='2')

ax.set_xlim(3600, 9200)
ax.set_ylim(-1.7, 0.7)
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(250))
ax.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=14)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=16)

ax2.set_xlim(3600, 9200)
ax2.set_yticklabels([])
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_major_locator(MultipleLocator(1000))
ax2.xaxis.set_minor_locator(MultipleLocator(250))
ax2.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=14)
ax2.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=14)
ax2.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)

fig.subplots_adjust(wspace=0.01)
fig.savefig('PLOT_CompDesModSpec.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Spectra In Comparison With Models From Jerkstrand
# ------------------------------------------------------------------------------------------------------------------- #

def plot_jerkspec(ax_obj, file_name, phase='0d', offset=0, color='k', smooth=True, sp=3, mass=r'$\rm 12\ M_{\odot}$'):
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
    data_df = data_df[(data_df['Wavelength'] > 3650) & (data_df['Wavelength'] < 9100)]
#     data_df['Flux'] = data_df['Flux'] / data_df['Flux'].mean()

    if file_name[0:3] != 'rfz':
        label = mass
        linestyle = '--'
        ax_obj.text(7450, offset + 0.8, s=label + ' (' + phase + ')', color=color, fontsize=10)
        data_df['Flux'] = data_df['Flux'].apply(lambda x: x * (0.044 / 0.062) * ((5.5 / dist_val) ** 2))
    else:
        label = name_SN
        linestyle = '-'
        ax_obj.text(7450, 1.2, s=label + ' (' + phase + ')', color=color, fontsize=10)

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    ax_obj.plot(data_df['Wavelength'], 17 + np.log10(data_df['Flux']) + offset, lw=1.2, ls=linestyle,
                c=color, alpha=0.9, label=label + ' +' + phase)


fig = plt.figure(figsize=(11, 9))
ax = fig.add_subplot(111)

plot_jerkspec(ax, 'rfz_2017-04-13_2016gfy.dat', '216 d')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams12_212d.dat', '212 d', -0.8, color='orange')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams15_212d.dat', '212 d', -1.6, color='r', mass=r'$\rm 15\ M_{\odot}$')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams19_212d.dat', '212 d', -2.4, color='b', mass=r'$\rm 19\ M_{\odot}$')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams25_212d.dat', '212 d', -3.2, color='g', mass=r'$\rm 25\ M_{\odot}$')

ax.set_xlim(3500, 9300)
ax.set_ylim(-4.2, 2.1)
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.2))
ax.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=16)
ax.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=16)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=16)

fig.subplots_adjust(wspace=0.01)
fig.savefig('PLOT_CompJerkModSpec.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Spectra In Comparison With Models From Jerkstrand
# ------------------------------------------------------------------------------------------------------------------- #

def plot_jerkspec2(ax_obj, file_name, phase, offset=0, color='k', smooth=False, sp=2, mass=r'$\rm 12\ M_{\odot}$'):
    """
    Plots the spectrum of the object along with model spectra from Jerkstrand et al.(2014).
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
    data_df = data_df[(data_df['Wavelength'] > 3650) & (data_df['Wavelength'] < 9100)]

    if file_name[0:3] != 'rfz':
        label = mass
        data_df['Flux'] = data_df['Flux'].apply(lambda x: x * (0.044 / 0.062) * ((5.5 / dist_val) ** 2))
    else:
        label = name_SN
        data_df['Flux'] = data_df['Flux'] - 3e-17
#         for index, row in data_df.iterrows():
#             if 6700 < row[0] < 7800:
#                 row[1] = row[1] - 3e-17
#             if 7900 < row[0]:
#                 row[1] = row[1] + 3e-17

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    ax_obj.plot(data_df['Wavelength'], data_df['Flux'] + offset, lw=1.2, c=color, alpha=0.9,
                label=label + ' (' + phase + ')')


fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)

plot_jerkspec(ax, 'rfz_2017-04-13_2016gfy.dat', '216 d')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams12_212d.dat', '212 d', color='orange', mass=r'$\rm 12\ M_{\odot}$')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams15_212d.dat', '212 d', color='r', mass=r'$\rm 15\ M_{\odot}$')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams19_212d.dat', '212 d', color='b', mass=r'$\rm 19\ M_{\odot}$')
plot_jerkspec(ax, DIR_SNe + 'ModJerk_Data/mzams25_212d.dat', '212 d', color='g', mass=r'$\rm 25\ M_{\odot}$')

ax.set_xlim(3600, 9200)
ax.legend(markerscale=2, frameon=False, fontsize=14)
# ax.set_ylim(-4.0, 1.1)
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=18)
ax.tick_params(which='minor', direction='in', length=4, width=0.9, labelsize=18)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=18)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$]', fontsize=18)

fig.subplots_adjust(wspace=0.01)
fig.savefig('PLOT_CompJerkModSpec2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
#  Plots The Spectra In Comparison With Spectra Of Other Well-Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def plot_compspec(ax_obj, path_name, phase='0d', offset=0, smooth=False, sp=1, clip_str='3600:9050'):
    """
    Plots the spectrum for comparison with other well-studied SNe.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        file_name   : Name of the 1-D Spectrum FITS file to be plotted
        phase       : Phase of the spectrum to be plotted
        offset      : Offset in Flux units to be applied to the spectra
        smooth      : Should the spectrum be smoothened before plotting?
        sp          : Smoothing parameter to be used for Gaussian Kernel
        clip_str    : Wavelength region to be used for plotting
    Returns:
        None
    """
    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'] / np.median(data_df['Flux'])

    lower_clip, upper_clip = clip_str.split(':')
    file_name = path_name.split('/')[-1]

    if file_name[0:3] != 'rfz':
        name = file_name.split('_')[0][2:]
        phase = file_name.split('_')[-1]
    else:
        name = name_SN[2:]
        copy_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]
        ax_obj.plot(copy_df['Wave'], copy_df['Flux'], lw=0.8, alpha=0.5, label='_nolegend')

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))
    if name in dict_redshift.keys():
        data_df['Wave'] = data_df['Wave'].apply(lambda x: (1 - dict_redshift[name]) * x)

    data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]
    data_df['Flux'] = data_df['Flux'] - offset

    ax_obj.plot(data_df['Wave'], data_df['Flux'], lw=1, label=phase + ' ' + name)
    ax_obj.text(data_df['Wave'].tolist()[-1] + 50, data_df['Flux'].tolist()[-1], s=name + ' [+' + phase + ']',
                fontsize=10)

    ax_obj.set_yticklabels([])
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(1000))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(100))
    ax_obj.yaxis.set_major_locator(MultipleLocator(2))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.4))
    ax_obj.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
    ax_obj.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)


sns.set_palette(sns.cubehelix_palette(7, start=1.5, rot=2, dark=0.2, light=.6, reverse=True))

dict_redshift = {'2007Y': 0.00463, '2007uy': 0.00644, '2008D': 0.00644, '2009jf': 0.007942,
                 '2012au': 0.004546, '13bvn': 0.005}

list_specmax = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/NearExp/*')
list_specearlyplat = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/EarlyPlat/*')
list_speclateplat = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/LatePlat/*')
list_specnebular = group_similar_files('', common_text=DIR_SNe + 'Spec_Data/Nebular/*')

fig1 = plt.figure(figsize=(9, 10))

ax1 = fig1.add_subplot(211)
plot_compspec(ax1, DIR_SPEC + 'rfz_2016-09-13_2016gfy.dat', phase='04d', smooth=True, sp=3)
for index, path_name in enumerate(list_specmax):
    plot_compspec(ax1, path_name, offset=index * 1.1 + 0.8, smooth=False)

ax1.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax1.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax1.text(8000, 3.1, s='Near Maximum', color='r', fontsize=16)

ax2 = fig1.add_subplot(212, sharex=ax1)
plot_compspec(ax2, DIR_SPEC + 'rfz_2016-10-04_2016gfy.dat', phase='25d', smooth=True, sp=3)
for index, path_name in enumerate(list_specearlyplat):
    plot_compspec(ax2, path_name, offset=index * 1.2 + 1.0, smooth=False)

ax2.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax2.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax2.text(7700, 2.5, s='Early Plateau Phase', color='r', fontsize=16)

ax2.set_xlim(3400, 9950)
ax2.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
fig1.text(0.09, 0.5, r'Normalized Flux + Const.', va='center',
          rotation='vertical', fontsize=15)
fig1.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig1.savefig('PLOT_SpecComp1.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig1)

fig2, (ax3, ax4) = plt.subplots(2, 1, figsize=(9, 11), gridspec_kw={'height_ratios': [4, 5]}, sharex=True)

plot_compspec(ax3, DIR_SPEC + 'rfz_2016-11-24_2016gfy.dat', phase='76d', smooth=True, sp=2.2, clip_str='3700:9100')
for index, path_name in enumerate(list_speclateplat):
    plot_compspec(ax3, path_name, offset=index * 1.4 + 1.4, smooth=True, sp=2, clip_str='3700:9100')

plot_compspec(ax4, DIR_SPEC + 'rfz_2017-04-13_2016gfy.dat', phase='216d', smooth=True, sp=3, clip_str='3700:9100')
for index, path_name in enumerate(list_specnebular):
    plot_compspec(ax4, path_name, offset=index * 2 + 2.5, smooth=True, sp=2, clip_str='3700:9100')

ax3.text(6777, 1.4, s=r'$\rm \oplus$', fontsize=10)
ax3.text(7520, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax3.text(8000, 3.3, s='Late Plateau Phase', color='r', fontsize=16)

ax4.text(6777, 2.6, s=r'$\rm \oplus$', fontsize=10)
ax4.text(7520, 2.4, s=r'$\rm \oplus$', fontsize=10)
ax4.text(8500, 10.0, s='Nebular Phase', color='r', fontsize=16)

ax4.set_ylim(-7.5, 12)
ax4.set_xlim(3500, 10250)
ax4.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
fig2.text(0.085, 0.5, r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', va='center',
          rotation='vertical', fontsize=15)

fig2.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig2.savefig('PLOT_SpecComp2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Halpha Early Evolution
# ------------------------------------------------------------------------------------------------------------------- #

def plot_linecomp(ax_obj, path_name, phase='0d', offset=0, smooth=False, sp=1, legend=True, clip_str='6000:7000'):
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
        clip_str    : Wavelength region to be used for plotting
    Returns:
        None
    """
    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'] / np.median(data_df['Flux'])

    lower_clip, upper_clip = clip_str.split(':')
    file_name = path_name.split('/')[-1]

    if file_name[0:3] != 'rfz':
        name = file_name.split('_')[0]
        phase = file_name.split('_')[-1]
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
        ax_obj.plot(data_df['Vel'], np.log10(data_df['Flux']), lw=1.5, label=' +' + phase + '')
    else:
        ax_obj.plot(data_df['Vel'], np.log10(data_df['Flux']), lw=1.5, label=name + ' [+' + phase + ']')

    if not legend:
        ax_obj.text(-27, np.log10(data_df['Flux']).tolist()[0] + 0.03, fontsize=12, s=name + ' [+' + phase + ']')


list_compHalpha = ['BoxyHalpha/2007od_5.5d', 'BoxyHalpha/2007od_9.2d', 'BoxyHalpha/2016esw_19.5d']
list_compCaII = ['BoxyHalpha/2007od_5.5d' + 'BoxyHalpha/2007od_5.5d', 'BoxyHalpha/2007od_5.5d']

sns.set_palette(sns.color_palette('colorblind', 10))
fig_ha, (ax_ha, ax_comp) = plt.subplots(2, 1, figsize=(8, 14), gridspec_kw={'height_ratios': [2, 3]}, sharex=True)

plot_linecomp(ax_ha, DIR_SPEC + 'rfz_2016-09-20_2016gfy.dat', phase='11d', smooth=True, sp=1)
plot_linecomp(ax_ha, DIR_SPEC + 'rfz_2016-09-27_2016gfy.dat', phase='18d', smooth=True, sp=1)
plot_linecomp(ax_ha, DIR_SPEC + 'rfz_2016-09-30_2016gfy.dat', phase='21d', smooth=True, sp=1)

ax_ha.set_yticklabels([])
ax_ha.set_ylim(-0.18, 0.44)
ax_ha.text(6.810, 0.05, s=r'$\rm \oplus$', fontsize=10)
ax_ha.text(-25, 0.38, s=r'SN 2016gfy', fontsize=14)
ax_ha.axvline(0, lw=1, ls='--', color='k')

ax_ha.legend(markerscale=2, frameon=False, fontsize=14)
ax_ha.yaxis.set_ticks_position('both')
ax_ha.xaxis.set_ticks_position('both')
ax_ha.xaxis.set_major_locator(MultipleLocator(10))
ax_ha.xaxis.set_minor_locator(MultipleLocator(1))
ax_ha.yaxis.set_major_locator(MultipleLocator(0.1))
ax_ha.yaxis.set_minor_locator(MultipleLocator(0.01))
ax_ha.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
ax_ha.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)
ax_ha.set_ylabel(r'Log $\rm\ F_{\lambda}$', fontsize=15)

plot_linecomp(ax_comp, DIR_SPEC + 'rfz_2016-09-20_2016gfy.dat', phase='11d', offset=0.1, legend=False, smooth=False)
for index, path_name in enumerate(list_compHalpha):
    plot_linecomp(ax_comp, DIR_SNe + 'Spec_Data/' + path_name, offset=index * 0.2 + 0.2, legend=False, smooth=False)

ax_comp.set_ylim(-0.66, 0.4)
ax_comp.set_yticklabels([])
ax_comp.axvline(0, lw=1, ls='--', color='k')

ax_comp.yaxis.set_ticks_position('both')
ax_comp.xaxis.set_ticks_position('both')
ax_comp.xaxis.set_major_locator(MultipleLocator(10))
ax_comp.xaxis.set_minor_locator(MultipleLocator(1))
ax_comp.yaxis.set_major_locator(MultipleLocator(0.2))
ax_comp.yaxis.set_minor_locator(MultipleLocator(0.02))
ax_comp.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
ax_comp.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)
ax_comp.set_xlabel(r'Velocity $\rm [\times 10^3\ km s^{-1}$]', fontsize=15)
ax_comp.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=15)

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
    ax_obj.plot(data_df['Wave'], np.log10(data_df['Flux']) + offset, lw=1.2, c='grey', alpha=0.5, label='_nolegend_')

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    lower_clip, upper_clip = clip_str.split(':')
    if int(upper_clip) > data_df['Wave'].max():
        data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] <= data_df['Wave'].max() - 50)]
    else:
        data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]

#     ax_obj.plot(data_df['Wave'], data_df['Flux'], linewidth=1.1, label='_nolegend_')
    ax_obj.plot(data_df['Wave'], np.log10(data_df['Flux']) + offset, lw=1.3, c='k', alpha=0.6, label='_nolegend_')


list_specgal = group_similar_files('', common_text=DIR_GAL + 'z_*' + name_hostgal + '*.dat')

fig_host = plt.figure(figsize=(13, 8))
ax_host = fig_host.add_subplot(111)

for index, path_name in enumerate(list_specgal):
    plot_galspec(ax_host, path_name, offset=index * -0.8, smooth=True, sp=2.5)

# ax_host.axvline(4861, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(4959, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(5007, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6548, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6563, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6584, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6717, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')
# ax_host.axvline(6731, color='k', ls='-.', lw=1, alpha=0.4, label='_nolegend_')

ax_host.text(4740, -15.1, color='k', s=r'$\rm H\beta$', fontsize=13)
ax_host.text(5040, -15.0, color='k', s=r'$\rm [O\,III]$', fontsize=13)
ax_host.text(6200, -14.9, color='k', s=r'$\rm H\alpha\,+\,[N\,II]$', fontsize=13)
ax_host.text(6750, -15.25, color='k', s=r'$\rm [S\,II]$', fontsize=13)
ax_host.text(3950, -15.05, color='r', s=r'HII Region', fontsize=15)
ax_host.text(3950, -15.75, color='b', s=r'Nucleus', fontsize=15)

ax_host.set_ylim(-16.1, -14.45)
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

ax_host.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)
ax_host.set_ylabel(r'Log $\rm F_{\lambda}$ + Const.', fontsize=16)

fig_host.savefig('PLOT_SpecHostGal.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_host)
# ------------------------------------------------------------------------------------------------------------------- #
