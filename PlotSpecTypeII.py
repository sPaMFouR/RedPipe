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
sns.set_palette(sns.cubehelix_palette(8, start=3, rot=0.5, light=0.6))
# sns.set_palette(sns.color_palette('rocket_r', 10))
# sns.set_palette(sns.dark_palette(color='muted blue', input='xkcd', n_colors=8))
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
redshift = 0.008
phase_plateauend = 110
phase_nebstart = 110
date_explosion = 2457644.60
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/IIP_Data/Spec_Data/"
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

def plot_velevol(ax_obj, file_df, index, wavelength=6563, offset=4.1e-16, smooth=False, sp=2, plot_label=True):
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
        file_df.loc[:, 'Flux'] += 1.1e-15
    if smooth:
        file_df.loc[:, 'Flux'] = convolve(file_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    file_df['Velocity'] = file_df['Wavelength'].apply(lambda x: ((x - wavelength) / wavelength) * light_speed / 1e3)
    file_df = file_df[(file_df['Velocity'] > -18) & (file_df['Velocity'] < 6)]
    file_df.loc[:, 'Flux'] -= offset * index

    if plot_label:
        ax_obj.text(x=file_df['Velocity'].values[-1] + 1.2, y=file_df['Flux'].values[-1],
                    s=evolution_df.loc[index, 'Label'], fontsize=12)

    ax_obj.plot(file_df['Velocity'], file_df['Flux'])
    ax_obj.axvline(0, linestyle='--', color='k', linewidth=0.8)

    ax_obj.set_xlim(-16, 7)
    ax_obj.set_ylim(-7.8e-15, 0.1e-15)
    ax_obj.set_yticklabels([])
#     ax_obj.grid(which='both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(5))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(1))
    ax_obj.yaxis.set_major_locator(MultipleLocator(1e-15))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(1e-16))
    ax_obj.tick_params(which='major', direction='in', length=8, width=1.6, labelsize=13)
    ax_obj.tick_params(which='minor', direction='in', length=4, width=0.8, labelsize=13)


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

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files',
#                         choices=['Yes', 'No'])
# ctext = eg.enterbox('Common Text Of Files To Be Flux Calibrated?', title='Flux Calibration', default='fz_*.fits')
# bool_smooth = eg.boolbox('Perform Smoothening Of Spectra?', title='Smoothening 1-D Spectra', choices=['Yes', 'No'])
# clip_str = eg.enterbox('Specify Clipping Section: ', title='Clipping Of 1-D Spectra', default='4000:8500')
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
early_df = spec_df[spec_df['Phase'] < 30].copy()
late_df = spec_df[(spec_df['Phase'] > 30) & (spec_df['Phase'] < phase_plateauend)].copy()
plateau_df = spec_df[spec_df['Phase'] < phase_plateauend].copy()
nebular_df = spec_df[spec_df['Phase'] > phase_nebstart].copy()

evolution_df = spec_df.copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
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
ax.axvline(4861, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(4959, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(5007, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6563, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6717, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6731, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')

for (line, [wavelength, width]) in dict_labelsplat.items():
    ax.text(wavelength - 30, 5.6e-15, line, rotation='vertical', fontsize=10)
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(4e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.4e-15))
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=12)
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Early Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 7))
ax = fig.add_subplot(111)

for index, file_name in early_df['TextFile'].items():
    plot_epoch(ax, file_name, index, early_df, offset=0.7e-15)

dict_labelsearly = {r'$\rm Ca\,II$ (H & K)': [3800, 60], r'$\rm H\delta$ 4102': [3960, 50],
                    r'$\rm H\gamma$ 4340': [4190, 50], r'$\rm H\beta$ 4861': [4700, 70],
                    r'$\rm Fe\,I$ 5169': [5070, 50],
                    r'$\rm He\,I$ 5876': [5660, 70], r'$\rm H\alpha$ 6563': [6360, 60],
                    r'$\rm \oplus[O_2$ 6867 B-band]': [6820, 40], r'$\rm \oplus[H_{2}O$ 7165]': [7115, 40],
                    r'$\rm \oplus[O_2$ 7620 A-band]': [7560, 50], r'$\rm \oplus[H_{2}O$ band]': [8100, 40], }

ax.set_ylim(-4.5e-15, 4e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)
ax.axvline(4861, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(4959, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(5007, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6563, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6717, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6731, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')

for (line, [wavelength, width]) in dict_labelsearly.items():
    ax.text(wavelength - 30, 3.6e-15, line, fontsize=9, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(2e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=12)
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecEarlyPlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Late Plateau Phase
# ------------------------------------------------------------------------------------------------------------------- #
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

ax.set_ylim(-17.4e-15, -1.5e-15)
ax.set_xlim(int(lower_lim) - 100, int(upper_lim) + 550)
ax.axvline(4861, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(4959, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(5007, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6563, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6717, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6731, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')

for (line, [wavelength, width]) in dict_labelslate.items():
    ax.text(wavelength - 30, -1.9e-15, line, fontsize=9, rotation='vertical')
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(2e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=12)
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecLatePlateau.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Nebular Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 10))
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
ax.axvline(4861, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(4959, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(5007, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6563, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6717, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')
ax.axvline(6731, color='k', linestyle='--', linewidth=0.8, alpha=0.4, label='_nolegend_')

for (line, [wavelength, width]) in dict_labelsneb.items():
    ax.text(wavelength - 30, -3.4e-15, line, rotation='vertical', fontsize=9)
    ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor='silver', alpha=0.7)

ax.set_yticklabels([])
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_major_locator(MultipleLocator(0.5e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.05e-15))
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=12)
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecNebular.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Balmer Features Evolution Feature Across Various Epochs
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('rocket_r', 6))
evolution_df = evolution_df[(evolution_df['Phase'] > 30) & (evolution_df['Phase'] < 110)]

fig, (ax_hbeta, ax_halpha) = plt.subplots(1, 2, figsize=(9, 16), sharey=True)

for index, file_name in evolution_df['TextFile'].items():
    data_df = pd.read_csv(file_name, names=['Wavelength', 'Flux'], sep='\s+', dtype='float64')
    data_df = data_df[data_df['Wavelength'] <= 6800]

    plot_velevol(ax_hbeta, data_df, index, wavelength=4861, plot_label=False, smooth=True)
    plot_velevol(ax_halpha, data_df, index, offset=4.7e-16, smooth=True)

dict_halpha = {(-9.5, -9.0): [-2.5e-15, -4.4e-15, '-.'], (-9.0, -8.9): [-4.4e-15, -4.96e-15, '-.'],
               (-8.9, -8.7): [-4.96e-15, -7.53e-15, '-.'], (-8.5, -7.8): [-1.65e-15, -3.15e-15, '-'],
               (-7.8, -7.0): [-3.15e-15, -3.65e-15, '-'], (-7.0, -6): [-3.65e-15, -7.64e-15, '-']}

dict_hbeta = {(-9.5, -9.0): [-2.7e-15, -7.47e-15, '-.'], (-6.9, -6.6): [-2.23e-15, -3.0e-15, '-'],
              (-6.6, -6.0): [-3.0e-15, -3.57e-15, '-'], (-6.0, -4.8): [-3.57e-15, -7.63e-15, '-']}

for ((xstart, xend), [ystart, yend, ls]) in dict_halpha.items():
    ax_halpha.plot([xstart, xend], [ystart, yend], linewidth=1, linestyle=ls, color='k')
for ((xstart, xend), [ystart, yend, ls]) in dict_hbeta.items():
    ax_hbeta.plot([xstart, xend], [ystart, yend], linewidth=1, linestyle=ls, color='k')

ax_halpha.text(x=-14, y=-0.3e-15, s=r'$\rm H_{\alpha}\ 6563\ \AA$', fontsize=18)
ax_hbeta.text(x=-14, y=-0.3e-15, s=r'$\rm H_{\beta}\ 4861\ \AA$', fontsize=18)

ax_halpha.set_xlabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=14)
ax_hbeta.set_xlabel(r'Velocity [$\rm x10^3\ km\ s^{-1}$]', fontsize=14)
ax_hbeta.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.subplots_adjust(wspace=0.03)
fig.savefig('PLOT_BalmerEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
