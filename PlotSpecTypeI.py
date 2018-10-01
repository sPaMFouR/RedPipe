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
import seaborn as sns
from pyraf import iraf
from jdcal import gcal2jd
from astropy.io import fits
import matplotlib.pyplot as plt
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel

sns.set_style('ticks')
sns.set_palette(sns.dark_palette((260, 75, 60), input='husl', n_colors=6))
# sns.set_palette(sns.color_palette('RdBu', 10))
# sns.set_palette(sns.color_palette('rocket_r', 10))
# sns.set_palette(sns.hls_palette(10, l=.4, s=0.9))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
JD_keyword = 'JD'
lower_lim = 3650
upper_lim = 9050
light_speed = 2.99792458e5
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Details Of SN In Study (2016gfy)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
phase_nebstart = 60
redshift = 0.006191
date_explosion = 2458096.2
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/Ib_Data/Spec_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2017iro/Spectroscopy/"
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

        output_spectrum = spectrum.rstrip('.fits') + out_ext
        spec_df.to_csv(output_spectrum, sep=' ', index=True, header=None)
        list_outspectra.append(output_spectrum)

    return list_outspectra

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Plot Individual Spectra
# ------------------------------------------------------------------------------------------------------------------- #

def plot_epoch(ax_obj, file_name, index, master_df, offset=0.6e-15, smooth=False, sp=2):
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
    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    if data_df['Wavelength'].max() >= int(upper_lim):
        data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    else:
        data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) &
                          (data_df['Wavelength'] <= data_df['Wavelength'].max() - 50)]

    data_df = data_df[(data_df['Wavelength'] >= int(lower_lim)) & (data_df['Wavelength'] <= int(upper_lim))]
    data_df['Flux'] -= offset * index

    if master_df.loc[index, 'Phase'] > 10:
        data_df = data_df[(data_df['Wavelength'] > 3800) & (data_df['Wavelength'] < 9020)]

    ax_obj.plot(data_df['Wavelength'], data_df['Flux'], linewidth=1, label='_nolegend_')
    ax_obj.text(x=data_df['Wavelength'].values[-1] + 50, y=data_df['Flux'].values[-1],
                s=master_df.loc[index, 'Label'], fontsize=10)


def plot_spec(ax_obj, path_name, phase='0d', offset=0, label_off=0.2, smooth=False, sp=1, clip_str='3600:9100'):
    file_name = path_name.split('/')[-1]

    if file_name[0:3] != 'rfz':
        name = file_name.split('_')[0]
        phase = file_name.split('_')[-1].split('.')[0]
    else:
        name = name_SN

    if name[0:2] == '20':
        name = name[2:]

    data_df = pd.read_csv(path_name, sep='\s+', names=['Wave', 'Flux'], header=None, comment='#')
    data_df['Flux'] = data_df['Flux'] / np.median(data_df['Flux'])

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))
    if name in dict_redshift.keys():
        data_df['Wave'] = data_df['Wave'].apply(lambda x: (1 - dict_redshift[name]) * x)

    if name == '09jf' and phase == '-7d':
        offset += -1.3
    elif name == '07uy' and phase == '-5d':
        offset += 0.5

    lower_clip, upper_clip = clip_str.split(':')
    data_df = data_df[(data_df['Wave'] > int(lower_clip)) & (data_df['Wave'] < (int(upper_clip)))]
    data_df['Flux'] = data_df['Flux'] - offset

    ax_obj.plot(data_df['Wave'], data_df['Flux'], linewidth=1, label=phase + ' ' + name)
    ax_obj.text(data_df['Wave'].tolist()[-1] + 10, data_df['Flux'].tolist()[-1] + label_off,
                s=name + ' (' + phase + ')', fontsize=10)

    ax_obj.set_yticklabels([])
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(1000))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(100))
    ax_obj.yaxis.set_major_locator(MultipleLocator(2))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax_obj.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
    ax_obj.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)

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
specpremax_df = spec_df[spec_df['Phase'] < 0].copy()
specpostmax_df = spec_df[(spec_df['Phase'] >= 0) & (spec_df['Phase'] < 10)].copy()

specmax_df = spec_df[spec_df['Phase'] < 10].copy()
specearly_df = spec_df[(spec_df['Phase'] > 10) & (spec_df['Phase'] < phase_nebstart)].copy()
specnebular_df = spec_df[(spec_df['Phase'] > phase_nebstart)].copy()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra Around The Maximum
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(13, 9))
ax = fig.add_subplot(111)

for index, file_name in specmax_df['TextFile'].items():
    plot_epoch(ax, file_name, index, specmax_df, offset=1.2e-15, smooth=True)

dict_labels = {r'$\rm Ca\,II$ (H & K)': [3760, 90, 'teal'], r'$\rm He\,I,\ Fe\,II,\ Mg\,II$': [4400, 80, 'teal'],
               r'$\rm Fe\,II\ 4924, 5018$': [4805, 0, 'teal'], '+        ': [4890, 170, 'teal'],
               r'$\rm Fe\,II\ 5169$  ': [4980, 0, 'teal'], r'$\rm He\,I$ 5876': [5690, 60, 'r'],
               r'$\rm Si\,II$ 6355': [6220, 80, 'teal'], r'$\rm He\,I$ 6678': [6500, 50, 'r'],
               r'$\rm He\,I$ 7065': [6910, 60, 'r'], r'$\rm Ca\,II\ NIR\ Triplet$': [8310, 150, 'teal']}

ax.set_xlim(lower_lim - 150, upper_lim + 450)
ax.set_ylim(-11.3e-15, 5.4e-15)
ax.set_yticklabels([])

ax.axvline(6563, color='k', linestyle='--', linewidth=1.2, alpha=0.4, label='_nolegend_')
ax.axvline(6583, color='k', linestyle='--', linewidth=1.2, alpha=0.4, label='_nolegend_')
ax.axvline(4861, color='k', linestyle='--', linewidth=1.2, alpha=0.4, label='_nolegend_')
ax.axvline(6721, color='k', linestyle='--', linewidth=1.2, alpha=0.4, label='_nolegend_')

ax.annotate(r'$\rm Na\,I$D Host', xy=(5896, 1.6e-15), xytext=(5856, 4.8e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')
ax.annotate(r'$\rm \oplus$', xy=(6827, 0.9e-15), xytext=(6787, 2.0e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')
ax.annotate(r'$\rm \oplus$', xy=(7570, 0.6e-15), xytext=(7530, 1.8e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')
ax.annotate(r'$\rm Fe\,II,\ Sc\,II$', xy=(5430, -8.9e-15), xytext=(5390, -9.6e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')

for (line, [wavelength, width, color]) in dict_labels.items():
    ax.text(wavelength - 30, 4.8e-15, line, rotation='vertical', fontsize=10)
    if width != 0:
        ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor=color, alpha=0.2)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(5e-15))
ax.yaxis.set_minor_locator(MultipleLocator(1e-15))
ax.tick_params(which='major', direction='in', length=8, width=1, labelsize=13)
ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=13)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=16)

fig.savefig('PLOT_SpecAroundMax.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Post Maximum Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111)

for index, file_name in specearly_df['TextFile'].items():
    plot_epoch(ax, file_name, index, specearly_df, offset=0.4e-15, smooth=True, sp=5)

ax.set_xlim(lower_lim, upper_lim + 450)
ax.set_ylim(-7e-15, -1.1e-15)
ax.set_yticklabels([])

dict_labels = {r'$\rm He\,I\ 4471\ +\ Mg\,II\ 4481$': [4365, 50, 'teal'],
               r'$\rm Fe\,II\ 4924, 5018$': [4880, 0, 'teal'], '+        ': [4950, 160, 'teal'],
               r'$\rm Fe\,II\ 5169$  ': [5040, 0, 'teal'], r'$\rm Fe\,II\ 5535\ +\ Sc\,II\ 5527$': [5410, 50, 'teal'],
               r'$\rm Sc\,II$ 5663 Multiplet': [5570, 40, 'teal'], r'$\rm He\,I$ 5876': [5705, 55, 'r'],
               r'$\rm He\,I$ 6678': [6500, 50, 'r'], r'$\rm He\,I$ 7065': [6860, 70, 'r'],
               r'$\rm He\,I$ 7281': [7140, 60, 'r'], r'$\rm O\,I\ 7774$': [7650, 60, 'teal'],
               r'$\rm Ca\,II\ NIR\ Triplet$': [8410, 160, 'teal']}

ax.axvline(6563, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(6583, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(4861, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(6721, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')

ax.annotate(r'$\rm Na\,I$D Host', xy=(5896, -2.3e-15), xytext=(5856, -1.4e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')
ax.annotate(r'$\rm \oplus$', xy=(6827, -3.1e-15), xytext=(6787, -2.5e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')
ax.annotate(r'$\rm \oplus$', xy=(7570, -3.3e-15), xytext=(7530, -2.8e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')

for (line, [wavelength, width, color]) in dict_labels.items():
    ax.text(wavelength - 30, -1.3e-15, line, rotation='vertical', fontsize=10)
    if width != 0:
        ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor=color, alpha=0.2)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(1e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.2e-15))
ax.tick_params(which='major', direction='in', length=8, width=1, labelsize=12)
ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecPostMax.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The 1-Dimensional Spectra During The Nebular Phase
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(13, 9))
ax = fig.add_subplot(111)

for index, file_name in specnebular_df['TextFile'].items():
    plot_epoch(ax, file_name, index, specnebular_df, offset=0.2e-15, smooth=True, sp=2)

ax.set_xlim(lower_lim, upper_lim + 500)
ax.set_ylim(-5.65e-15, -2.2e-15)
ax.set_yticklabels([])

dict_labels = {r'$\rm Fe\,II$ Blend': [5000, 90, 'teal'], r'$\rm [O\,I]$ 5577': [5540, 50, 'k'],
               r'$\rm He\,I$ 5876 + Na ID': [5870, 100, 'teal'], r'$\rm [O\,I]$ 6300, 6364': [6290, 60, 'k'],
               r'$\rm He\,I$ 7065': [7010, 60, 'r'], r'$\rm He\,I\ 7281\ +\ [Ca\,II]$ D': [7270, 100, 'teal'],
               r'$\rm O\,I\ 7774$': [7750, 60, 'teal'], r'$\rm Ca\,II$ NIR Triplet': [8470, 0, 'teal'],
               '+      ': [8560, 160, 'teal'], r'$\rm[C\,I]$ 8730  ': [8650, 0, 'teal']}

ax.axvline(6563, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(6583, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(4861, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(5007, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')
ax.axvline(6721, color='k', linestyle='--', linewidth=1.2, alpha=0.3, label='_nolegend_')

ax.annotate(r'$\rm Na\,I$D Host', xy=(5896, -3.1e-15), xytext=(6006, -2.8e-15), fontsize=9,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical', ha='left', va='bottom')
ax.annotate(r'$\rm \oplus$', xy=(6827, -3.3e-15), xytext=(6787, -2.9e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical', va='center')
ax.annotate(r'$\rm \oplus$', xy=(7570, -3.4e-15), xytext=(7530, -3.0e-15), fontsize=10,
            arrowprops=dict(facecolor='blue', arrowstyle='-'), rotation='vertical')

for (line, [wavelength, width, color]) in dict_labels.items():
    ax.text(wavelength - 30, -2.35e-15, line, rotation='vertical', fontsize=9)
    if width != 0:
        ax.fill_betweenx(ax.get_ybound(), wavelength - width - 10, wavelength + width + 20, facecolor=color, alpha=0.2)

ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(500))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(0.5e-15))
ax.yaxis.set_minor_locator(MultipleLocator(0.1e-15))
ax.tick_params(which='major', direction='in', length=8, width=1, labelsize=12)
ax.tick_params(which='minor', direction='in', length=4, width=1, labelsize=12)
ax.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=14)
ax.set_ylabel(r'Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', fontsize=14)

fig.savefig('PLOT_SpecNebular.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
#  Plots The Spectra In Comparison With Spectra Of Other Well-Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #
sns.set_palette(sns.color_palette('cubehelix', 6))

dict_redshift = {'2007Y': 0.00463, '2007uy': 0.00644, '2008D': 0.00644, '2009jf': 0.007942,
                 '2012au': 0.004546, '13bvn': 0.005}

list_specpremax = group_similar_files('', common_text=DIR_SNe + 'PreMax/*')
list_specpostmax = group_similar_files('', common_text=DIR_SNe + 'PostMax/*')
list_specmax = group_similar_files('', common_text=DIR_SNe + 'NearMax/*')
list_specnebular = group_similar_files('', common_text=DIR_SNe + 'Nebular/*')

fig1 = plt.figure(figsize=(9, 9))

ax1 = fig1.add_subplot(211)
plot_spec(ax1, DIR_SPEC + 'rfz_sn17iro_dec01.dat', phase='-7d', smooth=True, sp=3)
for index, path_name in enumerate(list_specpremax):
    plot_spec(ax1, path_name, offset=index * 1.2 + 0.8, smooth=True, sp=2)

ax1.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax1.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax1.text(9300, 2.1, s='A', color='r', fontsize=16)

ax2 = fig1.add_subplot(212, sharex=ax1)
plot_spec(ax2, DIR_SPEC + 'rfz_sn17iro_dec08.dat', phase='0.3d')
for index, path_name in enumerate(list_specmax):
    plot_spec(ax2, path_name, offset=index + 0.9, smooth=True)

ax2.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax2.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax2.text(9300, 2.0, s='B', color='r', fontsize=16)

ax1.set_ylim(-2.2, 2.8)
ax2.set_ylim(-2.7, 2.7)
ax2.set_xlim(3400, 9900)
ax2.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
fig1.text(0.09, 0.5, r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', va='center',
          rotation='vertical', fontsize=15)
fig1.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig1.savefig('PLOT_SpecComp1.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig1)

fig2 = plt.figure(figsize=(9, 9))

ax3 = fig2.add_subplot(211)
plot_spec(ax3, DIR_SPEC + 'rfz_sn17iro_jan08_18.dat', phase='31d', smooth=True, sp=2.2, clip_str='3700:9100')
for index, path_name in enumerate(list_specpostmax):
    plot_spec(ax3, path_name, offset=index * 0.9 + 1.1, smooth=True, sp=2, clip_str='3700:9100')

ax3.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax3.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax3.text(9300, 2.2, s='C', color='r', fontsize=16)

ax4 = fig2.add_subplot(212, sharex=ax3)
plot_spec(ax4, DIR_SPEC + 'rfz_sn17iro_mar15_2018.dat', phase='97d', smooth=True, sp=4, label_off=0.4,
          clip_str='3700:9100')
for index, path_name in enumerate(list_specnebular):
    plot_spec(ax4, path_name, offset=index * 1.1 + 1.2, smooth=True, sp=3, clip_str='3700:9100')

ax4.text(6777, 1.2, s=r'$\rm \oplus$', fontsize=10)
ax4.text(7520, 1.0, s=r'$\rm \oplus$', fontsize=10)
ax4.text(9300, 3.1, s='D', color='r', fontsize=16)

ax3.set_ylim(-2.9, 2.9)
ax4.set_ylim(-2.3, 3.9)
ax4.set_xlim(3500, 9900)
ax4.set_xlabel(r'Rest Wavelength [$\rm \AA$]', fontsize=15)
fig2.text(0.085, 0.5, r'Scaled Flux [$\rm erg\ cm^{-2}\ s^{-1}\ \AA^{-1}$] + Const.', va='center',
          rotation='vertical', fontsize=15)

fig2.subplots_adjust(hspace=0.01, top=0.9, right=0.95)
fig2.savefig('PLOT_SpecComp2.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #
