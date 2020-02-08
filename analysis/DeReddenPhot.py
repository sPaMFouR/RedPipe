#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxx-----------------Fit Planckian Function To A Spectrum---------------xxxxxxxxxxxxxxxxxxxxxxxxx #
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
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
from dust_extinction.dust_extinction import F99, CCM89
from astropy.modeling.blackbody import blackbody_lambda
from astropy.convolution import convolve, Gaussian1DKernel
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
guess_amp = 1e-22
guess_temp = 10000
lower_lim = 3650
upper_lim = 9100

EBV_mag = 0.21
EBV_err = 0.11
date_explosion = 2457641.4

JD_keyword = 'JD'
name_SN = '2016gfy'
file_appmag = '2016gfy_HCT.dat'

DIR_PHOT = '/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/'
DIR_CODE = '/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Data Containing Information On FILTERS and Other Type II SNe
# Extinction Coefficients For Different Photometric Bands (For Rv = 3.1, Fitzpatrick(1999))
# ------------------------------------------------------------------------------------------------------------------- #
filter_df = pd.read_csv(DIR_CODE + 'FILTERS.dat', sep='\s+')
filter_df = filter_df.replace('INDEF', np.nan).set_index(['FILTER', 'Marker', 'Color']).astype('float64')
filter_df = filter_df.reset_index().set_index('FILTER')
list_filters = filter_df.index.tolist()

for index, row in filter_df.iterrows():
    if len(index) == 3 and index[0:2] == 'uv':
        name = index[-1].upper()
    else:
        name = index
    if row['Offset'] > 0:
        filter_df.loc[index, 'Label'] = name + ' + ' + str(row['Offset'])
    elif row['Offset'] == 0:
        filter_df.loc[index, 'Label'] = name
    else:
        filter_df.loc[index, 'Label'] = name + ' - ' + str(abs(row['Offset']))
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
        ERROR : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print("ERROR : File '{0}' Not Found".format(text_list))
        sys.exit(1)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.images(_doprint=0)
iraf.crutil(_doprint=0)
iraf.ccdred.instrument = 'ccddb$kpno/camera.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def deredden(file_name, ebv, prefix_str='d'):
    """
    Corrects the 1-D spectrum 'file_name' for the specified value of reddening.
    Args:
        file_name       : Name of the 1-D spectrum to be dereddened
        ebv             : Color excess E(B-V) to be used for dereddening
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file    Returns:
    Returns:
        output_filename : Output file
    """
    task = iraf.noao.onedspec.deredden
    task.unlearn()

    output_filename = prefix_str + file_name
    remove_file(output_filename)
    task(input=file_name, output=output_filename, value=ebv)

    return output_filename

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Fitting BlackBody Curve
# ------------------------------------------------------------------------------------------------------------------- #

def magtoflux(mag, band):
    return 10 ** (-0.4 * (filter_df.loc[band, 'ZeroPoint'] - filter_df.loc[band, 'RLambda'] * EBV_mag +
                          float(mag) + 21.100))


def calc_flux(wave, amp, temp):
    """
    Calculates blackbody flux as a function of wavelength (um) and temperature (K).
    Args:
        wave    : Wavelength (In Angstroms)
        amp     : Amplitude of the blackbody flux
        temp    : Temperature (In Kelvin)
    Returns:
         units of erg/s/cm^2/Angstrom
    """
    return amp * blackbody_lambda(in_x=np.asarray(wave), temperature=temp).value

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Corrects Spectra For Reddening (DeRedden)
# ------------------------------------------------------------------------------------------------------------------- #
ccm = CCM89(Rv=Rv)
fpk = F99(Rv=Rv)

data_df = pd.read_csv(DIR_PHOT + file_appmag, sep='\s+').set_index('Phase').drop(['Date', 'JD'], axis=1)
data_df = data_df[[x for x in data_df.columns.values if 'Err' not in x]].replace('INDEF', np.nan)

for band in data_df.columns.values:
    data_df[band] = data_df[band].apply(lambda mag: magtoflux(mag, band))

data_df = data_df.T
data_df.index = [filter_df.loc[band, 'CentreWave'] for band in data_df.index.values]
data_df = data_df.sort_index()

wave_data = np.array(data_df.index.values)
flux_data = np.array(data_df[3.9].values)

ebv_array = [0.0, 0.2, 0.3, 0.4, 0.5]

list_ccmflux = []
list_fpkflux = []
for ebv in ebv_array:
    list_ccmflux.append(flux_data / ccm.extinguish(wave_data * u.AA, ebv))
    list_fpkflux.append(flux_data / fpk.extinguish(wave_data * u.AA, ebv))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To SED
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

ax.plot(wave_data, flux_data, ls='', marker='*', lw=1, c='r', label='Original SED')
for index, fpkflux in enumerate(list_fpkflux):
    ax.scatter(wave_data, fpkflux, marker='*', c='k', label='F99 - E(B-V)={0}'.format(ebv_array[index]))
    popt, pcov = curve_fit(calc_flux, wave_data, fpkflux, p0=[guess_amp, guess_temp])
    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/-{1:.2f}".format(popt[1], np.sqrt(np.diag(pcov)[1])))
    print ("Amp = {0:.2e}+/-{1:.2e}".format(popt[0], np.sqrt(np.diag(pcov)[0])))

    wavearr = np.linspace(np.min(wave_data), np.max(wave_data), 100)
    ax.plot(wavearr, calc_flux(wavearr, *popt), ls='-.', lw=1.5,
            label='Blackbody Fit')

# ax.set_yticklabels([])
# ax.set_ylim(0, 40)
ax.legend(frameon=False, markerscale=2, fontsize=14)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(100))
# ax.yaxis.set_major_locator(MultipleLocator(10))
# ax.yaxis.set_minor_locator(MultipleLocator(2))
ax.set_xlabel(r'Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Flux $\rm [erg\ s^{-1}\ cm^{-2}\ {\AA}^{-1}]$', fontsize=16)

ax.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
ax.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)

fig.savefig('PLOT_BlackbodyFit.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #
