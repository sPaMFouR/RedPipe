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
guess_amp = 10
guess_temp = 10000
norm_factor = 1e-16
lower_lim = 3650
upper_lim = 9100

Rv = 3.1
EBV_mag = 0.21
EBV_err = 0.11

JD_keyword = 'JD'
name_SN = '2016gfy'
file_earlyspec = '../rfz_2016-09-13_2016gfy.dat'

DIR_SPEC = '/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/Test/'
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


def read_1dspec(file_name):
    """
    Reads 1-D Spectra from a FITS file and returns wavelength and flux arrays.
    Args:
        file_name    : FITS file from which data has to be extracted
    Returns:
        wave_array   : Array containing wavelength values extracted from the 1-D spectra
        flux_array   : Array containing flux values extracted from the 1-D spectra
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


def fit_specarr((wave_data, flux_data), smooth=True, sp=3, write=False):
    """
    Fits a blackbody (planckian) to a wave_arr, flux_arr containing 1-D spectrum data.
    Args:
        wave_data   : Wavelength array of the 1-D spectrum to be fit
        flux_data   : Flux array of the 1-D spectrum to be fit
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
        write       : Whether to write the best fit data to a file 
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    comb_data = zip(wave_data, flux_data)                                            
    wave_data = [wave for wave, _ in sorted(comb_data, key=lambda x: x[0]) if upper_lim >= wave >= lower_lim]         
    flux_data = [flux for wave, flux in sorted(comb_data, key=lambda x: x[0]) if upper_lim >= wave >= lower_lim]             

    if smooth:
        flux_data = convolve(flux_data, Gaussian1DKernel(int(sp)))

    popt, pcov = curve_fit(calc_flux, wave_data, flux_data, p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], pcov[1, 1]))

    data_df = pd.DataFrame()
    data_df['Flux'] = flux_data
    data_df['Wavelength'] = wave_data
    data_df['BBFitFlux'] = calc_flux(data_df['Wavelength'], *popt)
    if write:
        data_df.to_csv('BBFit_SpecArr.dat', sep=' ', index=False, header=True, na_rep='INDEF')

    return data_df, popt, pcov


def fit_specfits(file_name, smooth=True, sp=3, write=False):
    """
    Fits a blackbody (planckian) to a FITS file containing 1-D spectra.
    Args:
        file_name   : Name of the 1-D spectrum FITS file to be fit
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
        write       : Whether to write the best fit data to a file 
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    wave_data, flux_data = read_1dspec(file_name)    
    comb_data = zip(wave_data, flux_data)                                            
    wave_data = [wave for wave, _ in sorted(comb_data, key=lambda x: x[0]) if upper_lim >= wave >= lower_lim]         
    flux_data = [flux for wave, flux in sorted(comb_data, key=lambda x: x[0]) if upper_lim >= wave >= lower_lim]             

    if smooth:
        flux_data = convolve(flux_data, Gaussian1DKernel(int(sp)))

    popt, pcov = curve_fit(calc_flux, wave_data, flux_data, p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], pcov[1, 1]))
    print ("Amp = {0:.2f}+/- {1:.2f}".format(popt[0], pcov[0, 0]))

    data_df = pd.DataFrame()
    data_df['Flux'] = flux_data
    data_df['Wavelength'] = wave_data
    data_df['BBFitFlux'] = calc_flux(data_df['Wavelength'], *popt)
    
    if write:
        data_df.to_csv('BBFit_' + file_name.split('.')[0] + '.dat', sep=' ', index=False, 
                       header=True, na_rep='INDEF')

    return data_df, popt, pcov


def fit_specdat(file_name, smooth=True, sp=3, write=False):
    """
    Fits a blackbody (planckian) to a '.dat' file containing 1-D spectra.
    Args:
        file_name   : Name of the 1-D spectrum .dat file to be fit
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
        write       : Whether to write the best fit data to a file 
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    data_df = pd.read_csv(file_name, sep='\s+', header=None, names=['Wavelength', 'Flux'])
    data_df = data_df[(data_df['Wavelength'] >= lower_lim) & (data_df['Wavelength'] <= upper_lim)]

    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))
        
    popt, pcov = curve_fit(calc_flux, data_df['Wavelength'].tolist(), data_df['Flux'].tolist(),
                           p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], pcov[1, 1]))

    data_df['BBFitFlux'] = calc_flux(data_df['Wavelength'], *popt)
    
    if write:
        data_df.to_csv('BBFit_' + file_name.split('.')[0] + '.dat', sep=' ', index=False, 
                       header=True, na_rep='INDEF')

    return data_df, popt, pcov


def plot_fit(file_name, fits=False, smooth=True, sp=3):
    """
    Args:
        file_name   : Name of the 1-D spectrum file to be fit
        fits        : Whether or not the file is a FITS file
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    if fits:
        data_df, popt, pcov = fit_specfits(file_name, smooth=smooth, sp=sp)
    else:
        data_df, popt, pcov = fit_specdat(file_name, smooth=smooth, sp=sp)

    temp_fit = popt[1]
    temp_err = pcov[1, 1]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    ax.plot(data_df['Wavelength'], data_df['Flux'], ls='steps-mid', lw=1, c='r', label=name_SN)
    ax.plot(data_df['Wavelength'], calc_flux(data_df['Wavelength'], *popt), c='g', ls='-.', lw=1.5,
            label='Blackbody Fit')

    ax.set_yticklabels([])
    ax.legend(frameon=False, markerscale=3, fontsize=14)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(0.2e-14))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05e-14))
    ax.set_xlabel(r'Wavelength [$\rm \AA$]', fontsize=16)
    ax.set_ylabel(r'Flux $\rm [erg\ s^{-1}\ cm^{-2}\ {\AA}^{-1}]$', fontsize=16)
    ax.set_title('Best Fit Temp = {0:.2f} +/- {1:.2f}'.format(popt[1], pcov[1, 1]))

    ax.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
    ax.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)

#     fig.savefig('PLOT_BlackbodyFit.pdf', format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Reads Data To Be Fit
# # ------------------------------------------------------------------------------------------------------------------- #
# ctext = 'rfz_*.fits'
# list_files = group_similar_files('', common_text=ctext)
# # ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Fit The Blackbody Function To All Spectra And Plot The Data
# # ------------------------------------------------------------------------------------------------------------------- #
# data_temp = pd.DataFrame(columns=['JD', 'Temp', 'TempErr'])

# for index, file_name in enumerate(list_files):
#     best_temp, sigma_temp = fit_specfits(file_name, sp=10)
#     data_temp.loc[index, 'JD'] = round(float(fits.getval(filename=file_name, keyword=JD_keyword)), 2)
#     data_temp.loc[index, 'Temp'] = float("{0:.1f}".format(best_temp))
#     data_temp.loc[index, 'TempErr'] = float("{0:.1f}".format(sigma_temp))

# data_temp.to_csv('OUTPUT_BlackBodyFit', sep=' ', index=False, header=True)

# for index, file_name in enumerate(list_files):
#     plot_fit(file_name, fits=False)

# # ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Corrects Spectra For Reddening (DeRedden)
# ------------------------------------------------------------------------------------------------------------------- #

ccm = CCM89(Rv=Rv)
fpk = F99(Rv=Rv)

data_df = pd.read_csv(DIR_SPEC + file_earlyspec, sep='\s+', names=['Wavelength', 'Flux'])
data_df = data_df[(data_df['Wavelength'] > lower_lim) & (data_df['Wavelength'] < upper_lim)]
data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(4)) / norm_factor

wave_data = np.array(data_df['Wavelength'])
flux_data = np.array(data_df['Flux'])

ebv_array = [0.0, 0.21, 0.3, 0.4]

list_ccmflux = []
list_fpkflux = []
for ebv in ebv_array:
    list_ccmflux.append(flux_data / ccm.extinguish(wave_data * u.AA, ebv))
    list_fpkflux.append(flux_data / fpk.extinguish(wave_data * u.AA, ebv))
#     deredden(ebv=)

# ------------------------------------------------------------------------------------------------------------------- #



# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To The Earliest Spectra
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

ax.plot(wave_data, data_df['Flux'], ls='-', lw=1, c='r', label='Original Spectra')
for index, fpkflux in enumerate(list_fpkflux):
    ax.plot(wave_data, fpkflux, ls='-.', lw=0.8, label='F99 - E(B-V)={0}'.format(ebv_array[index]))
    data_df, popt, pcov = fit_specarr((wave_data, fpkflux))
    ax.plot(data_df['Wavelength'], calc_flux(data_df['Wavelength'], *popt), ls='-.', lw=1.5,
            label='Blackbody Fit')

# ax.set_yticklabels([])
ax.set_ylim(0, 40)
ax.legend(frameon=False, markerscale=10, fontsize=14)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_minor_locator(MultipleLocator(2))
ax.set_xlabel(r'Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Flux $\rm [erg\ s^{-1}\ cm^{-2}\ {\AA}^{-1}]$', fontsize=16)

ax.tick_params(which='major', direction='in', width=1.2, length=7, labelsize=14)
ax.tick_params(which='minor', direction='in', width=0.8, length=4, labelsize=14)

fig.savefig('PLOT_BlackbodyFit.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #

