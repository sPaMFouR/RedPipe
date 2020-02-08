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
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import specutils.io.read_fits as spec
from matplotlib.ticker import MultipleLocator
from astropy.modeling.blackbody import blackbody_lambda
from astropy.convolution import convolve, Gaussian1DKernel
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
guess_amp = 1e-20
guess_temp = 16000
name_SN = '2017hcc'
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


def fit_specfits(file_name, smooth=True, sp=2):
    """
    Fits a blackbody (planckian) to a FITS file containing 1-D spectra.
    Args:
        file_name   : Name of the 1-D spectrum FITS file to fit
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    wave_data, flux_data = read_1dspec(file_name)
    if smooth:
        flux_data = convolve(flux_data, Gaussian1DKernel(int(sp)))

    popt, pcov = curve_fit(calc_flux, wave_data, flux_data, p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], pcov[1, 1]))

    data_df = pd.DataFrame()
    data_df['Flux'] = flux_data
    data_df['Wavelength'] = wave_data
    data_df['BBFitFlux'] = calc_flux(data_df['Wavelength'], *popt)
    data_df.to_csv(name_SN + '_BlackBodyFit.dat', sep=' ', index=False, header=True)

    return data_df, popt, pcov


def fit_specdat(file_name, smooth=True, sp=2):
    """
    Fits a blackbody (planckian) to a '.dat' file containing 1-D spectra.
    Args:
        file_name   : Name of the 1-D spectrum .dat file to be fit
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    data_df = pd.read_csv(file_name, sep='\s+', header=None, names=['Wavelength', 'Flux'])
    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    popt, pcov = curve_fit(calc_flux, data_df['Wavelength'].tolist(), data_df['Flux'].tolist(),
                           p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], pcov[1, 1]))

    data_df['BBFitFlux'] = calc_flux(data_df['Wavelength'], *popt)
    data_df.to_csv(name_SN + '_BlackBodyFit.dat', sep=' ', index=False, header=True)

    return data_df, popt, pcov


def plot_fit(file_name, fits=True):
    """
    Args:

    """
    if fits:
        data_df, popt, pcov = fit_specfits(file_name)
    else:
        data_df, popt, pcov = fit_specdat(file_name)

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

    fig.savefig('PLOT_BlackbodyFit.pdf', format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit The Blackbody Function To All Spectra And Plot The Data
# ------------------------------------------------------------------------------------------------------------------- #
file_name = '2017hcc_27d.txt'
plot_fit(file_name, fits=False)
# ------------------------------------------------------------------------------------------------------------------- #
