#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxx-----------------OBTAIN SED AND FIT BLACKBODY FROM PHOTOMETRY---------------xxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from dust_extinction.dust_extinction import CCM89, F99
from astropy.modeling.blackbody import blackbody_lambda
from astropy.convolution import convolve, Gaussian1DKernel
from scipy.interpolate import CubicSpline, InterpolatedUnivariateSpline
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
fmt_flt = '{0:>7.3f}'
fmt_exp = '{0:>7.4e}'
Rv = 3.1
solar_rad = 6.957e10

lower_lim = 3600
upper_lim = 9000
guess_amp = 1e-20
guess_temp = 12000
wave_plot = np.arange(1000, 9001, 1)
wave_bol = np.arange(1600, 9001, 1)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Details Of SN In Study
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017hcc'
JD_keyword = 'JD'
EBV_mag = 0.029
EBV_err = 0.001
dist_val = 73.0
dist_err = 2.5
redshift = 0.0168
date_explosion = 2458027.9
date_maximum = 2458054.5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Files And Directories
# ------------------------------------------------------------------------------------------------------------------- #
file_swift = '2017hcc_SWIFT.dat'
file_spec = '2017hcc_27d.dat'

DIR_PHOT = "/home/avinash/Supernovae_Data/2017hcc/"
DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
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
# Functions For Manipulating Pandas DataFrame Containing Data Of Other Well Studied SNe
# ------------------------------------------------------------------------------------------------------------------- #

def calc_radius(lum, temp):
        return (lum / (4 * np.pi * 5.67e-5 * (temp ** 4))) ** 0.5


def calc_objlum(flux):
    val = float(flux) * 4 * np.pi * (3.086e24 ** 2)
    lum = fmt_exp.format(val * dist_val ** 2)
    lumerr = fmt_exp.format(val * ((dist_val + dist_err) ** 2 - dist_val ** 2))

    return float(lum), float(lumerr)


def calc_magflux(mag, err, band):
    mag = float(mag)
    err = float(err)
    zp = filter_df.loc[band, 'ZeroPoint']
    rlambda = filter_df.loc[band, 'RLambda']

    distmod_mag = 5 * np.log10(dist_val * 1e6) - 5
    distmod_err = 5 * np.log10((dist_val + dist_err) * 1e6) - 5 - distmod_mag

    absmag = fmt_flt.format(mag - rlambda * EBV_mag - distmod_mag)
    abserr = fmt_flt.format((err ** 2 + (rlambda * EBV_err) ** 2 + distmod_err ** 2) ** 0.5)

    flux = float(fmt_exp.format(10 ** (-0.4 * (mag - rlambda * EBV_mag + zp + 21.100))))
    fluxerr = fmt_exp.format(abs(flux - 10 ** (-0.4 * (mag + err - rlambda * EBV_mag + zp + 21.100))))

    return float(absmag), float(abserr), float(flux), float(fluxerr)


def multicol_to_fluxdf(input_df):
    """
    Converts a column-wise magnitude Pandas DataFrame to a row-wise Pandas DataFrame with
    the flux values and the absolute magnitudes.
    Args:
        input_df    : Input Pandas DataFrame
    Returns:
        output_df   : Output Pandas DataFrame
    """
    input_df = input_df.set_index('JD').drop('Phase', axis=1)
    data_arr = input_df.as_matrix()

    size = data_arr.shape
    list_jd = np.repeat(input_df.index.values, (size[1] / 2))
    list_filters = [x for x in input_df.columns.values if 'Err' not in x]
    data_arr = np.reshape(data_arr, [size[0] * size[1] / 2, 2])

    output_df = pd.DataFrame(data_arr, index=list_jd, columns=['FMAG', 'FERR'])
    output_df.index.name = 'JD'
    output_df = output_df.reset_index(drop=False)

    output_df['FILTER'] = list_filters * size[0]
#     output_df['Phase'] = output_df['JD'] - date_maximum

    output_df = output_df.replace('INDEF', np.nan).dropna(axis=0, how='any')
    output_df = output_df[['JD', 'FILTER', 'FMAG', 'FERR']].reset_index(drop=True)
    output_df['ALambda'] = output_df['FILTER'].apply(lambda x: filter_df.loc[x, 'RLambda'] * EBV_mag)

    for index, band in output_df['FILTER'].items():
        magflux = calc_magflux(mag=output_df.loc[index, 'FMAG'], err=output_df.loc[index, 'FERR'],
                               band=band)
        output_df.loc[index, 'AbsMag'] = magflux[0]
        output_df.loc[index, 'AbsErr'] = magflux[1]
        output_df.loc[index, 'Flux'] = magflux[2]
        output_df.loc[index, 'FluxErr'] = magflux[3]

    return output_df


def calc_boldf(name, input_df, flux='Flux', fluxerr='FluxErr', plot=False):
    """
    Creates a Pandas DataFrame with Bolometric Fluxes.
    Args:
        name        : Name of the SNe whose bolometric magnitudes are to be computed
        input_df    : Input Pandas DataFrame containing individual band fluxes
        flux        : Name of the Flux column in the Pandas DataFrame
        fluxerr     : Name of the Flux Error column in the Pandas DataFrame
        plot        : Whether the spline fits to the fluxes should be plotted
    Returns:
        bol_df      : Output Pandas DataFrame containing bolometric fluxes
    """
    input_df = input_df.set_index('JD')

    dict_val = {}
    for index, row in input_df.iterrows():
        if index not in dict_val.keys():
            dict_val[index] = {}
        if row['FILTER'] not in dict_val[index]:
            dict_val[index][row['FILTER']] = []
            dict_val[index][row['FILTER'] + 'Err'] = []

        dict_val[index][row['FILTER']].append(row[flux])
        dict_val[index][row['FILTER'] + 'Err'].append(row[fluxerr])

    for (day, dict_jd) in dict_val.items():
        for (band, list_flux) in dict_jd.items():
            if 'Err' not in band:
                dict_val[day][band] = np.mean(list_flux)
            else:
                dict_val[day][band] = np.sqrt(np.sum(np.square(list_flux)))

    dict_flux = {}
    for (day, dict_date) in dict_val.items():
        if len(dict_date) > 2:
            if day not in dict_flux.keys():
                dict_flux[day] = {}
            for (band, flux) in dict_date.items():
                if 'Err' not in band:
                    dict_flux[day][filter_df.loc[band, 'CentreWave']] = flux
                else:
                    dict_flux[day][str(filter_df.loc[band.rstrip('Err'), 'CentreWave']) + 'Err'] = flux

    flux_df = pd.DataFrame(dict_flux).T
    flux_df.index.name = 'JD'
    flux_df = flux_df.interpolate(method='polynomial', order=1, limit=1).T

    dict_bolflux = {}
    for jd in flux_df.columns.values:
        series = flux_df[jd].dropna()

        mag = series.loc[[x for x in list(series.index) if type(x) != str]].drop(2228.1, axis=0)
        err = series.loc[[x for x in list(series.index) if type(x) == str]].drop('2228.1Err', axis=0)

        spline = CubicSpline(mag.index.values.tolist(), mag.values.tolist(), bc_type='natural', extrapolate='True')
        spline3 = InterpolatedUnivariateSpline(mag.index.values.tolist(), mag.values.tolist(), k=2)
        popt, pcov = curve_fit(calc_bbflux, mag.index.values.tolist(), mag.values.tolist(),
                               sigma=err.values.tolist(), p0=[1e-15, 15000])

        print ("JD = {0:.1f}".format(jd))
        print ("Best-Fit Temp [{0:.1f}]= {1:.2f}+/-{2:.2f}".format(jd - date_maximum, popt[1],
                                                                   np.sqrt(np.diag(pcov)[1])))

        wave_data = wave_bol
        flux_data = calc_bbflux(wave_data, *popt)
        flux_data[flux_data < 0] = 0
        netflux = np.trapz(flux_data, wave_data)

        lum, lumerr = calc_objlum(netflux)
        radius = calc_radius(lum, popt[1])
        dict_bolflux[jd] = {}
        dict_bolflux[jd]['Phase'] = jd - date_maximum
        dict_bolflux[jd]['Flux'] = netflux
        dict_bolflux[jd]['Lum'] = lum
        dict_bolflux[jd]['LumErr'] = lumerr
        dict_bolflux[jd]['Temp'] = popt[1]
        dict_bolflux[jd]['TempErr'] = np.sqrt(np.diag(pcov)[1])
        dict_bolflux[jd]['Rad'] = radius
        dict_bolflux[jd]['RadErr'] = (radius / 2) * ((lumerr / lum) - (4 * np.sqrt(np.diag(pcov)[1]) / popt[1]))

        if plot:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111)

            ax.loglog(mag.index.values, mag.values, 'ko', label='Observed Data')
            ax.loglog(wave_plot, spline(wave_plot), 'r-', label='CubicSpline Fit')
            ax.loglog(wave_plot, spline3(wave_plot), 'g-.', label='Interpolated UnivariateSpline Fit')
            ax.loglog(wave_plot, calc_bbflux(wave_plot, *popt), 'k--', label='Blackbody Fit')
            ax.axvline(1600, ls='--', c='k')
            ax.axvline(9000, ls='--', c='k')

            ax.legend()
            ax.grid(which='both')
            ax.set_ylabel(r'Apparent Flux [$\rm erg\ s^{-1}\ cm^{-2}\ \AA^{-1}$]')
            ax.set_xlabel(r'Wavelength [$\rm \AA$]')
            ax.set_title("Temperature = {0:.2f}+/- {1:.2f}".format(popt[1], np.sqrt(np.diag(pcov)[1])))

            plt.show()
            plt.close(fig)

    if len(dict_bolflux) != 0:
        bol_df = pd.DataFrame(dict_bolflux).T
        bol_df.index.name = 'JD'
        bol_df = bol_df.reset_index()
    else:
        bol_df = pd.DataFrame()

    return bol_df

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Fitting BlackBody Curve
# ------------------------------------------------------------------------------------------------------------------- #

def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


def fit_func(wave, temp):
    wave = np.asarray(wave)
    flux_bb = calc_bbflux(wave, temp)
    return flux_bb / np.mean(flux_bb)


def calc_bbflux(wave, amp, temp):
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


def fit_specarr((wave_data, flux_data), smooth=True, sp=1, write=False):
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

    popt, pcov = curve_fit(calc_bbflux, wave_data, flux_data, p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], np.sqrt(pcov[1, 1])))

    data_df = pd.DataFrame()
    data_df['Flux'] = flux_data
    data_df['Wavelength'] = wave_data
    data_df['BBFitFlux'] = calc_bbflux(data_df['Wavelength'], *popt)
    if write:
        data_df.to_csv(name_SN + '_BlackBodyFit.dat', sep=' ', index=False, header=True)

    return data_df, popt, pcov


def fit_specdat(file_name, smooth=True, sp=2, clip_wav=(3510, 9150)):
    """
    Fits a blackbody (planckian) to a '.dat' file containg 1-D spectra.
    Args:
        file_name   : 1-D Spectra .dat file
        smooth      : Whether or not to smooth the spectra
        sp          : Smoothing parameter to be used if smooth=True
        clip_wav    : Limits at which Wavelength needs to be clipped
    Returns:
        data_df     : Pandas DataFrame containing 1-D spectra
        popt        : Optimal fit parameters
        pcov        : Covariance of the fit parameters
    """
    data_df = pd.read_csv(file_name, sep='\s+', header=None, names=['Wavelength', 'Flux'])
    if smooth:
        data_df['Flux'] = convolve(data_df['Flux'].tolist(), Gaussian1DKernel(int(sp)))

    popt, pcov = curve_fit(calc_bbflux, data_df['Wavelength'].tolist(), data_df['Flux'].tolist(),
                           p0=[guess_amp, guess_temp])

    print ("\nBest-Fit Parameters:")
    print ("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], np.sqrt(pcov[1, 1])))

    data_df['BBFitFlux'] = calc_bbflux(data_df['Wavelength'], *popt)
    data_df.to_csv(name_SN + '_BlackBodyFit.dat', sep=' ', index=False, header=True)

    return data_df, popt, pcov

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Reads SWIFT Data To Be Fit For Determining Bolometric LC
# Modifies DataFrame For Temperature & Radius To Be Plotted
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_PHOT + file_swift, sep='\s+', comment='#')
data_df = data_df.replace('INDEF', np.nan).drop('Date', axis=1).astype('float64').dropna()
data_df['Phase'] = data_df['JD'] - date_maximum
fbolm_df = calc_boldf(name, multicol_to_fluxdf(data_df))

temprad_df = fbolm_df.copy()
temprad_df[['Temp', 'TempErr']] = temprad_df[['Temp', 'TempErr']] / 1000.
temprad_df[['Rad', 'RadErr']] = temprad_df[['Rad', 'RadErr']] / (1000 * solar_rad)

for col in ['Flux', 'Lum', 'LumErr']:
    temprad_df[col] = temprad_df[col].apply(lambda x: float("{0:.2e}".format(float(x))))

for col in ['Phase', 'Temp', 'TempErr', 'Rad', 'RadErr']:
    temprad_df[col] = temprad_df[col].apply(lambda x: float("{0:.2f}".format(float(x))))

temprad_df.to_csv('2017hcc_TempRad.dat', sep=' ', header=True, index=False)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Corrects Spectra For Reddening (DeRedden)
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_PHOT + file_spec, sep='\s+', names=['Wave', 'Flux'])
wave_data = np.array(data_df['Wave'])
flux_data = np.array(data_df['Flux'])

ccm = CCM89(Rv=Rv)
fpk = F99(Rv=Rv)
ccmflux_data = flux_data / ccm.extinguish(wave_data * u.AA, EBV_mag)
fpkflux_data = flux_data / fpk.extinguish(wave_data * u.AA, EBV_mag)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig_bol = plt.figure(figsize=(8, 8))
ax_bol = fig_bol.add_subplot(111)

ax_bol.semilogy(fbolm_df['Phase'], fbolm_df['Lum'], ls='-', color='k', marker='o', markerfacecolor='None',
                markeredgewidth=1, markersize=7, alpha=0.8, label=name_SN)
ax_bol.semilogy(fbolm_df['Phase'], fbolm_df['Lum'], ls='-', color='k', marker='o', markerfacecolor='None',
                markeredgewidth=1, markersize=3, alpha=0.8, label='_nolegend_')

ax_bol.set_ylim(4.5e43, 1.1e44)
ax_bol.set_xlim(-2, 50)
ax_bol.legend(fontsize=12, markerscale=2, loc=1, frameon=False)

ax_bol.yaxis.set_ticks_position('both')
ax_bol.xaxis.set_ticks_position('both')
ax_bol.xaxis.set_major_locator(MultipleLocator(10))
ax_bol.xaxis.set_minor_locator(MultipleLocator(1))
ax_bol.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=14)
ax_bol.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=14)
ax_bol.set_xlabel('Time Since Maximum [Days]', fontsize=16)
ax_bol.set_ylabel(r'Quasi-Bolometric Luminosity [$\rm erg\ s^{-1}$]', fontsize=16)

fig_bol.savefig('PLOT_BolometricLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_bol)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Temperature And Radius Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig_temp = plt.figure(figsize=(8, 8))
ax_temp = fig_temp.add_subplot(111)
ax_rad = ax_temp.twinx()

ax_temp.plot(temprad_df['Phase'], temprad_df['Temp'], ls='', lw=0.8, c='r', marker='*', ms=11, label='Temperature')
ax_temp.errorbar(temprad_df['Phase'], temprad_df['Temp'], yerr=temprad_df['TempErr'], c='grey', marker='o', capsize=3,
                 capthick=1, ms=2, alpha=1, elinewidth=1, label='_nolegend_')

ax_rad.plot(temprad_df['Phase'], temprad_df['Rad'], ls='', lw=0.8, c='blue', marker='o', ms=9,
            label='Radius [BB Fit]')
ax_rad.errorbar(temprad_df['Phase'], temprad_df['Rad'], yerr=temprad_df['RadErr'], c='grey', marker='o', capsize=3,
                capthick=1, ms=2, alpha=0.8, elinewidth=1, label='_nolegend_')

ax_temp.xaxis.set_ticks_position('both')
ax_temp.xaxis.set_major_locator(MultipleLocator(10))
ax_temp.xaxis.set_minor_locator(MultipleLocator(1))
ax_temp.yaxis.set_major_locator(MultipleLocator(2))
ax_temp.yaxis.set_minor_locator(MultipleLocator(0.4))
ax_rad.yaxis.set_major_locator(MultipleLocator(8))
ax_rad.yaxis.set_minor_locator(MultipleLocator(0.8))
ax_temp.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_temp.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)
ax_rad.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=15)
ax_rad.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=15)

ax_temp.set_xlim(-2, 50)
ax_temp.set_ylim(3.9, 18.1)
ax_temp.axvline(2458081.12 - date_maximum, ls='--', lw=1, c='k')
ax_temp.text(2458081.12 - date_maximum - 1.5, 17.5, s='Spectral Epoch', fontsize=12, rotation=90)
ax_temp.legend(fontsize=14, markerscale=1.5, loc=2)
ax_rad.legend(fontsize=14, markerscale=1.5, loc=4)

ax_temp.set_xlabel('Time Since Maximum [Days]', fontsize=16)
ax_temp.set_ylabel(r'Temperature [$\rm \times 1000\ K$]', fontsize=16)
ax_rad.set_ylabel(r'Radius [$\rm \times 1000\ R_{\odot}$]', fontsize=16)

fig_temp.savefig('PLOT_TempRadEvolution.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig_temp)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit The Blackbody Function To The DeReddened Spectrum
# Scale The Spectra
# ------------------------------------------------------------------------------------------------------------------- #
_, popt, pcov = fit_specarr((wave_data, fpkflux_data), write=True)
bbflux_data = calc_bbflux(wave_data, *popt)

# list_filters = ['U', 'B', 'V', 'R', 'I']
# list_centrewav = [filter_df.loc[band, 'CentreWave'] for band in list_filters]
list_centrewav = [3700, 4200, 4700, 5200, 5700, 6200, 6700, 7200, 7700, 8200, 8700]

specflux = []
bbflux = []
for wave in list_centrewav:
    specflux.append(fpkflux_data[find_nearest(wave_data, wave)])
    bbflux.append(bbflux_data[find_nearest(wave_data, wave)])

scale = [x / specflux[idx] for idx, x in enumerate(bbflux)]
spline = CubicSpline(list_centrewav, scale, bc_type='natural', extrapolate='True')

bbfit_df = pd.DataFrame()
bbfit_df['Wave'] = wave_data
bbfit_df['Flux'] = fpkflux_data * spline(wave_data)
bbfit_df.to_csv('2017hcc_27d_BBFit.dat', sep=' ', index=False, header=False)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Flux-Calibrated Spectra
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

# ax.plot(wave_data, flux_data, ls='-', lw=0.8, c='k', label='Original Spectra')
# ax.plot(wave_data, ccmflux_data, ls='-.', lw=0.8, c='b', label='CCM89')
ax.plot(wave_data, fpkflux_data * spline(wave_data), 'r--', lw=0.8, label='Dereddened Spectra [F99]')
ax.plot(wave_data, bbflux_data, 'k-', lw=1, label='Blackbody Fit')

ax.set_yticklabels([])
ax.legend(frameon=False, markerscale=10, fontsize=14)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(MultipleLocator(1000))
ax.xaxis.set_minor_locator(MultipleLocator(100))
ax.yaxis.set_major_locator(MultipleLocator(5e-14))
ax.yaxis.set_minor_locator(MultipleLocator(1e-14))
ax.set_xlabel(r'Wavelength [$\rm \AA$]', fontsize=16)
ax.set_ylabel(r'Flux $\rm [erg\ s^{-1}\ cm^{-2}\ {\AA}^{-1}]$', fontsize=16)
ax.set_title("Temp = {0:.2f}+/- {1:.2f}".format(popt[1], np.sqrt(pcov[1, 1])))
ax.tick_params(axis='both', which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(axis='both', which='minor', direction='in', length=4, width=0.8, labelsize=14)

fig.savefig('PLOT_BBFluxCalibSpec.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
