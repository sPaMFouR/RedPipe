#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxx-------------------LIGHT CURVE MODELLING-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from jdcal import jd2gcal, gcal2jd
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Input Data Files
# ------------------------------------------------------------------------------------------------------------------- #
obsvel_file = 'OUTPUT_PhotVel'
obsphot_file = 'OUTPUT_NetSNMag'
obsbol_file = 'OUTPUT_DateWiseSNBolFlux'
input_file = 'ASASSN-14dq.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_Model = "/home/avinash/Dropbox/ModelSNe/IIP/"
DIR_SNe = "/home/avinash/Supernovae_Data/2016gfy/Photometry/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
precision = 3
JD_offset = 2400000
lum_solar = 3.828e33
light_speed = 299792.458

vel_file = 'velo.dat'
phot_file = 'phot.dat'
velmod_file = 'model_velo.dat'
photmod_file = 'model_phot.dat'
bol_file = 'lbol_err.dat'

input_cols = ['Ref', 'DataSet', 'Band', 'JD', 'Val', 'Err']
vel_cols = ['JD', 'Vel', 'Err', 'FILTER', 'No.', 'dy', 'Radius', 'Weight']
phot_cols = ['JD', 'Mag', 'Err', 'FILTER', 'No.', 'dy', 'Radius', 'Tau', 'Weight']

velmod_cols = ['JD', 'Vel']
photmod_cols = ['JD', 0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
bol_cols = ['JD', 'Lum', 'Err1', 'Err2', 'TErr']

dict_bands = {'U': 0, 'B': 1, 'V': 2, 'R': 3, 'I': 4, 'J': 5, 'H': 6, 'K': 7, 'g': 9, 'r': 10, 'i': 11, 'z': 12, 
              'uvw2': 13, 'uvm2': 14, 'uvw1': 15, 'u': 16, 'b': 17, 'v': 18, 'Z': 19, 'Y': 20}
dict_codes = dict(zip(dict_bands.values(), dict_bands.keys()))

dict_markers = {0: [+3.0, 'o', 'k', ' + 3.0'], 1: [+1.5, 'D', 'b', ' + 1.5'], 2: [+0.5, 's', 'g', ' + 0.5'],
                3: [+0.0, '^', 'r', '      '], 4: [-1.0, 'p', 'm', ' - 1.0'], 9: [+1.0, '*', 'y', ' + 1.0'],
                10: [-1.0, 'P', 'c', ' - 1.0'], 11: [-2.0, 'X', 'indigo', ' - 2.0']}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Convert Julian Date Into Calendar Date In String Format & Vice Versa
# ------------------------------------------------------------------------------------------------------------------- #

def jd_to_cald(julian_day):
    """
    Converts julian day into calendar day in string format.
    Args:
        julian_day  : Julian day value to be converted to calendar day
    Returns:
        cal_date    : Calendar date corresponding to input julian day
    """
    time_tuple = jd2gcal(epoch, julian_day - epoch)
    cal_date = date(*time_tuple[0:3]).strftime("%Y-%m-%d")
    
    return cal_date


def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split("-")
    jd_tuple = gcal2jd(date_comp[0], date_comp[1], date_comp[2])
    julian_day = jd_tuple[0] + jd_tuple[1]

    return julian_day

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Plotting
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj):
    """
    Sets the plot parameters for displaying the Light Curves and Expansion Velocity Evolution.
    Args:
        ax_obj   : Axes object to be used for plotting and setting plot parameters
    Returns:
        None
    """
    ax_obj.legend()
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(50))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(10))
    ax_obj.set_xlabel('Time Since Explosion [Days]', fontsize=12)
    ax_obj.tick_params(which='both', direction='in', width=1, labelsize=12)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Generate An Input File To The Light Curve Fitting Code
# ------------------------------------------------------------------------------------------------------------------- #
obsvel_df = pd.read_csv(obsvel_file, sep='\s+', engine='python')[['Date', '5169', '5169Err']]
obsphot_df = pd.read_csv(obsphot_file, sep='\s+', engine='python')[['FILTER', 'JD', 'FMAG', 'FERR']]
obsbol_df = pd.read_csv(obsbol_file, sep='\s+', engine='python')[['Date', 'JD', 'Lum', 'Err']]

obsvel_df.columns = [['Date', 'Val', 'Err']]
obsphot_df.columns = [['FILTER', 'JD', 'Val', 'Err']]
obsbol_df = obsbol_df[(obsbol_df.Date != '2014-11-06') & (obsbol_df.Date != '2015-06-05') & 
                    (obsbol_df.Date != '2014-11-05') & (obsbol_df.Date != '2014-12-02') & 
                    (obsbol_df.Date != '2015-05-23')]

obsvel_df['Ref'] = 0
obsvel_df['DataSet'] = 0
obsvel_df['Band'] = 0
obsvel_df = obsvel_df.replace('INDEF', np.nan, regex=True)
obsvel_df['JD'] = obsvel_df['Date'].apply(lambda x: cald_to_jd(x) - JD_offset)
obsvel_df['Val'] = obsvel_df['Val'].apply(lambda x: (((5169 - float(x)) / 5169) * light_speed))
obsvel_df['Err'] = obsvel_df['Err'].apply(lambda x: (((float(x))/ 5169) * light_speed))
obsvel_df[['JD', 'Val', 'Err']] = obsvel_df[['JD', 'Val', 'Err']].round(precision)
obsvel_df = obsvel_df[input_cols].dropna(axis=0, how='any')

obsphot_df['Ref'] = 0
obsphot_df['DataSet'] = 1
obsphot_df['Band'] = obsphot_df['FILTER'].apply(lambda x: dict_bands[x])
obsphot_df['JD'] = obsphot_df['JD'].apply(lambda x: x - JD_offset)
obsphot_df[['JD', 'Val', 'Err']] = obsphot_df[['JD', 'Val', 'Err']].round(precision)
obsphot_df = obsphot_df[input_cols]

obsinput_df = pd.concat([obsvel_df, obsphot_df])
obsinput_df.to_csv(input_file, sep=' ', header=None, index=None)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Runs The Code & Reads I/O Data Files
# ------------------------------------------------------------------------------------------------------------------- #
os.system(DIR_Model + "./fit_single.exe " + input_file)

input_df = pd.read_csv(input_file, sep='\s+', header=None, names=input_cols, engine='python')
vel_df = pd.read_csv(vel_file, sep='\s+', header=None, names=vel_cols, engine='python')
phot_df = pd.read_csv(phot_file, sep='\s+', header=None, names=phot_cols, engine='python')

velmod_df = pd.read_csv(velmod_file, sep='\s+', header=None, names=velmod_cols, engine='python')
photmod_df = pd.read_csv(photmod_file, sep='\s+', header=None, names=photmod_cols, comment='#', engine='python')
bol_df = pd.read_csv(bol_file, sep='\s+', header=None, names=bol_cols, engine='python')
bol_df[['Lum', 'Err1', 'Err2', 'TErr']] = bol_df[['Lum', 'Err1', 'Err2', 'TErr']].multiply(lum_solar, axis='index')

JDvel_max = vel_df['JD'].max()
JDvel_min = vel_df['JD'].min()
JDphot_max = phot_df['JD'].max()
JDphot_min = phot_df['JD'].min()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To The Expansion Velocities
# ------------------------------------------------------------------------------------------------------------------- #
fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111)

velinp_df = input_df[input_df['DataSet'] == 0]
velinp_df = velinp_df[velinp_df['JD'] < JDvel_max]
velmod_df = velmod_df[velmod_df['JD'] < JDvel_max]

ax1.errorbar(vel_df['JD'], vel_df['Vel'], yerr=vel_df['Err'], c='k', marker='*', linestyle='', linewidth=0.5, 
            capsize=2, capthick=1, label=None)

ax1.plot(velmod_df['JD'], velmod_df['Vel'], marker=None, c='k', linestyle='--', label='_nolegend_')
    
set_plotparams(ax1)
ax1.set_xlim(JDvel_min - 20, JDvel_max + 20)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_major_locator(MultipleLocator(2000))
ax1.yaxis.set_minor_locator(MultipleLocator(500))
ax1.set_ylabel(r'Velocity [$\rm km\ s^{-1}$]', fontsize=12)

plt.show()
plt.close(fig1)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To The Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111)

photinp_df = input_df[input_df['DataSet'] == 1]
photinp_df = photinp_df[photinp_df['JD'] < JDphot_max]
photmod_df = photmod_df[photmod_df['JD'] < JDphot_max]

for band, band_df in photinp_df.groupby('Band'):
    ax2.scatter(band_df['JD'], band_df['Val'] + dict_markers[band][0], marker=dict_markers[band][1], 
               c=dict_markers[band][2], label=dict_codes[band] + str(dict_markers[band][3]))
    ax2.errorbar(band_df['JD'], band_df['Val'] + dict_markers[band][0], yerr=band_df['Err'], c=dict_markers[band][2], 
                fmt='', linestyle='', linewidth=0.5, capsize=2, capthick=1, label=None)

for band in photinp_df['Band'].unique():
    ax2.plot(photmod_df['JD'], photmod_df[band] + dict_markers[band][0], marker=None, c='k', linestyle='--', 
            label='_nolegend_')

set_plotparams(ax2)
ax2.invert_yaxis()
ax2.set_xlim(JDphot_min - 20, JDphot_max + 20)
ax2.yaxis.set_major_locator(MultipleLocator(2))
ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
ax2.set_ylabel('Apparent Magnitude [mag]', fontsize=12)

plt.show()
plt.close(fig1)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Fit To The Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig3 = plt.figure(figsize=(10, 7))
ax3 = fig3.add_subplot(111)

obsbol_df['JD'] = obsbol_df['JD'].apply(lambda x: x - JD_offset)
ax3.semilogy(obsbol_df['JD'], obsbol_df['Lum'], c='b', marker='*', linestyle='', linewidth=0.5, label='Observed')
ax3.semilogy(bol_df['JD'], bol_df['Lum'], c='k', linestyle='--', linewidth=0.5, label='Fit')
    
set_plotparams(ax3)
ax3.set_xlim(JDphot_min - 20, JDphot_max + 20)
ax3.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
ax3.set_ylabel(r'Quasi-Bolometric Luminosity [$\rm erg\ s^{-1}$]', fontsize=12)

plt.show()
plt.close(fig3)
# ------------------------------------------------------------------------------------------------------------------- #

