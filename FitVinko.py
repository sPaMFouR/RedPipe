#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxx----------------FITTING VINKO MODEL FOR A STRIPPED ENVELOPE SUPERNOVA--------------xxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables (& File Names)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
file_name = 'BolLC_2017iro.dat'
file_interp = 'InterpBolLC_2017iro.dat'
file_modA = '2017iro_VinkoA.csv'
file_modB = '2017iro_VinkoB.csv'
file_modC1 = '2017iro_VinkoC1.csv'
file_modC2 = '2017iro_VinkoC2.csv'
fit_offset = 15
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Calculating Energy Production
# ------------------------------------------------------------------------------------------------------------------- #

def energyrate(t, mni):
    return (6.45e43 * np.exp(-t / 8.8) + 1.45e43 * np.exp(-t / 111.3)) * mni * 1e-42


def totenergy(mni):
    return 1.885e50 * mni

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit A Spline To The Bolometric Light Curve And Record The Spline Fit Data In A Data File
# ------------------------------------------------------------------------------------------------------------------- #
bol_df = pd.read_csv(file_name, sep='\s+')
bol_df = bol_df[['JD', 'Lum']]

jdarray = np.arange(bol_df['JD'].min(), bol_df['JD'].max() + 0.25, 0.25)
spline = interpolate.UnivariateSpline(bol_df['JD'], bol_df['Lum'], k=3, s=0.1)
magarray = np.array([float(spline(jd)) for jd in jdarray])

# fig = plt.figure(figsize=(8, 8))
# ax = fig.add_subplot(111)
# ax.plot(jdarray, spline(jdarray), 'k-')
# ax.scatter(bol_df['JD'], bol_df['Lum'], marker='*', c='r')
# fig.savefig('PLOT_InterpBolLC.eps', format='eps', dpi=600, bbox_inches='tight')
# plt.show()
# plt.close(fig)

outbol_df = pd.DataFrame()
outbol_df['JD'] = jdarray
outbol_df['Lum'] = magarray
outbol_df = outbol_df.round(2)
outbol_df.to_csv(file_interp, sep=' ', index=False, header=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculating Date Of Explosion From The Interpolated Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
jdmax = outbol_df.loc[outbol_df['Lum'] == outbol_df['Lum'].max(), 'JD'].item()
print("Epoch of Bolometric Maximum: {0:>10.2f}".format(jdmax))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Observational Bolometric Light Curve Data
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(file_name, sep='\s+', engine='python')
# data_df = data_df[data_df['JD'] >= fit_offset + jdmax]

data_df['Phase'] = (data_df['JD'] - jdmax).round(2)
data_df['Time'] = data_df['Phase'].apply(lambda x: 86400 * x).round(2)
data_df['Lum'] = data_df['Lum'].apply(lambda x: x * 1e-42)
data_df = data_df[['Time', 'Phase', 'Lum']]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Observational_Data & Best Fit From Different Sub-Models Of The Vinko Model
# ------------------------------------------------------------------------------------------------------------------- #
data_modA = pd.read_csv(file_modA, names=['Time', 'Lum'], dtype='float64')
data_modB = pd.read_csv(file_modB, names=['Time', 'Lum'], dtype='float64')
data_modC1 = pd.read_csv(file_modC1, names=['Time', 'Lum'], dtype='float64')
data_modC2 = pd.read_csv(file_modC2, names=['Time', 'Lum'], dtype='float64')

for df in [data_modA, data_modB, data_modC1, data_modC2]:
    df['Phase'] = df['Time'].apply(lambda x: x / 86400.)
    df['LogLum'] = df['Lum'].apply(lambda x: np.log10(x))
    df['Lum'] = df['Lum'].apply(lambda x: x * 1e-42)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)

ax.scatter(data_df['Phase'], data_df['Lum'], marker='*', s=100, color='r', label='Observed Data')
ax.plot(data_modA['Phase'], data_modA['Lum'], linestyle='-', color='k', label='Model A, B, C1, C2')
# ax.plot(data_modB['Phase'], data_modB['Lum'], linestyle='--', color='r', label='Model B')
# ax.plot(data_modC1['Phase'], data_modC1['Lum'], linestyle='-.', color='g', label='Model C1')
# ax.plot(data_modC2['Phase'], data_modC2['Lum'], linestyle=':', color='orange', label='Model C2')

phasearray = np.linspace(data_df['Phase'].min(), data_df['Phase'].max(), 500)
ax.plot(phasearray, energyrate(phasearray, 0.04), linestyle='-.', color='b',
        label=r'$\rm M({^{56}Ni})=0.04\ M_{\odot}$')
ax.plot(phasearray, energyrate(phasearray, 0.05), linestyle='--', color='orange',
        label=r'$\rm M({^{56}Ni})=0.05\ M_{\odot}$')
ax.plot(phasearray, energyrate(phasearray, 0.06), linestyle=':', color='darkcyan',
        label=r'$\rm M({^{56}Ni})=0.06\ M_{\odot}$')
ax.plot(phasearray, energyrate(phasearray, 0.09), linestyle='-', color='k', label=r'$\rm M({^{56}Ni})=0.09\ M_{\odot}$')

ax.legend(markerscale=2, fontsize=14)
ax.set_ylim(0.2, 2.7)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.tick_params(which='major', direction='in', length=6, width=1, labelsize=14)
ax.tick_params(which='minor', direction='in', length=3, width=1, labelsize=14)

ax.set_ylabel(r'Bolometric Luminosity [$\rm x10^{42}\ erg\ s^{-1}$]', fontsize=16)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)

fig.savefig('PLOT_FitVinko2.eps', format='eps', dpi=1000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
