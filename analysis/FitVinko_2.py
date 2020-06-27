#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxx----------------FITTING VINKO MODEL FOR A STRIPPED ENVELOPE SUPERNOVA--------------xxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

sns.set_style('ticks')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables (& File Names)
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
fit_offset = 15
lum_solar = 3.828e33
mag_solar = 4.83
jd_explosion = 2458080.7

file_name = 'BolLC_2017iro.dat'
file_interp = 'InterpBolLC_2017iro.dat'
file_modA = '2017iro_VinkoA.csv'
file_modB = '2017iro_VinkoB.csv'
file_modC1 = '2017iro_VinkoC1.csv'
file_modC2 = '2017iro_VinkoC2.csv'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Calculating Energy Production
# ------------------------------------------------------------------------------------------------------------------- #

def totenergy(mni):
    return 1.885e50 * mni


def lumtomag(lum):
    return -2.5 * np.log10(lum / lum_solar) + mag_solar


def energyrate(t, mni):
    return lumtomag((6.45e43 * np.exp(-t / 8.8) + 1.45e43 * np.exp(-t / 111.3)) * mni)

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
# fig.savefig('PLOT_InterpBolLC.pdf', format='pdf', dpi=2000, bbox_inches='tight')
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
jd_max = outbol_df.loc[outbol_df['Lum'] == outbol_df['Lum'].max(), 'JD'].item()
print("Epoch of Bolometric Maximum: {0:>10.2f}".format(jd_max))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Observational Bolometric Light Curve Data
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(file_name, sep='\s+', engine='python')

data_df['PhaseMax'] = (data_df['JD'] - jd_max).round(2)
data_df['PhaseExp'] = (data_df['JD'] - jd_explosion).round(2)

data_df['AbsMag'] = data_df['Lum'].apply(lambda x: lumtomag(x))
data_df['LogLum'] = data_df['Lum'].apply(lambda x: np.log10(x))
data_df['Lum'] = data_df['Lum'].apply(lambda x: x * 1e-42)
data_df = data_df[['PhaseMax', 'PhaseExp', 'AbsMag', 'Lum', 'LogLum']]
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Observational_Data & Best Fit From Different Sub-Models Of The Vinko Model
# ------------------------------------------------------------------------------------------------------------------- #
data_modA = pd.read_csv(file_modA, names=['Time', 'Lum'], dtype='float64')
data_modB = pd.read_csv(file_modB, names=['Time', 'Lum'], dtype='float64')
data_modC1 = pd.read_csv(file_modC1, names=['Time', 'Lum'], dtype='float64')
data_modC2 = pd.read_csv(file_modC2, names=['Time', 'Lum'], dtype='float64')

for df in [data_modA, data_modB, data_modC1, data_modC2]:
    df['PhaseExp'] = df['Time'].apply(lambda x: x / 86400. + jd_max - jd_explosion)
    df['AbsMag'] = df['Lum'].apply(lambda x: lumtomag(x))
    df['LogLum'] = df['Lum'].apply(lambda x: np.log10(x))
    df['Lum'] = df['Lum'].apply(lambda x: x * 1e-42)

fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

ax.scatter(data_df['PhaseExp'], data_df['AbsMag'], marker='*', s=100, color='tomato', label='Observed Data')
ax.plot(data_modA['PhaseExp'], data_modA['AbsMag'], linestyle='-', color='k', label='Model A, B, C1, C2')
# ax.plot(data_modB['Phase'], data_modB['Lum'], linestyle='--', color='r', label='Model B')
# ax.plot(data_modC1['Phase'], data_modC1['Lum'], linestyle='-.', color='g', label='Model C1')
# ax.plot(data_modC2['Phase'], data_modC2['Lum'], linestyle=':', color='orange', label='Model C2')

phasearr = np.linspace(data_df['PhaseExp'].min(), data_df['PhaseExp'].max(), 500)
ax.plot(phasearr, energyrate(phasearr, 0.06), ls='-.', lw=1, c='b', label=r'$\rm M({^{56}Ni})=0.06\ M_{\odot}$')
ax.plot(phasearr, energyrate(phasearr, 0.09), ls='--', lw=1, c='orange', label=r'$\rm M({^{56}Ni})=0.09\ M_{\odot}$')
ax.plot(phasearr, energyrate(phasearr, 0.125), ls=':', lw=1, c='brown', label=r'$\rm M({^{56}Ni})=0.125\ M_{\odot}$')
ax.axvline(jd_max - jd_explosion, linestyle='-', linewidth=1, color='gray')
ax.text(jd_max - jd_explosion - 2, -16, s='Maximum Epoch', fontsize=10, rotation='vertical')

ax.legend(markerscale=2, frameon=False, fontsize=13)
ax.invert_yaxis()
ax.set_xlim(-2, 100)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

ax.set_ylabel('Bolometric Magnitude [mag]', fontsize=16)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)

fig.savefig('PLOT_FitVinko.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
