#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxx-------------FITTING ARNETT-VALENTI MODEL FOR A STRIPPED ENVELOPE SUPERNOVA--------------xxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.ticker import MultipleLocator
from lmfit import Minimizer, Parameters, report_fit
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
file_name = 'BolLC_2017iro.dat'
date_explosion = 2458083.1
fit_epoch = 2458120
guess_mni = 0.06
guess_taum = 1.5e6
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Constants
# ------------------------------------------------------------------------------------------------------------------- #
eni = 3.90e10               # Energy Production in 1s by 1g of 56Ni
eco = 6.78e9                # Energy Production in 1s by 1g of 56Co
tauni = 757728.0            # e-folding time of 56Ni
tauco = 9616320.0           # e-folding time of 56Co
mass_solar = 1.989e33       # Solar Mass
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Displaying Text
# ------------------------------------------------------------------------------------------------------------------- #

def display_text(text_to_display):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text_to_display : Text to be displayed
    Returns:
        None
    """
    print ("\n" + "# " + "-" * (12 + len(text_to_display)) + " #")
    print ("# " + "-" * 5 + " " + str(text_to_display) + " " + "-" * 5 + " #")
    print ("# " + "-" * (12 + len(text_to_display)) + " #" + "\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Estimating The Value of 'taum'
# ------------------------------------------------------------------------------------------------------------------- #
kappa = 0.06
beta = 13.8
c = 2.99792458e10
msolar = 1.99e33
vph = 9.5e8
mej = 4.0 * msolar
ek = 0.3 * (vph ** 2) * mej

est_taum = ((2 * kappa * mej) / (beta * c * vph)) ** 0.5
print("Estimated Value of Taum: {0:>5.2f}".format(est_taum))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculating Date Of Explosion From R-Band light Curve (Dessart et al. 2016)
# ------------------------------------------------------------------------------------------------------------------- #
datar = pd.read_csv('OUTPUT_InterpSNMag_R', sep='\s+')
datar['Phase'] = datar['JD'] - datar['JD'].min()

rmagmax = datar['R'].min()
jdmax = datar.loc[datar['R'] == datar['R'].min(), 'JD'].item()
rmag15 = datar.loc[datar['JD'] == jdmax + 15, 'R'].item()

delm = (rmag15 - rmagmax)
trise = 57.08 - 71.17 * delm + 32.98 * (delm ** 2)
est_dateexp = jdmax - trise
display_text("Dessart et al. (2016)")
print ("Epoch of R-Band Maximum: {0:>10.2f}".format(jdmax))
print ("Epoch of Explosion: {0:>10.2f}\n".format(est_dateexp))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Set Of Parameters
# ------------------------------------------------------------------------------------------------------------------- #
params = Parameters()
params.add('mni', value=guess_mni)
params.add('taum', value=guess_taum)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Calculating Luminosity From Valenti Model Based On Initial Values
# ------------------------------------------------------------------------------------------------------------------- #

def calc_lum(param, list_time):
    mni = param['mni'].value
    taum = param['taum'].value
    time_int = [float(time / taum) for time in list_time]

    y = float(taum / (2 * tauni))
    s = float(((taum * (tauco - tauni)) / (2 * tauco * tauni)))

    def funca(z, y):
        return 2 * z * np.exp((-2 * z * y) + (z ** 2))

    def funcb(z, y, s):
        return 2 * z * np.exp((-2 * z * y) + (2 * z * s) + (z ** 2))

    lph = np.zeros(len(list_time))
    lph_err = np.zeros(len(list_time))

    for index, time in enumerate(time_int):
        lph[index] = mass_solar * mni * np.exp(-time ** 2) * ((eni - eco) * quad(funca, 0, time, args=(y,))[0] +
                                                              eco * quad(funcb, 0, time, args=(y, s))[0])
        lph_err[index] = mass_solar * mni * np.exp(-time ** 2) * ((eni - eco) * quad(funca, 0, time, args=(y,))[1] +
                                                                  eco * quad(funcb, 0, time, args=(y, s))[1])
    return [lph, lph_err]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Defines Objective Function For Best Fit & Returns Numpy Array To Be Minimized
# ------------------------------------------------------------------------------------------------------------------- #

def func_min(param, list_time, list_flux):
    lph = calc_lum(param, list_time)[0]

    delta = np.zeros(len(list_time))
    delta_err = np.zeros(len(list_time))
    for index, flux in enumerate(list_flux):
        delta[index] = float(lph[index]) - float(flux)
        delta_err[index] = (float(lph[index]) - float(flux)) / float(flux)

    return delta

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Observational Bolometric Light Curve Data
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(file_name, sep='\s+', engine='python')
data_df['Phase'] = (data_df['JD'] - date_explosion).round(2)
data_df['Time'] = data_df['Phase'].apply(lambda x: 86400 * x).round(2)

data_df = data_df[data_df['JD'] <= fit_epoch]
data_df = data_df[['JD', 'Phase', 'Time', 'Lum']].dropna()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Minimize Using 'Least Squares' Technique
# ------------------------------------------------------------------------------------------------------------------- #
minner_object = Minimizer(func_min, params, fcn_args=(data_df['Time'], data_df['Lum']), nan_policy='omit')
result_fit = minner_object.minimize(method='leastsq', params=params, **{'maxfev': 100})
data_df['FitLum'] = data_df['Lum'] + result_fit.residual

report_fit(result_fit)
xaxis = np.linspace(data_df['Time'].min(), data_df['Time'].max(), 1000)
fit_lum, fit_lumerr = calc_lum(result_fit.params, xaxis)
datemax = xaxis[fit_lum.argmax()] / 86400.
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Observational_Data & Best_Fit From The Valenti Model
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

ax.scatter(data_df['Phase'], data_df['Lum'].apply(lambda x: x * 1e-42), marker='*', color='k', label='Observed Data')
ax.plot([time / 86400 for time in xaxis], fit_lum * 1e-42, linestyle='--', color='r', label='Best Fit')

ax.axvline(datemax, linestyle='--', color='k')
ax.text(datemax + 1, 1.5, r'Maximum Epoch', rotation=90, fontsize=12)

ax.set_xlim(0, 40)
ax.legend(frameon=False, markerscale=2, fontsize=12)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.xaxis.set_major_locator(MultipleLocator(10))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax.set_ylabel(r'Bolometric Luminosity [x$\rm 10^{42}\ erg\ s^{-1}$]', fontsize=16)

fig.savefig('PLOT_FitValenti.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #
