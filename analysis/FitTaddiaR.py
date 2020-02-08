#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx---------------FITTING TADDIA MODEL FOR A STRIPPED ENVELOPE SUPERNOVAE--------------xxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2017iro'
file_applc = 'AppLC_2017iro.dat'
file_bollc = 'BolLC_2017iro.dat'
absmag_solar = 4.83
lum_solar = 3.828e33
jdpirooffset = 2458050.00
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Fitting And Conversions
# ------------------------------------------------------------------------------------------------------------------- #

def lumtomag(lum):
    return absmag_solar - 2.5 * np.log10(lum / lum_solar)


def pirofunc(t, a, t0):
    return a + 0.78 * np.log10(t - t0)


def latedecayfunc(t, y0, m, t0):
    return y0 + m * (t - t0)


def gaussfunc(t, t0, g0, sigma0):
    return g0 * np.exp(-((t - t0) ** 2) / (2 * (sigma0 ** 2)))


def exprisefunc(t, tau, theta):
    return 1 - np.exp((tau - t) / theta)


def tadfunc(t, y0, m, t0, g0, sigma0, tau, theta):
    return (y0 + m * (t - t0) + g0 * np.exp(-((t - t0) ** 2) / (2 * (sigma0 ** 2)))) / 1 - np.exp((tau - t) / theta)

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
print ("Epoch of Explosion: {0:>10.2f}".format(est_dateexp))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Observational R-Band Light Curve
# Fit Taddia et al. (2018) Model
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(file_applc, sep='\s+', engine='python')
data_df = data_df.replace('INDEF', np.nan).astype('float64').round(2)
data_df['Phase'] = (data_df['JD'] - jdmax).round(2)
data_df = data_df[['JD', 'Phase', 'R', 'RErr']].dropna()
data_df['R'] = data_df['R'] - rmagmax

# y0,  m,   t0,  g0, sigma0, tau, theta
# 1, 0.02, jdmax, 5, 20,    jdmax, 10
opt, cov = curve_fit(tadfunc, data_df['JD'], data_df['R'], p0=[1, 0.02, jdmax, 5, 20, jdmax, 10])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The Observed Bolometric Light Curve
# Fit Piro & Nakar (2013) Power Law Model
# ------------------------------------------------------------------------------------------------------------------- #
bol_df = pd.read_csv(file_bollc, sep='\s+', engine='python')
bol_df['Phase'] = bol_df['JD'].apply(lambda x: x - jdpirooffset)
bol_df = bol_df[['Phase', 'JD', 'Lum', 'LogLum']]
bol_df = bol_df[bol_df['JD'] <= bol_df.loc[bol_df['Lum'].idxmax(), 'JD']]

opt2, cov2 = curve_fit(pirofunc, bol_df['Phase'], bol_df['LogLum'], p0=[40, 0])

display_text("Piro & Nakar (2013) Model")
print ("Rising Rate in Bolometric Light Curve, A = {0:.4f}+/-{1:.4f}".format(opt2[0], cov2[0, 0]))
print ("Epoch of Explosion, t0 = {0:.2f}+/-{1:.2f}\n".format(opt2[1] + jdpirooffset, cov2[1, 1]))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Power Law Fit To The Pre-Maximum Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

xaxis = np.linspace(bol_df['Phase'].min(), bol_df['Phase'].max(), 100)
ax.scatter(bol_df['Phase'] - opt2[1], bol_df['LogLum'], marker='*', c='k', s=50, label='Observed Data')
ax.plot(xaxis - opt2[1], pirofunc(xaxis, *opt2), linestyle='--', c='r', label='Model Fit')

ax.legend(markerscale=2, fontsize=12)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_major_locator(MultipleLocator(0.05))
ax.yaxis.set_minor_locator(MultipleLocator(0.01))
ax.xaxis.set_major_locator(MultipleLocator(2))
ax.xaxis.set_minor_locator(MultipleLocator(0.5))
ax.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
ax.set_ylabel('Log [Bolometric Luminosity]', fontsize=16)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)

fig.savefig('PLOT_FitPiroNakar.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots 'Observational_Data' And The 'Best_Fit' Taddia Model
# ------------------------------------------------------------------------------------------------------------------- #
fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(9, 12), sharex=True)

ax1.scatter(data_df['Phase'], data_df['R'], c='k', marker='*', s=50, label='Observed Data')
ax1.plot(data_df['Phase'], tadfunc(data_df['JD'], *opt), c='orange', lw=2, label='Best Fit')

xaxis = np.linspace(opt[2], data_df['JD'].max(), 1000)
ax1.plot(xaxis - jdmax, exprisefunc(xaxis, *[opt[i] for i in [5, 6]]), c='r', lw=2, ls='--', label='Exponential Rise')
ax1.plot(xaxis - jdmax, gaussfunc(xaxis, *[opt[i] for i in [2, 3, 4]]), c='g', lw=2, ls='-.', label='Gaussian Peak')
ax1.plot(xaxis - jdmax, latedecayfunc(xaxis, *opt[0:3]), c='b', lw=2, ls=':', label='Linear Decay')

ax1.axvline(opt[2] - jdmax, color='k', ls='--')
ax1.text(opt[2] - jdmax + 1, 1.25, r'Explosion Epoch', rotation=90, fontsize=12)

handles, labels = ax1.get_legend_handles_labels()
handles = [handles[4], handles[0], handles[1], handles[2], handles[3]]
labels = [labels[4], labels[0], labels[1], labels[2], labels[3]]
ax1.legend(handles, labels, markerscale=2, fontsize=12)

ax1.invert_yaxis()
ax1.set_ylim(4.5, -1)
ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_major_locator(MultipleLocator(1))
ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
ax1.xaxis.set_major_locator(MultipleLocator(20))
ax1.xaxis.set_minor_locator(MultipleLocator(4))
ax1.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax1.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)
ax1.set_ylabel(r'$\rm R-R_{Max}$ [mag]', fontsize=16)

ax2.scatter(data_df['Phase'], data_df['R'] - tadfunc(data_df['JD'], *opt), c='k', marker='^', label='_nolegend_')
ax2.axvline(opt[2] - jdmax, c='k', ls='--')

ax2.set_ylim(-0.15, 0.15)
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.02))
ax2.xaxis.set_major_locator(MultipleLocator(20))
ax2.xaxis.set_minor_locator(MultipleLocator(4))
ax2.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax2.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

ax2.set_xlabel('Time Since R-Band Maximum [Days]', fontsize=16)
ax2.set_ylabel(r'Residuals [Mag]', fontsize=16)

fig.subplots_adjust(hspace=0.01)
fig.savefig('PLOT_FitTaddia.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #
