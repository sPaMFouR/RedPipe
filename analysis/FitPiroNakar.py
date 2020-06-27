#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx----------------FITTING TADDIA MODEL FOR STRIPPED ENVELOPE SUPERNOVAE---------------xxxxxxxxxxxxxxxxx #
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
file_applc = '2017iro_TempSubAppLC.dat'
file_bollc = '2017iro_TempSubBolLC.dat'
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
    """
    Piro & Nakar (2013)
    """
    return a + 0.78 * np.log10(t - t0)


def caofunc(t, a, b, t0):
    """
    Cao et al. (2013)
    """
    return a + b * np.log10(t - t0)


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
# Function For Setting Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj, xticks=(50, 5), yticks=(1, 0.1), grid=True, fs=14, ms=1.5):
    """
    Sets plot parameters to the axes object 'ax_obj'.
    Args:
        ax_obj      : Axes object to be used for plotting and setting plot parameters
        xticks      : X-Axis Major and Minor tick intervals
        yticks      : Y-Axis Major and Minor tick intervals
        grid        : Boolean stating whether to enable Grid in the plot
        fs          : Font of the labels in the legend
        ms          : Scaled marker size to be used in the legend
    Returns:
    """
    if grid:
        ax_obj.grid(True, which='major', ls='--', lw=1)

    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(xticks[0]))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(xticks[1]))
    ax_obj.yaxis.set_major_locator(MultipleLocator(yticks[0]))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(yticks[1]))
    ax_obj.tick_params(axis='both', which='major', direction='in', width=1.6, length=9, color='k',
                       labelcolor='k', labelsize=fs)
    ax_obj.tick_params(axis='both', which='minor', direction='in', width=1.0, length=5, color='k',
                       labelcolor='k', labelsize=fs)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculating Date Of Explosion From R-Band light Curve (Dessart et al. 2016)
# ------------------------------------------------------------------------------------------------------------------- #
datar = pd.read_csv('OUTPUT_InterpSNMag_R', sep='\s+')
datar['Phase'] = datar['JD'] - datar['JD'].min()

rmagmax = datar['R'].min()
jdmax = datar.loc[datar['R'] == rmagmax, 'JD'].item()
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
bol_df = pd.read_csv(file_bollc, sep='\s+')
bol_df['LogLumErr'] = bol_df['LumErr'] / (np.log(10) * bol_df['Lum'])
bol_df['Phase'] = bol_df['JD'].apply(lambda x: x - jdpirooffset)
bol_df = bol_df[['Phase', 'JD', 'Lum', 'LogLum', 'LogLumErr']]
bol_df = bol_df[bol_df['JD'] <= bol_df.loc[bol_df['Lum'].idxmax(), 'JD']]

opt2, cov2 = curve_fit(pirofunc, bol_df['Phase'], bol_df['LogLum'], sigma=bol_df['LogLumErr'], p0=[40, 0])
opt3, cov3 = curve_fit(caofunc, bol_df['Phase'], bol_df['LogLum'], sigma=bol_df['LogLumErr'], p0=[41, 1, 0])
std2 = np.sqrt(np.diag(cov2))
std3 = np.sqrt(np.diag(cov3))

display_text("Piro & Nakar (2013) Model")
print ("Rising Rate in Bolometric Light Curve, A = {0:.4f}+/-{1:.4f}".format(opt2[0], std2[0]))
print ("Epoch of Explosion, t0 = {0:.2f}+/-{1:.2f}\n".format(opt2[1] + jdpirooffset, std2[1]))

display_text("Cao et al. (2013) Model - Free Power Law Index")
print ("Rising Rate in Bolometric Light Curve, A = {0:.4f}+/-{1:.4f}".format(opt3[0], std3[0]))
print ("Power Law Index, B = {0:.2f}+/-{1:.2f}".format(opt3[1], std3[1]))
print ("Epoch of Explosion, t0 = {0:.2f}+/-{1:.2f}\n".format(opt3[2] + jdpirooffset, std3[2]))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Power Law Fit To The Pre-Maximum Bolometric Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

xaxis = np.linspace(0, bol_df['Phase'].max(), 100)
# ax.plot(bol_df['Phase'] - opt2[1], bol_df['LogLum'], marker='*', c='k', ls='', ms=15, label='Observed Data')
ax.plot(xaxis - opt3[2], pirofunc(xaxis, *opt2), linestyle='--', c='r', label='Piro & Nakar Fit')
ax.errorbar(bol_df['Phase'] - opt3[2], bol_df['LogLum'], yerr=bol_df['LogLumErr'], marker='', c='dimgrey', ls='',
            capsize=3, capthick=0.5, label='_nolegend_')
ax.plot(bol_df['Phase'] - opt3[2], bol_df['LogLum'], marker='*', c='k', ls='', ms=15, label='Observed Data')
ax.plot(xaxis - opt3[2], caofunc(xaxis, *opt3), linestyle='-', c='navy', label='Piro & Nakar (Free Power Law) Fit')

ax.set_xlim(-1, 10)
ax.set_ylim(42.3, 42.75)
ax.legend(markerscale=1.5, fontsize=14)
set_plotparams(ax, xticks=(2, 0.5), yticks=(0.1, 0.01))
ax.set_ylabel('Log [Bolometric Luminosity]', fontsize=16)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)

fig.savefig('PLOT_FitPiroNakar.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #
