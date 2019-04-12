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
# Read Interpolated Bolometric Light Curve Data
# Fit The Bolometric Light Curve With Taddia et al. (2018) Model
# ------------------------------------------------------------------------------------------------------------------- #
interpbol_df = pd.read_csv('OUTPUT_InterpSNMag_Bol', sep='\s+')
interpbol_df['BolMag'] = interpbol_df['Lum'].apply(lumtomag)
bolmagmax = interpbol_df['BolMag'].min()
boljdmax = interpbol_df.loc[interpbol_df['BolMag'] == interpbol_df['BolMag'].min(), 'JD'].values[0]

bol_df = pd.read_csv(file_bollc, sep='\s+', engine='python')
bol_df['BolMag'] = bol_df['Lum'].apply(lambda x: lumtomag(x) - bolmagmax)
opt, cov = curve_fit(tadfunc, bol_df['JD'], bol_df['BolMag'], p0=[1, 0.015, boljdmax - 20, 10, 10, boljdmax, 20])

boljdarrexp = np.round(np.arange(opt[2], bol_df['JD'].max(), 0.1), 1)
boljdarr = np.round(np.arange(bol_df['JD'].min(), bol_df['JD'].max(), 0.1), 1)
bolfitarr = tadfunc(boljdarr, *opt) + bolmagmax

fitbolmagmax = np.min(bolfitarr)
fitboljdmax = boljdarr[np.where(bolfitarr == np.min(bolfitarr))][0]
fitbolmag15 = bolfitarr[np.where(boljdarr == fitboljdmax + 15)][0]
fitbolmag40 = bolfitarr[np.where(boljdarr == fitboljdmax + 40)][0]

boljdpremax = boljdarr[np.where(boljdarr <= fitboljdmax)]
boljdpostmax = boljdarr[np.where(boljdarr >= fitboljdmax)]
bolfitpremax = bolfitarr[np.where(boljdarr <= fitboljdmax)]
bolfitpostmax = bolfitarr[np.where(boljdarr >= fitboljdmax)]
fitboljdl25 = boljdpremax[np.abs(bolfitpremax - fitbolmagmax - 0.25).argmin()]
fitboljdh25 = boljdpostmax[np.abs(bolfitpostmax - fitbolmagmax - 0.25).argmin()]

fitbolerrmax = (abs(bolfitarr[np.where(np.roll(boljdarr, -5) == fitboljdmax)][0] - fitbolmagmax) +
                abs(bolfitarr[np.where(np.roll(boljdarr, 5) == fitboljdmax)][0] - fitbolmagmax)) / 2.
fitbolerr15 = (abs(bolfitarr[np.where(np.roll(boljdarr, -5) == fitboljdmax + 15)][0] - fitbolmag15) +
               abs(bolfitarr[np.where(np.roll(boljdarr, 5) == fitboljdmax + 15)][0] - fitbolmag15)) / 2.
fitbolerr40 = (abs(bolfitarr[np.where(np.roll(boljdarr, -5) == fitboljdmax + 40)][0] - fitbolmag40) +
               abs(bolfitarr[np.where(np.roll(boljdarr, 5) == fitboljdmax + 40)][0] - fitbolmag40)) / 2.

boldeltam15 = fitbolmag15 - fitbolmagmax
boldeltam40 = fitbolmag40 - fitbolmagmax
boldeltad25 = fitboljdh25 - fitboljdl25

boldeltam15err = (fitbolerrmax ** 2 + fitbolerr15 ** 2) ** 0.5
boldeltam40err = (fitbolerrmax ** 2 + fitbolerr40 ** 2) ** 0.5

print ("#" + "-" * 35 + "#")
print ("Epoch of Bolometric Maximum: {0:>9.1f}+/-{1:.1f}".format(fitboljdmax, 0.1))
print ("Epoch of Explosion: {0:>9.1f}+/-{1:.1f}".format(opt[2], np.sqrt(np.diag(cov)[2])))
print ("Maximum Bolometric Magnitude: {0:>5.2f}+/-{1:4.2f}".format(fitbolmagmax, fitbolerrmax))
print (r"$\rm \delta m_{15}$" + ": {0:4.2f}+/-{1:4.2f}".format(boldeltam15, boldeltam15err))
print (r"$\rm \delta m_{40}$" + ": {0:4.2f}+/-{1:4.2f}".format(boldeltam40, boldeltam40err))
print (r"$\rm \delta d_{0.25}$" + ": {0:4.2f}+/-0.14".format(boldeltad25))
print ("Late Time Decay Rate: {0:>6.4f}+/-{1:6.4f}".format(opt[1], cov[1, 1]))
print ("#" + "-" * 35 + "#")

fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(9, 12), sharex=True)

ax1.scatter(bol_df['JD'] - fitboljdmax, bol_df['BolMag'], marker='*', s=50, c='k', label='Observed Data')
ax1.plot(boljdarr - fitboljdmax, tadfunc(boljdarr, *opt), c='orange', lw=2, label='Best Fit')
ax1.plot(boljdarrexp - fitboljdmax, exprisefunc(boljdarrexp, *[opt[i] for i in [5, 6]]), c='r', lw=2, ls='--',
         label='Exponential Rise')
ax1.plot(boljdarrexp - fitboljdmax, gaussfunc(boljdarrexp, *[opt[i] for i in [2, 3, 4]]), c='g', lw=2, ls='-.',
         label='Gaussian Peak')
ax1.plot(boljdarrexp - fitboljdmax, latedecayfunc(boljdarrexp, *opt[0:3]), c='b', lw=2, ls=':',
         label='Linear Decay')

ax1.axvline(opt[2] - fitboljdmax, c='k', ls='--')
ax1.text(opt[2] - fitboljdmax + 1, 1.25, r'Explosion Epoch', rotation=90, fontsize=12)

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
ax1.set_ylabel(r'$\rm M-M_{Max}$ [mag]', fontsize=16)

ax2.scatter(bol_df['JD'] - fitboljdmax, bol_df['BolMag'] - tadfunc(bol_df['JD'], *opt), marker='^', c='k', label=None)
ax2.axvline(opt[2] - fitboljdmax, ls='--', c='k')

ax2.set_ylim(-0.15, 0.15)
ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.yaxis.set_major_locator(MultipleLocator(0.1))
ax2.yaxis.set_minor_locator(MultipleLocator(0.02))
ax2.xaxis.set_major_locator(MultipleLocator(20))
ax2.xaxis.set_minor_locator(MultipleLocator(4))
ax2.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
ax2.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

ax2.set_xlabel('Time Since Bolometric Maximum [Days]', fontsize=16)
ax2.set_ylabel(r'Residuals [Mag]', fontsize=16)

fig.subplots_adjust(hspace=0.01)
fig.savefig('PLOT_FitTaddiaBol.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Apparent Light Curve Data
# Fit The Apparent Light Curves With Taddia et al. (2018) Model
# ------------------------------------------------------------------------------------------------------------------- #
dict_guess = {'U': [2.5, 0.0015, 2458080, 20, 10, 2458100, 10],
              'B': [1.5, 0.005, 2458080, 20, 10, 2458100, 10],
              'V': [1, 0.015, 2458080, 20, 10, 2458100, 10],
              'R': [0.5, 0.015, 2458080, 10, 10, 2458100, 10],
              'I': [0.5, 0.015, 2458080, 10, 10, 2458100, 10]}

app_df = pd.read_csv(file_applc, sep='\s+', engine='python')
app_df = app_df.replace('INDEF', np.nan).astype('float64').round(2)
list_bands = [val for val in app_df.columns.values[0:] if 'Err' not in val and val not in ['JD', 'Phase']]

for band in list_bands:
    data_interpband = pd.read_csv('OUTPUT_InterpSNMag_' + band, sep='\s+')
    bandmagmax = data_interpband[band].min()
    bandjdmax = data_interpband.loc[data_interpband[band] == data_interpband[band].min(), 'JD'].values[0]

    band_df = app_df[['JD', band, band + 'Err']].copy().dropna()
    band_df[band] = band_df[band] - bandmagmax

#     opt, cov = curve_fit(tadfunc, band_df['JD'], band_df[band], sigma=band_df[band + 'Err'],
#                          p0=[1, 0.015, bandjdmax - 20, 10, 10, bandjdmax, 20])
    opt, cov = curve_fit(tadfunc, band_df['JD'], band_df[band], sigma=band_df[band + 'Err'],
                         p0=dict_guess[band])
    print opt
    bandjdarrexp = np.round(np.arange(opt[2], band_df['JD'].max(), 0.1), 1)
    bandjdarr = np.round(np.arange(band_df['JD'].min(), band_df['JD'].max(), 0.1), 1)
    bandfitarr = tadfunc(bandjdarr, *opt) + bandmagmax

    fitbandmagmax = np.min(bandfitarr)
    fitbandjdmax = bandjdarr[np.where(bandfitarr == np.min(bandfitarr))][0]
    fitbandmag15 = bandfitarr[np.where(bandjdarr == fitbandjdmax + 15)][0]
    fitbandmag40 = bandfitarr[np.where(bandjdarr == fitbandjdmax + 40)][0]

    bandjdpremax = bandjdarr[np.where(bandjdarr <= fitbandjdmax)]
    bandjdpostmax = bandjdarr[np.where(bandjdarr >= fitbandjdmax)]
    bandfitpremax = bandfitarr[np.where(bandjdarr <= fitbandjdmax)]
    bandfitpostmax = bandfitarr[np.where(bandjdarr >= fitbandjdmax)]
    fitbandjdl25 = bandjdpremax[np.abs(bandfitpremax - fitbandmagmax - 0.25).argmin()]
    fitbandjdh25 = bandjdpostmax[np.abs(bandfitpostmax - fitbandmagmax - 0.25).argmin()]

    fitbanderrmax = (abs(bandfitarr[np.where(np.roll(bandjdarr, -5) == fitbandjdmax)][0] - fitbandmagmax) +
                     abs(bandfitarr[np.where(np.roll(bandjdarr, 5) == fitbandjdmax)][0] - fitbandmagmax)) / 2.
    fitbanderr15 = (abs(bandfitarr[np.where(np.roll(bandjdarr, -5) == fitbandjdmax + 15)][0] - fitbandmag15) +
                    abs(bandfitarr[np.where(np.roll(bandjdarr, 5) == fitbandjdmax + 15)][0] - fitbandmag15)) / 2.
    fitbanderr40 = (abs(bandfitarr[np.where(np.roll(bandjdarr, -5) == fitbandjdmax + 40)][0] - fitbandmag40) +
                    abs(bandfitarr[np.where(np.roll(bandjdarr, 5) == fitbandjdmax + 40)][0] - fitbandmag40)) / 2.

    banddeltam15 = fitbandmag15 - fitbandmagmax
    banddeltam40 = fitbandmag40 - fitbandmagmax
    banddeltad25 = fitbandjdh25 - fitbandjdl25
    banddeltam15err = (fitbanderrmax ** 2 + fitbanderr15 ** 2) ** 0.5
    banddeltam40err = (fitbanderrmax ** 2 + fitbanderr40 ** 2) ** 0.5

    print ("#" + "-" * 45 + "#")
    print ("Epoch of {0}-Band Maximum: {1:>10.2f}+/-0.1".format(band, fitbandjdmax))
    print ("Epoch of Explosion: {0:>10.2f}+/-{1:4.2f}".format(opt[2], np.sqrt(np.diag(cov)[2])))
    print ("Maximum {0}-Band Magnitude: {1:>5.2f}+/-{2:4.2f}".format(band, fitbandmagmax, fitbanderrmax))
    print (r"$\rm \delta m_{15}$" + ": {0:4.2f}+/-{1:4.2f}".format(banddeltam15, banddeltam15err))
    print (r"$\rm \delta m_{40}$" + ": {0:4.2f}+/-{1:4.2f}".format(banddeltam40, banddeltam40err))
    print (r"$\rm \delta d_{0.25}$" + ": {0:4.2f}+/-0.14".format(banddeltad25))
    print ("Late Time Decay Rate : {0:>6.4f}+/-{1:6.4f}".format(opt[1], cov[1, 1]))
    print ("#" + "-" * 45 + "#")

    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(9, 12), sharex=True)

    ax1.scatter(band_df['JD'] - fitbandjdmax, band_df[band], marker='*', s=50, c='k', label='Observed Data')
    ax1.plot(bandjdarr - fitbandjdmax, tadfunc(bandjdarr, *opt), c='orange', lw=2, label='Best Fit')
    ax1.plot(bandjdarrexp - fitbandjdmax, exprisefunc(bandjdarrexp, *[opt[i] for i in [5, 6]]), c='r', lw=2, ls='--',
             label='Exponential Rise')
    ax1.plot(bandjdarrexp - fitbandjdmax, gaussfunc(bandjdarrexp, *[opt[i] for i in [2, 3, 4]]), c='g', lw=2, ls='-.',
             label='Gaussian Peak')
    ax1.plot(bandjdarrexp - fitbandjdmax, latedecayfunc(bandjdarrexp, *opt[0:3]), c='b', lw=2, ls=':',
             label='Linear Decay')

    ax1.axvline(opt[2] - bandjdmax, ls='--', c='k')
    ax1.text(opt[2] - bandjdmax + 1, 1.25, r'Explosion Epoch', rotation=90, fontsize=12)

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
    ax1.set_ylabel('{0}-{0}'.format(band) + r'$\rm_{Max}$ [mag]', fontsize=16)

    ax2.scatter(band_df['JD'] - fitbandjdmax, band_df[band] - tadfunc(band_df['JD'], *opt), marker='^', c='k',
                label=None)
    ax2.axvline(opt[2] - bandjdmax, ls='--', c='k')

    ax2.set_ylim(-0.15, 0.15)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.yaxis.set_major_locator(MultipleLocator(0.1))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.02))
    ax2.xaxis.set_major_locator(MultipleLocator(20))
    ax2.xaxis.set_minor_locator(MultipleLocator(4))
    ax2.tick_params(which='major', direction='in', length=8, width=1.4, labelsize=14)
    ax2.tick_params(which='minor', direction='in', length=4, width=0.7, labelsize=14)

    ax2.set_xlabel('Time Since ' + band + '-Band Maximum [Days]', fontsize=16)
    ax2.set_ylabel(r'Residuals [Mag]', fontsize=16)

    fig.subplots_adjust(hspace=0.01)
    fig.savefig('PLOT_FitTaddia{0}.pdf'.format(band), format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------- #
