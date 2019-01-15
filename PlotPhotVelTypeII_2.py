#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxx---------------PLOT PHOTOSPHERIC VELOCITY EVOLUTION IN TYPE II SNe----------------xxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import pandas as pd
from jdcal import gcal2jd
import uncertainties as unc
import matplotlib.pyplot as plt
from matplotlib import container
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
# plt.rc('font', family='serif', serif='Times')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
file_photwave = 'OUTPUT_SpecPhotMinima'
plot_epoch = 145
date_explosion = 2457641.4
light_speed = 2.99792458e5        # In Km/s
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Paths Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
DIR_CURNT = os.getcwd()
DIR_SNe = "/home/avinash/Dropbox/SNData/IIP_Data/"
DIR_SPEC = "/home/avinash/Supernovae_Data/2016gfy/Spectroscopy/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def expfunc(t, A, t0, alpha):
    return A * ((t - t0) ** alpha)


def powerlaw(t, A, B):
    return A * (t ** B)


def redchisq(ydata, ymod, sd, n=8):
    """
    Args:
        ydata    : Observed data
        ymod     : Model data
          sd     : Uncertainties in ydata
          n      : Number of free parameters in the model
    Returns:
        RedChiSq : Reduced Chi-Square
    """
    ydata = np.array(ydata)
    ymod = np.array(ymod)
    sd = np.array(sd)

    if not sd.any():
        chisq = np.sum((ydata - ymod) ** 2)
    else:
        chisq = np.sum(((ydata - ymod) / sd) ** 2)

    nu = ydata.size - 1 - n
    redchisq = chisq / nu

    return redchisq

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Confidence Intervals
# ------------------------------------------------------------------------------------------------------------------- #

def plot_confintervals(ax_obj, optpar, covpar, xarr, func=powerlaw, fcolor='blue'):
    """
    Plots 3-Sigma Confidence Intervals in Fits of SN Parameters.
    Args:
        ax_obj  : Axes object on which the confidence interval is to be plotted
        optpar  : Optimised Parameters of the Fit
        covpar  : Covariance Parameters of the Fit
        xarr    : Array of X-Values over which confidence interval will be plotted
        func    : Function used for fitting
        fcolor  : Fill color for the confidence intervals
    Returns:
        None
    """
    if func == powerlaw:
        A, B = unc.correlated_values(optpar, covpar)
        val = powerlaw(xarr, A / 1000, B)
    elif func == expfunc:
        A, t0, alpha = unc.correlated_values(optpar, covpar)
        val = expfunc(xarr, A / 1000, t0, alpha)
    else:
        print ("ERROR: Invalid Function Type Given As Input")

    fit = unp.nominal_values(val)
    sigma = unp.std_devs(val)

    if func == expfunc:
        fit_df = pd.DataFrame(index=np.round(xarr, 1))
        fit_df.index.name = 'Phase'
        fit_df = fit_df.reset_index(drop=False)
        fit_df['Vel'] = np.round(1000 * fit, 1)
        fit_df['VelErr'] = np.round(1000 * sigma, 1)
        fit_df.to_csv('OUTPUT_InterpFeII5169Vel', sep=' ', index=False, header=True)

    fitlow = fit - 3 * sigma
    fithigh = fit + 3 * sigma

    ax_obj.plot(xarr, fit, ls='--', c='k', lw=2, alpha=0.6, label='_nolegend_')
    ax_obj.plot(xarr, fitlow, ls='-.', c='k', lw=1, alpha=0.3, label='_nolegend_')
    ax_obj.plot(xarr, fithigh, ls='-.', c='k', lw=1, alpha=0.3, label='_nolegend_')
    ax_obj.fill_between(xarr, fitlow, fithigh, facecolor=fcolor, alpha=0.3)

    return fit, sigma

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

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

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def cald_to_jd(cal_date):
    """
    Converts calendar date into julian day.
    Args:
        cal_date    : Calendar date corresponding to input julian day
    Returns:
        julian_day  : Julian day value to be converted to calendar day
    """
    date_comp = cal_date.split('-')
    jd_tuple = gcal2jd(date_comp[0], date_comp[1], str(int(float(date_comp[2])) + 1))
    julian_day = jd_tuple[0] + jd_tuple[1]

    return julian_day


def calc_meanstd(str_num):
    """
    Calculates mean and standard deviation of a list of values.
    Args:
        str_num : Input string containing numbers separated by a comma
    Returns:
        mean    : Mean of the numbers in the input string
        std     : Standard deviation of the numbers in the input string
    """
    if type(str_num) != float:
        list_val = [float(val) for val in str_num.split(',')]
        return np.mean(list_val), np.std(list_val)
    else:
        return np.nan, np.nan

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj, ylims=(1.6, 16)):
    """
    This function sets the plot parameters for Photospheric Velocity evolution plots using "ax_obj" of the plot.
    Args:
        ax_obj   : Axes object to which the plot parameters are to be applied
        ylims    : Upper limit on Y-axis
    Returns:
        None
    """
    ax_obj.set_ylim(ylims[0], ylims[1])
    ax_obj.yaxis.set_ticks_position('both')
    ax_obj.xaxis.set_ticks_position('both')
    ax_obj.xaxis.set_major_locator(MultipleLocator(30))
    ax_obj.xaxis.set_minor_locator(MultipleLocator(3))
    ax_obj.yaxis.set_major_locator(MultipleLocator(2))
    ax_obj.yaxis.set_minor_locator(MultipleLocator(0.2))
    ax_obj.tick_params(which='major', direction='in', width=1.4, length=8, labelsize=16)
    ax_obj.tick_params(which='minor', direction='in', width=0.7, length=4, labelsize=16)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Pandas DataFrame With Velocity Information From Different Spectral Lines
# Fit Power Law To Line Evolution Of FeII5169
# ------------------------------------------------------------------------------------------------------------------- #
data_df = pd.read_csv(DIR_SPEC + file_photwave, sep='\s+')
data_df = data_df.set_index('Date').replace('INDEF', np.nan)
wave_df = pd.DataFrame(index=data_df.index.values)
wave_df.index.name = 'Date'

for column in data_df.columns.values:
    wave_df[column] = data_df[column].apply(lambda x: calc_meanstd(x)[0])
    wave_df[column + 'Err'] = data_df[column].apply(lambda x: calc_meanstd(x)[1])

wave_df = wave_df.reset_index()
wave_df['JD'] = wave_df['Date'].apply(lambda x: cald_to_jd(x))
wave_df['Phase'] = wave_df['JD'] - date_explosion
wave_df = wave_df.set_index('Phase', drop=True)

dict_markers = {'6563': ['*', 'b', r'H$\alpha$'], '4861': ['v', 'r', r'H$\beta$'],
                '5169': ['D', 'k', r'$\rm Fe\,II$ 5169'], '5018': ['P', 'darkgreen', r'$\rm Fe\,II\ 5018$'],
                '4924': ['s', 'darkorange', r'$\rm Fe\,II$ 4924'], '4340': ['o', 'm', r'H$\gamma$'],
                '5876': ['X', 'brown', r'$\rm He\,I$ 5876']}
# 6142: ['o', 'orange', r'Ba$\rm \,II$ 6142'], 6246: ['*', 'teal', r'Sc$\rm \,II$ 6246']

vel_df = pd.DataFrame()
vel_df['Date'] = wave_df['Date'].copy()

for wave in dict_markers.keys():
    vel_df[wave] = wave_df[wave].apply(lambda x: (float(wave) - float(x)) * light_speed / float(wave)).round(0)
    vel_df[wave + 'Err'] = wave_df[wave + 'Err'].apply(lambda x: (float(x) + float(wave) / 1000) *
                                                       light_speed / float(wave)).round(0)

vel_df.round(1).to_csv('OUTPUT_PhotVel', sep=' ', index=True, header=True, na_rep='INDEF')

temp_5169 = vel_df[['5169', '5169Err']].copy().dropna()
temp_6563 = vel_df[['6563', '6563Err']].copy().dropna()

opt, cov = curve_fit(expfunc, temp_5169.index.values, temp_5169['5169'], sigma=temp_5169['5169Err'], p0=[1e5, 1, -0.5])
opt2, cov2 = curve_fit(powerlaw, temp_6563.index.values, temp_6563['6563'], sigma=temp_6563['6563Err'], p0=[1e4, -0.4])
opt3, cov3 = curve_fit(powerlaw, temp_5169.index.values, temp_5169['5169'], sigma=temp_5169['5169Err'], p0=[1e4, -0.4])
err = np.sqrt(np.diag(cov))
err2 = np.sqrt(np.diag(cov2))
err3 = np.sqrt(np.diag(cov3))

print np.round(opt, 2), np.round(np.sqrt(np.diag(cov)), 2)
print np.round(opt2, 2), np.round(np.sqrt(np.diag(cov2)), 2)
print np.round(opt3, 2), np.round(np.sqrt(np.diag(cov3)), 2)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Photospheric Velocity Evolution
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111)

for wave in ['6563', '4861', '4340', '4924', '5018', '5169', '5876']:
    temp_series = vel_df[wave].dropna()
    ax.errorbar(temp_series.index.values, temp_series.values / 1000, yerr=vel_df[wave + 'Err'].dropna() / 1000,
                ls='--', lw=1, marker=dict_markers[wave][0], ms=9, capsize=3, capthick=1,
                c=dict_markers[wave][1], alpha=0.7, label=dict_markers[wave][2])

set_plotparams(ax, ylims=(1.3, 17))
ax.set_xlim(-5, plot_epoch + 20)
ax.set_xlabel('Time Since Explosion [Days]', fontsize=16)
ax.set_ylabel(r'Velocity [$\rm \times\ 10^3\ km\ s^{-1}$]', fontsize=16)

axins = zoomed_inset_axes(ax, 1.4, loc=1)
axins.errorbar(temp_5169.index.values, temp_5169['5169'] / 1000, yerr=temp_5169['5169Err'] / 1000, ls='',
               alpha=0.7, capsize=3, ms=8, marker=dict_markers['5169'][0], capthick=1, elinewidth=1,
               c=dict_markers['5169'][1], label=dict_markers['5169'][2])

handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
ax.legend(handles, labels, fontsize=14, markerscale=1.5, frameon=False, loc=3)

set_plotparams(axins, ylims=(2.1, 7.6))
axins.set_xlim(20, 110)
axins.yaxis.set_major_locator(MultipleLocator(1))
axins.text(65, 7.1, s=r'$\rm Fe\,II\ 5169\ Exponential\ Fit$', fontsize=14)

plot_confintervals(axins, opt, cov, np.arange(20, 120, 0.5), func=expfunc, fcolor='chocolate')

fig.savefig('PLOT_PhotVel.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Photospheric Velocity Inferred From Halpha Line & FeII Line (For Type II SNe)
# ------------------------------------------------------------------------------------------------------------------- #
dict_snmark = {'1999em': ['o', 'r'], '2012aw': ['D', 'b'], '2013ej': ['s', 'g'], '1999gi': ['^', 'r'],
               '2004et': ['X', 'm'], '2005cs': ['*', 'y'], '2013ab': ['P', 'darkorange'], '2009bw': ['h', 'violet'],
               '2006bp': ['+', 'lime']}

list_halpha = group_similar_files('', DIR_SNe + 'PhotVel_Data/*Halpha.txt')
list_phot = group_similar_files('', DIR_SNe + 'PhotVel_Data/*Phot.txt')

fig2, (ax2, ax3) = plt.subplots(2, 1, figsize=(8, 14), sharex=True)

for file_name in list_halpha:
    name = file_name.split('/')[-1].split('.')[0].split('_')[0]
    if name in ['2012aw', '2013ab', '2013ej', '2004et', '1999em']:
        temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', names=['Phase', 'Vel'])
        temp_df = temp_df[temp_df['Phase'] < plot_epoch]
        ax2.plot(temp_df['Phase'], temp_df['Vel'] / 1000, linestyle='--', marker=dict_snmark[name][0], ms=9,
                 c=dict_snmark[name][1], alpha=0.7, label=name)

for file_name in list_phot:
    name = file_name.split('/')[-1].split('.')[0].split('_')[0]
    if name in ['2012aw', '2013ab', '2013ej', '2004et', '1999em']:
        temp_df = pd.read_csv(file_name, sep='\s+', header=None, comment='#', names=['Phase', 'Vel'])
        temp_df = temp_df[temp_df['Phase'] < plot_epoch]
        ax3.plot(temp_df['Phase'], temp_df['Vel'] / 1000, linestyle='--', marker=dict_snmark[name][0], ms=9,
                 c=dict_snmark[name][1], alpha=0.7, label=name)

ax2.plot(vel_df.index.values, vel_df['6563'] / 1000, marker='*', ms=13, linestyle='-', color='k', label=name_SN)
ax3.plot(vel_df.index.values, vel_df['5169'] / 1000, marker='*', ms=13, linestyle='-', color='k', label=name_SN)

set_plotparams(ax2)
set_plotparams(ax3)
ax3.set_ylim(2, 11)
ax2.set_ylim(3.5, 15.5)
ax2.set_xlim(-3, plot_epoch)
ax2.legend(fontsize=13, markerscale=1.5, loc=1, frameon=False)
plot_confintervals(ax2, opt2, cov2, np.arange(0.5, 140, 0.5), fcolor='chocolate')
plot_confintervals(ax3, opt3, cov3, np.arange(20, 120, 0.5), fcolor='chocolate')

ax2.text(15, 14.6, s=r'$\rm Power\ Law\ Exponent: {0:.3f}\pm{1:.3f}$'.format(opt2[1], err2[1]), fontsize=14)
ax3.text(50, 10.3, s=r'$\rm Power\ Law\ Exponent: {0:.3f}\pm{1:.3f}$'.format(opt3[1], err3[1]), fontsize=14)

ax2.set_ylabel(r'$\rm V_{H\alpha}\ [\times\ 10^3\ km\ s^{-1}]$', fontsize=16)
ax3.set_ylabel(r'$\rm V_{FeII}\ [\times\ 10^3\ km\ s^{-1}]$', fontsize=16)
ax3.set_xlabel('Time Since Explosion [Days]', fontsize=16)

fig2.subplots_adjust(hspace=0.01)
fig2.savefig('PLOT_CompPhotVel.pdf', format='pdf', dpi=2000, bbox_inches='tight')
plt.show()
plt.close(fig2)
# ------------------------------------------------------------------------------------------------------------------- #
