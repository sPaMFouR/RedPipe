#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxx--------------DETERMINE THE STEEPNESS PARAMETER FOR TYPE IIP SUPERNOVA---------------xxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
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
# Useful Functions For Fitting And Conversions
# ------------------------------------------------------------------------------------------------------------------- #

def line(x, m, c):
    return m * x + c


def elmfunc(t, a, b, t0, p, q):
    return a * ((t / t0) ** p / (1 + (t / t0) ** q)) + b * np.exp(-t / 111.26)


def valfunc(t, a0, tpt, w0, p0, m0):
    return (-a0 / (1 + np.exp((t - tpt) / w0))) + (p0 * t) + m0


def magtoflux(mag):
    return 10 ** (-0.4 * (mag + 21.1))


def fluxtomag(flux):
    return -2.5 * np.log10(flux) - 21.1

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Reading And Fitting V-Band Light Curve Data
# ------------------------------------------------------------------------------------------------------------------- #

def read_data(file_name):
    """
    Read V-band data from the file 'file_name' and output it to a pandas DataFrame.
    Args:
        file_name   : Name of the file from which V-band magnitudes are to be extracted
    Returns:
        data_df     : Pandas DataFrame containing V-band data
    """
    name = file_name.split('/')[-1].split('.')[0]
    data_df = pd.read_csv(file_name, sep='\s+', comment='#', engine='python')

    if 'Phase' not in data_df.columns.values:
        data_df['Phase'] = data_df['JD'] - date_explosion

    data_df = data_df[['Phase', 'V', 'VErr']].sort_values(by='Phase')
    data_df = data_df.replace('INDEF', np.nan).astype('float64').dropna()
    data_df['Flux'] = data_df['V'].apply(magtoflux)
    data_df['FluxErr'] = data_df['Flux'] - (data_df['V'] + data_df['VErr']).apply(magtoflux)
    data_df = data_df[(data_df['Phase'] > dict_epoch[name][0]) & (data_df['Phase'] < dict_epoch[name][1])]

    return data_df


def fit_data(name, data_df, xdata, epoch_tran):
    """
    Read V-band data from the file 'file_name' and output it to a pandas DataFrame.
    Args:
        name        : Name of the SNe to which the function has to be fit
        data_df     : Pandas DataFrame containing V-band data
        xdata       : Linearly spaced array over which the fit to the SN is to be examined
        epoch_tran  : Rough value of the epoch of transition since the time of the explosion
    Returns:
        mag         : List of V-band magnitudes of the Functional fit to the SN
        grad        : First derivative of the functional fit to the SN
        steepness   : Maximum slope in the transition region of the SN
        inflection  : Epoch in the transition region where the light curve has the maximum slope
    """
    if name not in ['ASASSN-14ha', '2005af', '1992ba', '2007it', '2016X']:
        opt, cov = curve_fit(elmfunc, data_df['Phase'], data_df['Flux'], sigma=data_df['FluxErr'],
                             p0=[data_df['Flux'].mean(), data_df['Flux'].mean(), epoch_tran, 2, 10])
    else:
        opt, cov = curve_fit(elmfunc, data_df['Phase'], data_df['Flux'],
                             p0=[data_df['Flux'].mean(), data_df['Flux'].mean(), epoch_tran, 2, 10])

    mag = [fluxtomag(flux) for flux in elmfunc(xdata, *opt)]
    grad = np.gradient(mag, xdata[1] - xdata[0])
    steepness = grad.max()
    inflection = xdata[np.where(grad == grad.max())][0]

    return mag, grad, steepness, inflection


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function To Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plotparams(ax_obj, xticks=(50, 10), yticks=(1, 0.1), grid=True, fs=14, ms=1.5):
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


def plot_steepness(list_sne, xlims=[30, 210], ylim=0.45, out_suffix=1):
    """
    Plots the Functional fit to the V-band light curve of type II SNe and determines the steepness parameter
    and the point of inflection.
    Args:
        list_sne   : List of SNe for which steepness is to be determined
        xlim       : X-Axis limits for plotting epoch
        ylim       : Higher Y-Axis limit for steepness
        out_suffix : Suffix to be added to the plot name
    Returns:
        None
    """
    fig, axarr = plt.subplots(2, 5, gridspec_kw={'height_ratios': [2, 1]}, figsize=(20, 10), sharex=True)

    for count, file_name in enumerate(list_sne):
        name = file_name.split('/')[-1].split('.')[0]
        data_df = read_data(file_name)
        xaxis = np.linspace(data_df['Phase'].min(), data_df['Phase'].max(), 1000)
        mag, grad, steepness, inflection = fit_data(name, data_df, xaxis, epoch_tran=dict_epoch[name][2])

        axarr[0, count].scatter(data_df['Phase'], data_df['V'], marker='*', c='k', label=name)
        axarr[0, count].plot(xaxis, mag, color='r', linestyle='--', label='_nolegend_')
        axarr[1, count].plot(xaxis, grad, color='r', label='_nolegend_')

        set_plotparams(axarr[0, count])
        set_plotparams(axarr[1, count], yticks=(0.1, 0.02))

        axarr[0, count].invert_yaxis()
        axarr[0, count].set_xlim(xlims[0], xlims[1])
        axarr[1, count].set_ylim(-0.02, ylim)
        axarr[0, count].axvline(inflection, linestyle='--', color='r', linewidth=1)
        axarr[1, count].axvline(inflection, linestyle='--', color='r', linewidth=1)
        axarr[1, count].legend(markerscale=0, markerfirst=False, fontsize=14, frameon=False, loc=3)
        axarr[1, count].text(inflection + 10, ylim / 1.5, r'$t_i$={0:5.1f} d'.format(inflection), fontsize=14)
        axarr[1, count].legend(markerscale=0, markerfirst=False, fontsize=14, frameon=False, loc=3)

        if out_suffix == 2:
            axarr[1, count].text(inflection * 0.32, ylim / 2., r'$S$={0:.3f}'.format(steepness), fontsize=14)
        elif out_suffix == 6 or out_suffix == 4:
            axarr[1, count].text(inflection * 0.40, ylim / 2., r'$S$={0:.3f}'.format(steepness), fontsize=14)
        else:
            axarr[1, count].text(inflection * 0.45, ylim / 2., r'$S$={0:.3f}'.format(steepness), fontsize=14)

        if count == 0:
            axarr[0, count].set_ylabel('Apparent Magnitude [mag]', fontsize=16)
            axarr[1, count].set_ylabel(r'$|\rm dM_V/dt|\ [mag\ d^{-1}]$', fontsize=16)
            axarr[1, count].tick_params(axis='y', which='both', direction='in', width=1, labelsize=14)
        elif count == 2:
            axarr[1, count].set_xlabel('Time Since Explosion [Days]', fontsize=16)
            axarr[1, count].tick_params(axis='y', which='both', direction='in', width=1, labelleft='off')
        else:
            axarr[1, count].tick_params(axis='y', which='both', direction='in', width=1, labelleft='off')

    fig.subplots_adjust(hspace=0, wspace=0)
    fig.savefig('PLOT_Steepness' + str(out_suffix) + '.pdf', format='', dpi=600, bbox_inches='tight')
    plt.show()
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Steepness Parameter For SNe In Study
# ------------------------------------------------------------------------------------------------------------------- #

dict_epoch = {'1992af': [50, 110, 80], '1992ba': [60, 210, 125], '1999em': [70, 170, 115], '1999gi': [70, 170, 120],
              '2002hx': [20, 140, 75], '2003gd': [50, 180, 125], '2003hn': [50, 150, 95], '2004dj': [95, 180, 125],
              '2004ej': [80, 160, 115], '2004et': [70, 220, 125], '2004fx': [70, 140, 105], '2005af': [50, 150, 110],
              '2005cs': [70, 180, 125], '2007it': [15, 160, 90], '2008gz': [80, 170, 120], '2008in': [60, 160, 110],
              '2009N': [85, 140, 110], '2009bw': [80, 190, 135], '2009ib': [100, 200, 140], '2009md': [90, 160, 120],
              '2012A': [80, 170, 105], '2012aw': [30, 300, 130], '2012ec': [60, 160, 105], '2013K': [90, 230, 130],
              '2013ab': [70, 150, 100], '2013ej': [70, 130, 100], '2013by': [60, 155, 85], 'LSQ13dpa': [50, 180, 130],
              '2013hj': [70, 180, 105], '2014G': [50, 130, 90], '2014cx': [70, 150, 110], '2014dw': [50, 150, 90],
              'ASASSN-14dq': [60, 150, 100], 'ASASSN-14ha': [110, 170, 135], '2016X': [60, 135, 95],
              '2016bkv': [50, 240, 130], '2016gfy': [70, 170, 110], '2017eaw': [90, 170, 120],
              '2017gmr': [50, 140, 95], '2018zd': [65, 185, 125], '2018gj': [45, 125, 85]}

list1 = ['2005af.asc', '2009N.asc', '2013ab.asc', '2013ej.asc', '2014cx.asc']
list2 = ['2004et.asc', '2009ib.asc', '2012aw.asc', '2013K.asc', '2013hj.asc']
list3 = ['2003gd.asc', '2004dj.asc', '2005cs.asc', '2009md.asc', 'ASASSN-14ha.asc']
list4 = ['2008gz.asc', '2007it.asc', '2012A.asc', 'LSQ13dpa.asc', '2017eaw.asc']
list5 = ['2003hn.asc', '2004ej.asc', '2004fx.asc', '2008in.asc', '2012ec.asc']
list6 = ['2002hx.asc', '2013by.asc', '2014G.asc', '2014dw.asc', '2016X.asc']

# plot_steepness(list1, [40, 180], 0.32, 1)
# plot_steepness(list2, [20, 280], 0.25, 2)
# plot_steepness(list3, [40, 210], 0.55, 3)
# plot_steepness(list4, [10, 210], 0.36, 4)
# plot_steepness(list5, [40, 180], 0.36, 5)
# plot_steepness(list6, [20, 180], 0.28, 6)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Steepness Parameter For SNe In Study
# ------------------------------------------------------------------------------------------------------------------- #
date_explosion = 2458127.8

for file_name in ['2018gj.dat']:
    name = file_name.split('/')[-1].split('.')[0]
    data_df = read_data(file_name)
    xaxis = np.linspace(data_df['Phase'].min(), data_df['Phase'].max(), 1000)
    elmmag, gradelm, steepnesselm, inflectionelm = fit_data(name, data_df, xaxis, epoch_tran=dict_epoch[name][2])

    fig = plt.figure(figsize=(7, 12))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212, sharex=ax1)

    ax1.plot(data_df['Phase'], data_df['V'], marker='*', c='k', ls='', markerfacecolor='dodgerblue', ms=10, label=name)
    ax1.plot(xaxis, elmmag, color='orangered', linestyle='--', label='_nolegend_')
    ax2.plot(xaxis, gradelm, color='orangered')

    set_plotparams(ax1)
    set_plotparams(ax2, yticks=(0.1, 0.02))

    ax1.invert_yaxis()
    ax1.set_xlim(35, 120)
    ax2.set_ylim(0.00, 0.17)
    ax1.set_ylim(18.2, 14.9)
    ax1.legend(markerscale=0, markerfirst=False, fontsize=14, frameon=False, loc=3)

    ax1.axvline(inflectionelm, linestyle='--', color='k', linewidth=1)
    ax2.axvline(inflectionelm, linestyle='--', color='k', linewidth=1)
    ax2.axhline(steepnesselm, linestyle='--', color='dodgerblue', linewidth=1)
    ax2.text(inflectionelm - 3, 0.08, r'$\rm t_i$={0:5.1f} d'.format(inflectionelm), rotation=90, fontsize=14)
    ax2.text(inflectionelm + 10, steepnesselm + 0.003, 'Steepness={0:.3f}'.format(steepnesselm), fontsize=14)

#     ax2.text(inflectionelm * 1.1, steepnesselm * 1.1, r'$t_i$={0:5.1f} d'.format(inflectionelm), fontsize=14)
#     ax2.text(inflectionelm * 0.8, steepnesselm * 1.1, r'$S_e$={0:.3f}'.format(steepnesselm), fontsize=14)
#     valopt, valcov = curve_fit(valfunc, data_df['Phase'], data_df['V'], p0=[2, 100, 4, 0.01, 10])
#     valmag = valfunc(xaxis, *valopt)
#     fit_chisqval = chisquare(data_df['Phase'], elmfunc(data_df['Phase'], *elmopt))
#     gradval = np.gradient(valmag, dx)
#     steepnessval = gradval.max()
#     inflectionval = xaxis[np.where(gradval == gradval.max())][0]

#     ax1.plot(xaxis, valmag, color='c', linestyle='--', label='ValentiFunc')
#     ax2.plot(xaxis, gradval, color='c')
#     ax1.axvline(inflectionval, linestyle='--',  color='c', linewidth=1)
#     ax2.axvline(inflectionval, linestyle='--',  color='c', linewidth=1)
#     ax2.text(inflectionelm * 1.1, steepnesselm * 1.3, r'$t_i$={0:5.1f} d'.format(inflectionval), fontsize=14)
#     ax2.text(inflectionelm * 0.8, steepnesselm * 1.3, r'$S_v$={0:.3f}'.format(steepnessval), fontsize=14)
#     print ('Chi-Square [ValentiFunc]: {0}'.format(fit_chisqval[0]))
#     print ('Plateau Length : {0}'.format(valopt[1]))

    ax2.yaxis.set_major_locator(MultipleLocator(0.05))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.005))
    ax2.xaxis.set_major_locator(MultipleLocator(20))
    ax2.xaxis.set_minor_locator(MultipleLocator(2))
    ax1.set_ylabel(r'Apparent Magnitude, $\rm m_V$ [mag]', fontsize=14)
    ax2.set_ylabel(r'$|\rm dm_V/dt|\ [mag\ d^{-1}]$', fontsize=14)
    ax2.set_xlabel('Time Since Explosion [Days]', fontsize=14)

    fig.subplots_adjust(hspace=0.01)
    fig.savefig('PLOT_2018gjSteepness.pdf', format='pdf', dpi=2000, bbox_inches='tight')
    plt.show()
    plt.close(fig)

    lognimass = -3.5024 * steepnesselm - 1.0167
    lognierr = abs(((0.0960 * steepnesselm) ** 2 + (0.0034 ** 2)) ** 0.5)
    nimass = 10 ** lognimass
    nierr = 10 ** (lognimass + lognierr) - nimass

    print ('Nickel Mass = {0:.3f}+/-{1:.3f}'.format(nimass, nierr))
# ------------------------------------------------------------------------------------------------------------------- #
