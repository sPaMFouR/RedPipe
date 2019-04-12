#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxx-------------------PHOTOMETRY OF GALAXIES-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import glob
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from photutils import CircularAperture, SkyEllipticalAperture, SkyCircularAperture, aperture_photometry
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# Galaxy Details
# ------------------------------------------------------------------------------------------------------------------- #
name_glx = 'NGC 2276'
RA_glx = '07:27:14.36'
DEC_glx = '+85:45:16.4'
a_glx = 169.10
b_glx = 161.49
pa_glx = 23.0
dist_glx = 31.3       # In Mpc
AB_glx = 0.62         # Self Extinction
EBV_self = AB_glx / 4.06
EBV_los = 0.0859
light_speed = 3e10

dict_cpstoflux = {'FUV': 1.4e-15, 'NUV': 2.06e-16}
dict_zpAB = {'FUV': 18.82, 'NUV': 20.08}
dict_zpfluxVegajy = {'IRACI1': 277.22, 'IRACI2': 179.04, 'IRACI3': 113.85, 'IRACI4': 62.00}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# Location Of Directories
# ------------------------------------------------------------------------------------------------------------------- #
EXP_keyword = 'EXPTIME'
ZP_ps1comp = 'ZPT_'
ZP_ps1 = 'HIERARCH FPA.ZP'
ZP_2mass = 'MAGZP'
FLUXCONV_spitz = 'FLUXCONV'
xsize_spitz = 'PXSCAL1'
ysize_spitz = 'PXSCAL2'

DIR_CODE = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/"
DIR_GAL = "/home/avinash/Supernovae_Data/2016gfy/NGC2276/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# DataFrame Containing Information On FILTERS and Extinction Coefficients (For Rv = 3.1, Fitzpatrick(1999))
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
# Function For File Handling & Displaying Text
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


def display_text(text_to_display):
    """
    Displays text mentioned in the string 'text_to_display'
    Args:
        text_to_display : Text to be displayed
    Returns:
        None
    """
    print("\n" + "# " + "-" * (12 + len(text_to_display)) + " #")
    print("# " + "-" * 5 + " " + str(text_to_display) + " " + "-" * 5 + " #")
    print("# " + "-" * (12 + len(text_to_display)) + " #" + "\n")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions For Performing Photometry
# ------------------------------------------------------------------------------------------------------------------- #

def read_wcs(file_name):
    hdulist = fits.open(file_name)[0]
    wcs = WCS(hdulist.header)
    return wcs


def get_spitzer_zp(header):
    list_zp = []
    for val in header.keys():
        if ZP_ps1comp in val:
            list_zp.append(float(header[val]))
    return np.mean(list_zp)


def photgalexAB(file_name, aperture, fuv=True):
    data = fits.getdata(file_name, 0)
    counts = float(aperture_photometry(data, aperture)['aperture_sum'])

    if fuv:
        mag = -2.5 * np.log10(counts) + dict_zpAB['FUV']
        extcmag = mag - filter_df.loc['FUV', 'RLambda'] * (EBV_los + EBV_self)
        sfr = 10 ** (2.78 - 0.4 * extcmag + 2 * np.log10(dist_glx))
        print "FUV : CPS = {0:0.4f}, ABMag = {1:0.2f}, ExtCorMag = {2:0.2f}".format(counts, mag, extcmag)
        print "SFR : {0:0.2f}".format(sfr)
    else:
        mag = -2.5 * np.log10(counts) + dict_zpAB['NUV']
        extcmag = mag - filter_df.loc['FUV', 'RLambda'] * EBV_los
        print "NUV : CPS = {0:0.4f}, ABMag = {1:0.2f}, ExtCorMag = {2:0.2f}".format(counts, mag, extcmag)


def perform_photometry(file_path, aperture, dict_mag, skycoord=False):
    hdu = fits.open(file_name)[0]
    header = hdu.header
    band = file_path.split('.')[0].split('_')[-1]

#     data = np.array(hdu.data)
#     data[data == 'NULL'] = 0
#     data = data.astype('float64')
    if not skycoord:
        counts = float(aperture_photometry(hdu.data, aperture)['aperture_sum'])
    else:
        counts = float(aperture_photometry(hdu.data, aperture, wcs=read_wcs(file_path))['aperture_sum'])

    if 'Spitzer' in file_name:
        fluxjy = counts * 1e6 * ((np.pi / (180 * 60 * 60)) ** 2) * abs(header[xsize_spitz] * header[xsize_spitz])
        mag_vega = -2.5 * np.log10(fluxjy / dict_zpfluxVegajy[band])
        flux_vega = 10 ** (-0.4 * (mag_vega + 21.1 + filter_df.loc[band, 'ZeroPoint']))
    elif '2MASS' in file_name:
        mag_vega = -2.5 * np.log10(counts) + header[ZP_2mass]
        flux_vega = 10 ** (-0.4 * (mag_vega + 21.1 + filter_df.loc[band, 'ZeroPoint']))
    elif 'Pan-Starrs' in file_name:
        exptime = header[EXP_keyword]
        cps = counts / exptime
        mag_vega = -2.5 * np.log10(cps) + header[ZP_ps1]
        flux_vega = 10 ** (-0.4 * (mag_vega + 21.1 + filter_df.loc[band, 'ZeroPoint']))
    elif 'GalEX' in file_name:
        flux_vega = dict_cpstoflux[band] * counts
        mag_vega = -2.5 * np.log10(flux_vega) - filter_df.loc[band, 'ZeroPoint'] - 21.1
    else:
        print ("ERROR : Unidentified DataSet")

    flux_ab = 1e8 * flux_vega * ((filter_df.loc[band, 'CentreWave'] * 1e-8) ** 2) / light_speed
    mag_ab = -2.5 * np.log10(flux_ab) - filter_df.loc[band, 'ZeroPoint'] - 48.6

    cormag_vega = mag_vega - filter_df.loc[band, 'RLambda'] * EBV_los
    corflux_vega = 10 ** (-0.4 * (cormag_vega + 21.1 + filter_df.loc[band, 'ZeroPoint']))
    corflux_ab = 1e8 * corflux_vega * ((filter_df.loc[band, 'CentreWave'] * 1e-8) ** 2) / light_speed
    cormag_ab = -2.5 * np.log10(corflux_ab) - filter_df.loc[band, 'ZeroPoint'] - 48.6
    corflux_jy = 1e23 * corflux_ab

    dict_mag[band] = {}
    dict_mag[band] = {'CentreWave': filter_df.loc[band, 'CentreWave'], 'AppMagVega': mag_vega, 'AppFluxVega': flux_vega,
                      'AppMagAB': mag_ab, 'AppFluxAB': flux_ab, 'CorMagVega': cormag_vega, 'CorFluxVega': corflux_vega,
                      'CorMagAB': cormag_ab, 'CorFluxAB': corflux_ab, 'CorFluxJy': corflux_jy}

    return dict_mag

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Group Files From Different Telescopes
# ------------------------------------------------------------------------------------------------------------------- #
list_galex = group_similar_files('', common_text=DIR_GAL + 'GalEX/BS_*.fits')
list_ps1 = group_similar_files('', common_text=DIR_GAL + 'Pan-Starrs/BS_*.fits')
list_spitz = group_similar_files('', common_text=DIR_GAL + 'Spitzer/BS_*.fits')
list_2mass = group_similar_files('', common_text=DIR_GAL + '2MASS/BS_*.fits')

list_files = list_galex + list_ps1 + list_spitz + list_2mass
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Apertures For Different Datasets In Pixel Coordinates
# ------------------------------------------------------------------------------------------------------------------- #
aperture_final = SkyCircularAperture(SkyCoord('7:27:33.040', '+85:45:28.078', unit=(u.hourangle, u.deg), frame='fk5'),
                                     r=79.914 * u.arcsec)

aperture_ps1 = CircularAperture((511.99622, 635.88791), r=319.65599)
aperture_galex = CircularAperture((1230.4901, 2327.5328), r=53.276)
aperture_spitz = CircularAperture((1096.5952, 544.47686), r=133.18973)
aperture_2mass = CircularAperture((106.53127, 74.079436), r=79.913998)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate SFR Of The Host Galaxy Using GalEX FUV Magnitude (In AB System)
# ------------------------------------------------------------------------------------------------------------------- #
display_text("Calculating Star Formation Rate in {0}".format(name_glx))

for file_name in list_galex:
    if 'FUV' in file_name:
        photgalexAB(file_name, aperture_galex)
    else:
        photgalexAB(file_name, aperture_galex, fuv=False)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Photometry Using SkyCoord Object With Photutils
# Calculates Photometric Magnitudes And Fluxes In Different Bands
# ------------------------------------------------------------------------------------------------------------------- #
dict_mag = {}
for file_name in list_files:
    dict_mag = perform_photometry(file_name, aperture_final, dict_mag, skycoord=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs Photometry Using PixelCoord Object With Photutils
# Calculates Photometric Magnitudes And Fluxes In Different Bands
# ------------------------------------------------------------------------------------------------------------------- #
# dict_mag = {}
# for file_name in list_galex:
#     dict_mag = perform_photometry(file_name, aperture_galex, dict_mag)

# for file_name in list_ps1:
#     dict_mag = perform_photometry(file_name, aperture_ps1, dict_mag)

# for file_name in list_2mass:
#     dict_mag = perform_photometry(file_name, aperture_2mass, dict_mag)

# for file_name in list_spitz:
#     dict_mag = perform_photometry(file_name, aperture_spitz, dict_mag)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Outputs The Photometric Information Onto A Data File
# ------------------------------------------------------------------------------------------------------------------- #
display_text("Magnitudes & Fluxes of {0}".format(name_glx))

data_df = pd.DataFrame(dict_mag).T
data_df.index.name = 'FILTER'
data_df = data_df.reset_index().sort_values(by='CentreWave')


def fmt(col):
    if 'Flux' in col:
        return "{0:.3e}"
    elif 'Mag' in col:
        return "{0:.2f}"
    else:
        return None


data_tab = Table.from_pandas(data_df)
data_tab.write('NGC2276_Mag.dat', format='ascii.fixed_width', delimiter=' ',
               formats={x: fmt(x) for x in data_df.columns}, overwrite=True)
print data_tab
# ------------------------------------------------------------------------------------------------------------------- #
