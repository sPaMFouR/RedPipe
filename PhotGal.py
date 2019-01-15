#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxx-------------------PHOTOMETRY OF GALAXIES-----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import glob
import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils import SkyCircularAperture, EllipticalAperture, SkyEllipticalAperture, aperture_photometry
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variable
# ------------------------------------------------------------------------------------------------------------------- #
RA_glx = '07:27:14.36'
DEC_glx = '+85:45:16.4'
a_glx = 169.10
b_glx = 161.49
pa_glx = 23.0
dist_glx = 31.3
name_glx = 'NGC 2276'

A_B = 0.62
EBV_mag = A_B / 4.06
R_FUV = 7.45
R_NUV = 8.36
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
# Function For Performing Photometry
# ------------------------------------------------------------------------------------------------------------------- #

def photgalex(file_name, aperture, fuv=True):
    data = fits.getdata(file_name, 0)
    counts = float(aperture_photometry(data, aperture)['aperture_sum'])

    if fuv:
        mag = -2.5 * np.log10(counts) + 18.82
        extcmag = mag - R_FUV * EBV_mag
        sfr = 10 ** (2.78 - 0.4 * extcmag + 2 * np.log10(dist_glx))
        print "FUV : Counts = {0:0.4f}, Mag = {1:0.2f}, CorMag = {2:0.2f}".format(counts, mag, extcmag)
        print "SFR : {0:0.2f}".format(sfr)
    else:
        mag = -2.5 * np.log10(counts) + 20.08
        extcmag = mag - R_NUV * EBV_mag
        print "NUV : Counts = {0:0.4f}, Mag = {1:0.2f}, CorMag = {2:0.2f}".format(counts, mag, extcmag)


def phot_skycoord(file_name, aperture):
    hdulist = fits.open(file_name)[0]
    phottab = aperture_photometry(hdulist, aperture)
    return phottab


def phot_pixelcoord(file_name, aperture):
    data = fits.open(file_name)[0].data
    phottab = aperture_photometry(data, aperture)
    return phottab

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Group Files From Different Telescopes
# ------------------------------------------------------------------------------------------------------------------- #
list_galex = group_similar_files('', common_text='BS_*UV.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Apertures For Different Telescopes In Pixel Coordinates
# ------------------------------------------------------------------------------------------------------------------- #
centre = SkyCoord(RA_glx, DEC_glx, unit=(u.hourangle, u.deg), frame='fk5')
aperture_wcs = SkyEllipticalAperture(centre, a=a_glx * u.arcsec, b=b_glx * u.arcsec, theta=(pa_glx - 90) * u.deg)
aperture_galex = EllipticalAperture((2816.5582, 873.09807), a=112.73333, b=107.66, theta=23)
aperture_galex2 = EllipticalAperture((2803, 879), a=52, b=54, theta=23)

for file_name in list_galex:
    if 'FUV' in file_name:
        photgalex(file_name, aperture_galex)
        photgalex(file_name, aperture_galex2)
    else:
        photgalex(file_name, aperture_galex, fuv=False)
        photgalex(file_name, aperture_galex2, fuv=False)

# ------------------------------------------------------------------------------------------------------------------- #
