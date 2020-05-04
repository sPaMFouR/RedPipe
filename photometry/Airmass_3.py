#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------CALCULATION OF AIRMASS----------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import glob
import math
import ephem
import easygui
import datetime
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = "Indian Astronomical Observatory, Hanle"
OBS_LONG = '78:57:51'
OBS_LAT = '32:46:46'
OBS_ALT = 4486
OBS_TIMEZONE = +5.5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
RA_keyword = 'RA'
DEC_keyword = 'DEC'
UT_keyword = 'UT'
DATE_keyword = 'DATE-OBS'
DATEAVG_keyword = 'DATE-OBS'
OBJECT_keyword = 'OBJECT'
EXPTIME_keyword = 'EXPTIME'
AIRMASS_keyword = 'AIRMASS'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_RA = '21:57:59.9'
OBJECT_DEC = '+24:16:08.1'
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
# Functions For Calculating AIRMASS, JD etc.
# ------------------------------------------------------------------------------------------------------------------- #

def add_radec(file_name):
    """
    Adds the RA & DEC of the image observed to the header of the file 'file_name'.
    Args:
        file_name : FITS file to which header detail is to be added
    Returns:
        None
    """
    global OBJECT_RA, OBJECT_DEC

    file_header = fits.getheader(filename=str(file_name))
    date_obs = file_header[str(DATE_keyword)]

    file_header.set('DATE_OBS', str(date_obs))
    file_header.set(str(RA_keyword), OBJECT_RA)
    file_header.set(str(DEC_keyword), OBJECT_DEC)


def calculate_airmass(file_name):
    """
    Calculates AIRMASS for the FITS file and appends respective details in the header of the file 'file_name'
    Args:
        file_name : FITS file whose header has to be edited
    Returns:
        None
    """
    for file_name in list_files:
        hdulist = fits.open(file_name, mode='update')
        file_header = hdulist[0].header

        if str(RA_keyword) in file_header.keys():
            object_ra = file_header[str(RA_keyword)]
        else:
            object_ra = OBJECT_RA

        if str(DEC_keyword) in file_header.keys():
            object_dec = file_header[str(DEC_keyword)]
        else:
            object_dec = OBJECT_DEC

        date_avg = file_header[str(DATEAVG_keyword)]
        date_obs, time_utc = date_avg.split('T')

        datetime_utc = str(date_obs) + ' ' + str(time_utc)
        julian_day = ephem.julian_date(datetime_utc)

        telescope = ephem.Observer()
        telescope.lon = OBS_LONG
        telescope.lat = OBS_LAT
        telescope.elevation = OBS_ALT
        telescope.pressure = 0
        telescope.epoch = ephem.J2000
        telescope.date = datetime_utc
        time_sidereal = telescope.sidereal_time()

        object_pos = ephem.FixedBody()
        object_pos._ra = object_ra
        object_pos._dec = object_dec
        object_pos._epoch = ephem.J2000
        object_pos.compute(telescope)

        object_alt = Angle(str(object_pos.alt) + ' degrees').degree
        airmass = 1 / math.cos(math.radians(90 - object_alt))
        list_keywords = ['LAT', 'LONG', 'ALT', 'TIMEZONE', RA_keyword, DEC_keyword, UT_keyword,
                         DATE_keyword, 'JD', 'ST', 'ELE', 'AZ', 'AIRMASS']
        dict_header = {'LAT': str(OBS_LAT), 'LONG': str(OBS_LONG), 'ALT': str(OBS_ALT), 'TIMEZONE': str(OBS_TIMEZONE),
                       RA_keyword: str(object_ra), DEC_keyword: str(object_dec), DATE_keyword: str(date_obs),
                       UT_keyword: str(time_utc), 'JD': str(julian_day), 'ST': str(time_sidereal),
                       'ELE': str(object_pos.alt), 'AZ': str(object_pos.az), 'AIRMASS': str(airmass)}

        for keyword in list_keywords:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all=True)
            file_header.append(card=(keyword, dict_header[keyword]))

        hdulist.flush()
        hdulist.close()


# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
ctext = easygui.enterbox(msg='Enter The Common Text Of Files For Which Airmass Is To Be Calculated?',
                         title='Airmass Calculation', default='*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates AIRMASS etc. & Appends Respective Details In The Header
# ------------------------------------------------------------------------------------------------------------------- #
for file_name in group_similar_files("", common_text=str(ctext)):
    calculate_airmass(file_name=file_name)
# ------------------------------------------------------------------------------------------------------------------- #
