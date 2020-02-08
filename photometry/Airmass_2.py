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
date_keyword = 'DATE-OBS'
object_keyword = 'OBJECT'
airmass_keyword = 'AIRMASS'
time_start_keyword = 'TM_START'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
RA_object = '21:57:59.9'
DEC_object = '+24:16:08.1'
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
    global RA_object, DEC_object

    file_header = fits.getheader(filename=str(file_name))
    date_obs = file_header[str(date_keyword)]

    file_header.set('DATE_OBS', str(date_obs))
    file_header.set(str(RA_keyword), RA_object)
    file_header.set(str(DEC_keyword), DEC_object)


def calculate_airmass(file_name):
    """
    Calculates AIRMASS for the FITS file and appends respective details in the header of the file 'file_name'
    Args:
        file_name : FITS file whose header has to be edited
    Returns:
        None
    """
    hdulist = fits.open(file_name, mode='update')
    file_header = hdulist[0].header
    date_obs = file_header[str(date_keyword)]
    time_start = file_header[str(time_start_keyword)]

    if str(RA_keyword) in file_header.keys():
        object_ra = file_header[str(RA_keyword)]
    else:
        object_ra = RA_object

    if str(DEC_keyword) in file_header.keys():
        object_dec = file_header[str(DEC_keyword)]
    else:
        object_dec = DEC_object

    time_utc = str(datetime.timedelta(seconds=int(time_start)))
    datetime_utc = str(date_obs) + ' ' + str(time_utc)
    julian_day = ephem.julian_date(datetime_utc)

    telescope = ephem.Observer()
    telescope.lon = OBS_LONG
    telescope.lat = OBS_LAT
    telescope.elevation = OBS_ALT
    telescope.pressure = 0
    telescope.epoch = ephem.J2000
    telescope.date = datetime_utc

    obj_pos = ephem.FixedBody()
    obj_pos._ra = object_ra
    obj_pos._dec = object_dec
    obj_pos._epoch = ephem.J2000
    obj_pos.compute(telescope)

    time_sidereal = telescope.sidereal_time()
    object_alt = Angle(str(obj_pos.alt) + ' degrees').degree
    airmass = 1 / math.cos(math.radians(90 - object_alt))

    list_keywords = ['OBSERVAT', 'OBS_LAT', 'OBS_LONG', 'OBS_ALT', 'TIMEZONE', 'DATE_OBS', 'UT', 'JD', 'ST', 'RA',
                     'DEC', 'ALT', 'AZ', 'AIRMASS']

    dict_header = {'OBSERVAT': str(OBS_NAME), 'OBS_LAT': str(OBS_LAT), 'OBS_LONG': str(OBS_LONG),
                   'OBS_ALT': str(OBS_ALT), 'TIMEZONE': str(OBS_TIMEZONE), 'DATE_OBS': str(date_obs),
                   'UT': str(time_utc), 'JD': str(julian_day), 'ST': str(time_sidereal), 'RA': str(object_ra),
                   'DEC': str(object_dec), 'ALT': str(obj_pos.alt), 'AZ': str(obj_pos.az), 'AIRMASS': str(airmass)}

    for keyword in list_keywords:
        if keyword in file_header.keys():
            file_header.remove(str(keyword), remove_all=True)
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
