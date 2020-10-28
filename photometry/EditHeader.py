#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxx---------MODIFIES HEADER IN THE HFOSC2 FITS FILES AND APPENDS AIRMASS DETAILS------------xxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

__version__ = "---"
__author__ = "Avinash Singh"

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import re
import glob
import math
import ephem
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_LONG = '79:41:06'
OBS_LAT = '29:21:42'
OBS_ALT = 2450
OBS_TIMEZONE = +5.5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_keyword = 'OBJECT'
EXPTIME_keyword = 'EXPTIME'
UT_keyword = 'UT'
RA_keyword = 'OBJRA'
DEC_keyword = 'OBJDEC'
DATE_keyword = 'DATE-OBS'
DATEAVG_keyword = 'DATE-AVG'
BIASSEC_keyword = 'BIASSEC'
CCDSEC_keyword = 'CCDSEC'
TRIMSEC_keyword = 'TRIMSEC'

modRA_keyword = 'RA'
modDEC_keyword = 'DEC'

OBJ_RA = '12:26:12.05'
OBJ_DEC = '58:18:51.10'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def group_similar_files(text_list, common_text, exceptions=''):
    """
    Groups similar files based on the string "common_text". Writes the similar files
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
                test = re.search(str(text), file_name)
                print(file_name)
                if test:
                    try:
                        list_files.remove(file_name)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(str(text_list), "w") as f:
            for index in range(0, len(list_files)):
                f.write(str(list_files[index]) + "\n")

    return list_files

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates AIRMASS And Appends Respective Details In The Header
# ------------------------------------------------------------------------------------------------------------------- #
# date_mod = raw_input("Enter the prefix of recently observed files (Example: h2abj01 - 2018/10/01): ")
list_files = group_similar_files('', common_text='*.fits')

for file_name in list_files:
    hdulist = fits.open(file_name, mode='update')
    file_header = hdulist[0].header
    print file_name
    if str(RA_keyword) in file_header.keys():
        object_ra = file_header[str(RA_keyword)]
    else:
        object_ra = OBJ_RA

    if str(DEC_keyword) in file_header.keys():
        object_dec = file_header[str(DEC_keyword)]
    else:
        object_dec = OBJ_DEC

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
    list_keywords = ['LAT', 'LONG', 'ALT', 'TIMEZONE', modRA_keyword, modDEC_keyword, UT_keyword,
                     DATE_keyword, 'JD', 'ST', 'ELE', 'AZ', 'AIRMASS']
    dict_header = {'LAT': str(OBS_LAT), 'LONG': str(OBS_LONG), 'ALT': str(OBS_ALT), 'TIMEZONE': str(OBS_TIMEZONE),
                   modRA_keyword: OBJ_RA, modDEC_keyword: OBJ_DEC, DATE_keyword: str(date_obs),
                   UT_keyword: str(time_utc), 'JD': str(julian_day), 'ST': str(time_sidereal),
                   'ELE': str(object_pos.alt), 'AZ': str(object_pos.az), 'AIRMASS': str(airmass)}

    if object_ra == '' or object_dec == '':
        list_keywords = list_keywords[:-3]

    for keyword in [TRIMSEC_keyword, BIASSEC_keyword, CCDSEC_keyword]:
        if keyword in file_header.keys():
            file_header.remove(keyword, remove_all=True)

    for keyword in list_keywords:
        if keyword in file_header.keys():
            file_header.remove(keyword, remove_all=True)
        file_header.append(card=(keyword, dict_header[keyword]))

    hdulist.close()

# ------------------------------------------------------------------------------------------------------------------- #
