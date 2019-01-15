#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------CALCULATION OF AIRMASS----------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import glob
import re
import math
import ephem
import datetime
from pyraf import iraf
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory Site Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = "Indian Astronomical Observatory, Hanle"
# OBS_LONG = 78.9642
# OBS_LAT = 32.7794
OBS_LONG = '78:57:51'
OBS_LAT = '32:46:46'
OBS_ALT = 4486
OBS_TIMEZONE = +5.5
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Telescope CCD Specifications
# ------------------------------------------------------------------------------------------------------------------- #
read_noise = 4.87
ccd_gain = 1.22
data_max = 55000
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
RA_keyword = 'RA'
DEC_keyword = 'DEC'
date_keyword = 'DATE-OBS'
grism_keyword = 'GRISM'
filter_keyword = 'IFILTER'
object_keyword = 'OBJECT'
airmass_keyword = 'AIRMASS'
exptime_keyword = 'EXPTIME'
time_start_keyword = 'TM_START'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Object Details
# ------------------------------------------------------------------------------------------------------------------- #
name_object = 'ASASSN14dq'
RA_object = '21:57:59.9'
DEC_object = '+24:16:08.1'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Location Of AIRMASS File
# ------------------------------------------------------------------------------------------------------------------- #
FILE_AIRMASS = '/home/avinash/PyCharmProjects/Reduction_Pipeline/Airmass.asc'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Standard Stars Details
# ------------------------------------------------------------------------------------------------------------------- #
mag_Feige110 = 11.83
RA_Feige110 = '23:19:58.4'
DEC_Feige110 = '-05:09:56.2'
mag_Feige34 = 11.18
RA_Feige34 = '10:39:36.7'
DEC_Feige34 = '+43:06:09.3'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.astutil(_doprint=0)
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
                print file_name
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


def text_list_to_python_list(text_list):
    """
    Returns data in the file 'text_list' as a python_list.
    Args:
        text_list   : Input file containing filenames
    Returns:
        python_list : List of all the elements in the file 'text_list'
    Raises:
        Error : File 'text_list 'Not Found
    """
    if os.path.isfile(text_list):
        with open(text_list, "r+") as f:
            python_list = f.read().split()
            return python_list
    else:
        print "Error : File " + str(text_list) + " Not Found"

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def asthedit(text_list_edit, file_airmass):
    """
    Calculates AIRMASS for the files in the list 'text_list_edit' using the commands file 'file_airmass'.
    Also, appends AIRMASS, RA, DEC, JD etc. keywords to the header of the file 'file_name'.
    Args:
        text_list_edit : Text list of FITS files whose headers have to be appended with AIRMASS
        file_airmass   : AIRMASS file which has to be edited for date of observation
    Returns:
        None
    """
    list_edit = text_list_to_python_list(str(text_list_edit))

    task = iraf.noao.astutil.asthedit
    task.unlearn()

    task.verbose = "no"                                 # Verbose Output?
    task.update = "yes"                                 # Update Image Header?
    task.oldstyle = "no"                                # Use Old Style Format?

    for file_name in list_edit:
        edit_header(file_name)
        edit_airmass_file(file_name, FILE_AIRMASS)
        task(images=str(file_name), commands=str(file_airmass))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Airmass Correction
# ------------------------------------------------------------------------------------------------------------------- #

def edit_airmass_file(file_name, file_airmass):
    """
    Edits the 'DATE_OBS' variable in file 'file_airmass' based on the date of observation of the file 'file_name'.
    Args:
        file_name    : FITS file from which date has to be extracted
        file_airmass : AIRMASS file which has to be edited for date of observation
    Returns:
        None
    """
    with open(str(file_airmass), 'r') as f:
        f.readline()
        rest_data = f.read()

    # DATE_OBS = '2014-09-16'
    file_header = fits.getheader(file_name)
    date_obs = file_header[str(date_keyword)]

    new_date = "DATE_OBS = " + "'" + str(date_obs) + "'" + "\n"
    with open(str(file_airmass), 'w') as f:
        f.write(new_date)
        f.write(rest_data)


def edit_header(file_name):
    """
    Edits the header of the file 'file_name' based on the object_keyword.
    Args:
        file_name : FITS file whose header has to be edited
    Returns:
        None
    """
    (file_data, file_header) = fits.getdata(filename=str(file_name), header=True, ext=0)
    object_value = file_header[str(object_keyword)]
    date_obs = file_header[str(date_keyword)]

    object_ra = RA_object
    object_dec = DEC_object

    if object_value == 'Feige110':
        object_ra = RA_Feige110
        object_dec = DEC_Feige110
    elif object_value == 'Feige34':
        object_ra = RA_Feige34
        object_dec = DEC_Feige34

    file_header.set('DATE_OBS', str(date_obs))
    file_header.set(str(RA_keyword), object_ra)
    file_header.set(str(DEC_keyword), object_dec)

    fits.writeto(filename=str(file_name), data=file_data, header=file_header, clobber=True)


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
    object_ra = file_header[str(RA_keyword)]
    object_dec = file_header[str(DEC_keyword)]

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

    object_pos = ephem.FixedBody()
    object_pos._ra = object_ra
    object_pos._dec = object_dec
    object_pos._epoch = ephem.J2000
    object_pos.compute(telescope)

    time_sidereal = telescope.sidereal_time()
    object_alt = Angle(str(object_pos.alt) + ' degrees').degree
    airmass = 1 / math.cos(math.radians(90 - object_alt))

    file_header.append('OBSERVAT', str(OBS_NAME))
    file_header.append('LAT', str(OBS_LAT))
    file_header.append('LONG', str(OBS_LONG))
    file_header.append('ALT', str(OBS_ALT))
    file_header.append('TIMEZONE', str(OBS_TIMEZONE))
    file_header.append('DATE_OBS', str(date_obs))
    file_header.append('UT', str(time_utc))
    file_header.append('JD', str(julian_day))
    file_header.append('ST', str(time_sidereal))
    file_header.append('RA', str(object_ra))
    file_header.append('DEC', str(object_dec))
    file_header.append('ALT', str(object_pos.alt))
    file_header.append('AZ', str(object_pos.az))
    file_header.append('AIRMASS', str(airmass))

    hdulist.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates AIRMASS And Appends Respective Details In The Header
# ------------------------------------------------------------------------------------------------------------------- #
list_files = group_similar_files("list_files", common_text="*.fits")

for file_name in list_files:
    calculate_airmass(file_name=file_name)

# asthedit("list_files", file_airmass=FILE_AIRMASS)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# AIRMASS Test Not Verified
# ------------------------------------------------------------------------------------------------------------------- #
# from PyAstronomy import pyasl
# from astropy.coordinates import SkyCoord, EarthLocation, AltAz
# from astropy.utils import iers
# iers.conf.auto_download = False
#
# lat_deg = Angle(OBS_LAT + ' degrees').degree
# long_deg = Angle(OBS_LONG + ' degrees').degree
# ra_deg, dec_deg = pyasl.coordsSexaToDeg(str(object_ra) + ' ' + str(object_dec))
#
# object_radec = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
# telescope_loc = EarthLocation(lat=lat_deg*u.deg, lon=-long_deg*u.deg, height=OBS_ALT*u.m)
# utcoffset = OBS_TIMEZONE*u.hour
# time_utc = Time(str(date_obs) + ' ' + str(time_obs))
#
# object_altaz = object_radec.transform_to(AltAz(obstime=time_utc, location=telescope_loc))
# airmass = 1 / math.cos(math.radians(90 - object_altaz[0]))
# print object_radec
# print object_altaz
# print airmass
#
# from sidereal import JulianDate, SiderealTime
#
# JulianDate.fromDatetime(date_obs)
# SiderealTime
# ------------------------------------------------------------------------------------------------------------------- #
