#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------CALCULATION OF AIRMASS----------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import math
import ephem
import easygui
import datetime
import pandas as pd
import easygui as eg
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Files to Be Read
# ------------------------------------------------------------------------------------------------------------------- #
file_telescopes = 'TelescopeList.dat'
file_headers = 'HeaderInfo.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_keyword = 'OBJECT'
DATE_keyword = 'DATE-OBS'
RA_keyword = 'RA'
DEC_keyword = 'DEC'
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
# Functions For Formatting Header and Calculating AIRMASS, JD etc.
# ------------------------------------------------------------------------------------------------------------------- #

def init_telescope(telescopename, telescope_df):
    """
    Defines a Telescope object in Ephem based on the data of telescopes in 'telescope_df'.
    Args:
        telescopename : Name of the Telescope from which the data was observed
        telescope_df  : Pandas DataFrame containing data of telescopes
    Returns:
        telescope     : Ephem.Observer object containing Telescope site details
    """
    _, OBS_LONG, OBS_LAT, OBS_ALT, _, _, _ = telescope_df.loc[telescopename].values

    telescope = ephem.Observer()
    telescope.lon = OBS_LONG
    telescope.lat = OBS_LAT
    telescope.elevation = OBS_ALT
    telescope.pressure = 0
    telescope.epoch = ephem.J2000

    return telescope


def identify_imagetype(filename):
    """
    Identify the Object Type of the file 'filename'.
    Args:
        filename : FITS file for which the object type has to be appended
    Returns:
        True     : When the file is an object
        False    : When the file is a bias or a flat or a lamp 
    """
    with fits.open(filename, mode='update') as hdulist:
        header = hdulist[0].header

        if header[OBJECT_keyword] in ['BIAS', 'FLAT', 'BIASSPEC', 'FeAr', 'FeNe']:
            return False
        else:
            return True


def update_header(filepath, header_df):
    """
    Update the header details of the file indicated by path 'filepath' with the updated header
    information in Pandas DataFrame 'header_df'.
    Args:
        filepath   : Path to the FITS file whose header details have to be updated
        header_df  : Pandas DataFrame containing updated header of files
    Returns:
        None
    """
    filename = filepath.split('/')[-1]
    if filename in header_df.index:
        with fits.open(filepath, mode='update') as hdulist:
            header = hdulist[0].header
            columns = header_df.columns
            for column in columns:
                val = header_df.loc[header_df.index == filename, column].values[0]
                if column.upper() in header.keys():
                    header.set(column.upper(), str(val).upper(), '')
    else:
        display_text("ERROR: File '{0}' is not logged in '{1}'".format(filepath, file_headers))
        display_text("ERROR: Run 'PrintHeaderInfo.py' before running this script")
        pass


def append_details(filename, telescopename, telescope_df, extn=0):
    """
    Calculates Header Info (JD, AIRMASS etc.) for the FITS file 'filename'.
    Args:
        filename      : FITS file for which header details have to be appended
        telescopename : Name of the Telescope from which the data was observed
        telescope_df  : Pandas DataFrame containing data of telescopes
        extn          : Extension of the FITS file which stores the complete header
    Returns:
        None
    """
    OBS_NAME, OBS_LONG, OBS_LAT, OBS_ALT, OBS_TIMEZONE, _, _ = telescope_df.loc[telescopename].values

    telescope = ephem.Observer()
    telescope.lon = OBS_LONG
    telescope.lat = OBS_LAT
    telescope.elevation = OBS_ALT
    telescope.pressure = 0
    telescope.epoch = ephem.J2000

    with fits.open(filename, mode='update') as hdulist:
        file_header = hdulist[extn].header

        date_avg = file_header[str(DATE_keyword)]
        date_obs, time_utc = date_avg.split('T')
        datetime_utc = str(date_obs) + ' ' + str(time_utc)
        julian_day = round(ephem.julian_date(datetime_utc), 3)

        telescope.date = datetime_utc
        time_sidereal = telescope.sidereal_time()                                                                    

        source = ephem.FixedBody()
        source._ra = file_header[RA_keyword]
        source._dec = file_header[DEC_keyword]
        source._epoch = ephem.J2000
        source.compute(telescope)

        object_alt = Angle(str(source.alt) + ' degrees').degree
        airmass = round(1 / math.cos(math.radians(90 - object_alt)), 3)

        # List of Header Keywords to be Appended
        dict_header = {'TELESCOP': OBS_NAME, 'GEOLONG': OBS_LONG, 'GEOLAT': OBS_LAT, 'GEOELEV': OBS_ALT,
                        'TIMEZONE': OBS_TIMEZONE, DATE_keyword: date_obs, 'UT': time_utc, 'JD': julian_day,
                        'ST': time_sidereal, 'ALTITUDE': source.alt, 'AZIMUTH': source.az, 'AIRMASS': airmass}

        for keyword, value in dict_header.items():
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all=True)
            file_header.append(card=(keyword, str(value)))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Main Function
# ------------------------------------------------------------------------------------------------------------------- #

def main():
    """
    Step 1: GUI Code for User Input
    Step 2: Group FITS files and print the Header Info
    """
    # GUI Code for User Input
    DIR_FILES = eg.enterbox('Enter the directory in which headers of files have to be updated:',
                            title='Enter the Directory Path', default=[os.path.join(os.getcwd(), 'preprocessed')])
    
    telescopename = eg.enterbox('Enter the Name of the Telescope from which the data was observed:',
                                title='Enter the Name of the Telescope', default=['HCT'])

    instrument = eg.enterbox('Enter the Instrument from which the data was observed:',
                                 title='Enter the Short Name of the Instrument', default=['HFOSC2'])

    input_file = eg.enterbox('Enter the name of the output file containing the header info:',
                             title='Enter the Name of the Output File', default=['HeaderInfo.dat'])

    # Group FITS Files whose header are to be Updated + Read Input File
    if instrument == 'HFOSC2':
        ctext = 'fix*.fits'
    else:
        ctext = 'raw*.fits'
    list_files = group_similar_files('', common_text=os.path.join(DIR_FILES, ctext))
    header_df = pd.read_csv(input_file, sep='\s*[|]\s*', comment='#', dtype='string', engine='python')
    header_df = header_df.set_index('FILENAME')

    # Initialises the Telescope Object using Ephem after Extracting Details from 'file_telescopes'
    telescope_df = pd.read_csv(file_telescopes, sep='\s+', comment='#').set_index('ShortName')

    # Calculates AIRMASS etc. & Appends Respective Details In The Header
    for filename in list_files:
        update_header(filename, header_df)
        if identify_imagetype(filename):
            append_details(filename, telescopename, telescope_df)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Execute the Standalone Code
# ------------------------------------------------------------------------------------------------------------------- #
if __name__ == '__main__':
    main()
# ------------------------------------------------------------------------------------------------------------------- #