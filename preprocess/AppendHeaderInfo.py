#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx---------------UPDATES AND APPENDS HEADER INFORMATION----------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import math
import ephem
import shutil
import pandas as pd
import easygui as eg
from astropy.io import fits
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #

"""
Things to keep in mind before running this script:
--> This script is based on the output file 'HeaderInfo.dat' generated by the code 'PrintHeaderInfo.py'
--> File 'HeaderInfo.dat' should be formatted to account for errors in 'OBJECT', 'FILTER1' and 'FILTER2' keywords
--> 'OBJECT' keyword should be formatted in capital letters
    'BIAS' for bias frames
    'FLAT' for flat frames
    <OBJECT_NAME> for object frames
--> Observation in Bessel 'UBVRI' filters should be denoted in capital letters in 'FILTER1' column
--> Observation in SDSS 'ugriz' filters should be denoted in small letters in 'FILTER2' column
--> Make sure one among the two filter columns i.e. 'FILTER1' and 'FILTER2' denotes 'Clear'
"""

# ------------------------------------------------------------------------------------------------------------------- #
# Files to Be Read
# ------------------------------------------------------------------------------------------------------------------- #
file_telescopes = 'TelescopeList.dat'
file_headerinfo = 'HeaderInfo.dat'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_keyword = 'OBJECT'
DATE_keyword = 'DATE-OBS'
FILTER_keyword = 'FILTER'
GRISM_keyword = 'GRISM'
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


def file_rename(filepath, header_df):
    """
    Rename the file according to the format 'date_object-filter/grism' with the details 'filepath' with the updated header
    information in Pandas DataFrame 'header_df'.
    Args:
        filepath   : Path to the FITS file whose header details have to be updated
        header_df  : Pandas DataFrame containing updated header of files
    Returns:
        None
    """
    dirname, filename = os.path.split(filepath)

    if filename in header_df.index:
        file_data = header_df.loc[filename]
        date = file_data[DATE_keyword].split('T')[0]
        objname = file_data[OBJECT_keyword]
        bandpass = file_data[FILTER_keyword]
        grism = file_data[GRISM_keyword]

        if bandpass != 'Free':
            name_prefix = os.path.join(dirname, date + '_' + objname + '-' + bandpass)
        elif grism != 'Free':
            name_prefix = os.path.join(dirname, date + '_' + objname + '-' + bandpass)
        else:
            display_text("ERROR: Both 'Filter' and 'Grism' Keywords are 'Free'")
        
        count = len(glob.glob(name_prefix + '*.fits'))
        new_filename = name_prefix + str(count + 1) + '.fits'

        if os.path.exists(new_filename):
            os.remove(new_filename)
        shutil.copy(filepath, new_filename)
    else:
        pass

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
            file_header = hdulist[0].header
            columns = header_df.columns
            for column in columns:
                val = header_df.loc[header_df.index == filename, column].values[0]
                if column.upper() in file_header.keys():
                    file_header.set(column.upper(), str(val).upper(), '')
    else:
        display_text("ERROR: File '{0}' is not logged in '{1}'".format(filepath, file_headerinfo))
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

        date_avg = file_header[DATE_keyword]
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
    Step 2: Group FITS Files whose header are to be Updated + Read Input File
    Step 3: Extract the details from 'file_telescopes' in a Pandas DataFrame
    Step 4: Updates the Header with changes in 'HeaderInfo.dat' & Appends AIRMASS etc. Details to the Header
    """
    # GUI Code for User Input
    DIR_FILES = eg.diropenbox('Enter the directory in which headers of files have to be updated:',
                              title='Path to the Directory', default=[os.path.join(os.getcwd(), 'preprocessed')])
    
    telescopename = eg.enterbox('Enter the Name of the Telescope from which the data was observed:',
                                title='Name of the Telescope', default=['HCT'])

    instrument = eg.enterbox('Enter the Instrument from which the data was observed:',
                                 title='Short Name of the Instrument', default=['HFOSC2'])

    input_file = eg.enterbox('Enter the name of the output file containing the header info:',
                             title='Path of the Input File', default=[os.path.join(os.getcwd(), file_headerinfo)])

    # Group FITS Files whose header are to be Updated + Read Input File
    if instrument == 'HFOSC2':
        prefix = 'fix'
    else:
        prefix = 'raw'

    list_files = group_similar_files('', common_text=os.path.join(DIR_FILES, prefix + '*.fits'))
    list_remove = group_similar_files('', common_text=os.path.join(DIR_FILES, '*.fits'), exceptions=prefix)
    for filename in list_remove:
        os.remove(filename)

    header_df = pd.read_csv(input_file, sep='\s*[|]\s*', comment='#', dtype='string', engine='python')
    header_df = header_df.set_index('FILENAME')

    # Extract the details from 'file_telescopes' in a Pandas DataFrame
    telescope_df = pd.read_csv(file_telescopes, sep='\s+', comment='#').set_index('ShortName')

    # Updates the Header with changes in 'HeaderInfo.dat' & Appends AIRMASS etc. Details to the Header
    for filename in list_files:
        update_header(filename, header_df)
        if identify_imagetype(filename):
            append_details(filename, telescopename, telescope_df)
        file_rename(filename, header_df)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Execute the Standalone Code
# ------------------------------------------------------------------------------------------------------------------- #
if __name__ == '__main__':
    main()
# ------------------------------------------------------------------------------------------------------------------- #