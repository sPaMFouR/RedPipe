#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxx-------------------------PRINT HEADER INFORMATION----------------------xxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import shutil
import numpy as np
import pandas as pd
import easygui as eg
from astropy.io import fits
from datetime import datetime
from astropy.table import Table
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
DATE_keyword = 'DATE-OBS'
NAXIS1_keyword = 'NAXIS1'
NAXIS2_keyword = 'NAXIS2'

OBJECT_keyword = 'OBJECT'
RA_keyword = 'RA'
DEC_keyword = 'DEC'
FILTER_keyword = 'FILTER'
GRISM_keyword = 'GRISM'
EXPTIME_keyword = 'EXPTIME'

UT_keyword = 'UT'
JD_keyword = 'JD'
AIRMASS_keyword = 'AIRMASS'

header_keys = ['FILENAME', DATE_keyword, NAXIS1_keyword, NAXIS2_keyword, OBJECT_keyword, RA_keyword, DEC_keyword,
               FILTER_keyword, GRISM_keyword, EXPTIME_keyword, UT_keyword, JD_keyword, AIRMASS_keyword]
header_format = {JD_keyword: '{0:>.5f}', AIRMASS_keyword: '{0:>.3f}'}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def remove_similar_files(common_text):
    """
    Removes similar files based on the string 'common_text'.
    Args:
        common_text : String containing partial name of the files to be deleted
    Returns:
        None
    """
    for filename in glob.glob(common_text):
        try:
            os.remove(filename)
        except OSError:
            pass


def create_dir(dir_path):
    """
    Creates a directory specified by 'dir_path' if the directory does not exist.
    Args:
        dir_source  : Path of the directory to be created
    Returns:
        None
    """
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
    else:
        remove_similar_files(os.path.join(dir_path, '*'))


def copy_files(list_files, dir_source, prefix='raw_'):
    """
    Copies list of files specified by 'list_files' onto the directory 'based on the string 'common_text'.
    Args:
        list_files  : Python list containing the names of the grouped files
        dir_source  : Path of the source directory where the files have to be copied
    Returns:
        None
    """
    for filename in list_files:
        print (filename, os.path.join(dir_source, prefix + filename.split('/')[-1]))
        shutil.copy(filename, os.path.join(dir_source, prefix + filename.split('/')[-1]))
        

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
        for filename in glob.glob(common_text):
            for text in list_exception:
                test = re.search(text, filename)
                if test:
                    try:
                        list_files.remove(filename)
                    except ValueError:
                        pass

    list_files.sort()
    if len(text_list) != 0:
        with open(text_list, 'w') as f:
            for filename in list_files:
                f.write(filename + '\n')

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
# Function for Formatting and Reading the Header of the FITS file
# ------------------------------------------------------------------------------------------------------------------- #

def format_dateobs(dateobs):
    """
    Formats the value of DATE_keyword in the FITS header to account for time in milliseconds.
    Args:
        dateobs : Value of the DATE_keyword that is to be modified
    Returns:
        datenew : Modified value of the DATE_keyword which accounts for time in milliseconds.
    """
    datetime_master = datetime.strptime(dateobs, '%Y-%m-%dT%H:%M:%S.%f')
    datenew = datetime_master.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]
    return datenew


def format_header(list_files, instrument):
    """
    Formats and homogenizes the header of the files in the list 'list_files'.
    Args:
        list_files : List of FITS files whose header details have to be formatted
        instrument : Name of the Instrument from which the data was observed
    Returns:
        None
    """
    for filename in list_files:
        with fits.open(filename, mode='update') as hdulist:
            header = hdulist[0].header

            if instrument == 'HFOSC2':
                header[DATE_keyword] = format_dateobs(header['DATE-AVG'])
            elif DATE_keyword in header.keys():
                header[DATE_keyword] = format_dateobs(header[DATE_keyword])
            else:
                display_text("ERROR: DATE keyword not found in the Header")
                sys.exit(1)


def fix_header(list_files, dir_source, extnheader=0, extndata=1, prefix='fix_'):
    """
    FIX multi-extension headers of FITS files specified by the list 'list_files' as in the case of 'HFOSC2'.
    The extension of the header and the data is specified by the keywords 'extnheader' and 'extndata'.
    The new file will have the prefix 'prefix' appended and copied to the directory 'dir_source'.
    Args:
        list_files : List of FITS files whose header details have to be read
        dir_source : Path of the source directory where the files have to be copied
        extnheader : Extension of the FITS files which stores the complete header
        extndata   : Extension of the FITS files which stores the data
        prefix     : Prefix to be appended to the Fixed FITS files
    Returns:
        None
    """
    for filename in list_files:
        with fits.open(filename, mode='readonly') as hdulist:
            hdulist.info()
            file_header = hdulist[extnheader].header
            file_data = hdulist[extndata].data

            hdunew = fits.PrimaryHDU(data=file_data, header=file_header)
            hdunew.writeto(os.path.join(dir_source, prefix + filename.split('/')[-1]))

    display_text("Hurray: FITS Files have been Successfully Fixed to a Single Extension Format.")


def read_headerinfo(list_files, outfile, extn=0):
    """
    Read headers of the FITS files in the list 'list_files' from the extension specified by 'extn'.
    Args:
        list_files : List of FITS files whose header details have to be read
        outfile    : Name of the output file containing the header info
        extn       : Extension of the FITS files which stores the complete header
    Returns:
        None
    """
    header_df = pd.DataFrame(index=range(len(list_files)), columns=header_keys)

    for idx, filename in enumerate(list_files):
        with fits.open(filename, mode='readonly') as hdulist:
            file_header = hdulist[extn].header
            for header in header_keys:
                if header in file_header.keys():
                    if header == 'FILENAME':
                        header_df.loc[idx, header] = filename.split('/')[-1]
                    elif header not in header_format:
                        header_df.loc[idx, header] = file_header[header]
                    else:
                        header_df.loc[idx, header] = header_format[header].format(float(file_header[header]))
    print (header_df)
    header_df = header_df.replace(np.nan, 'INDEF')
    tab_header = Table.from_pandas(header_df)
    tab_header.write(outfile, format='ascii.fixed_width', delimiter=' | ', overwrite=True)

    print (" ")
    print (tab_header)
    display_text("Hurray: Header Info Sucessfully Printed.")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Main Function
# ------------------------------------------------------------------------------------------------------------------- #

def main():
    """
    Step 1: GUI Code for User Input
    Step 2: Group FITS Files whose header are to be read
    Step 3: Create a Directory where Pre-Processed Data will be saved
    Step 4: Perform Fixing of File Header Data Unit (Only for HFOSC2)
    Step 5: Print Header Info of the files in the directory
    """
    # GUI Code for User Input
    DIR_FILES = eg.diropenbox('Enter the directory from which header of files are to be read:',
                              title='Path to the Directory', default=[os.getcwd()])

    common_text = eg.enterbox('Enter the common text of Files to be read:',
                              title='Common Text of Files', default=['*.fits'])

    instrument = eg.enterbox('Enter the Instrument from which the data was observed:',
                                 title='Short Name of the Instrument', default=['HFOSC2'])

    output_file = eg.enterbox('Enter the name of the output file containing the header info:',
                              title='Name of the Output File', default=['HeaderInfo.dat'])

    # Group FITS Files whose header are to be read
    list_files = group_similar_files('', common_text=os.path.join(DIR_FILES, common_text))
    if len(list_files) == 0:
        display_text("ERROR: No FITS files found in the Specified Directory")
        sys.exit(1)
    
    # Create a Directory where Pre-Processed Data will be saved
    DIR_SAVE = os.path.join(DIR_FILES, 'preprocessed')
    create_dir(DIR_SAVE)
    copy_files(list_files, DIR_SAVE, prefix='raw_')

    # Perform Fixing of File Header Data Unit (Only for HFOSC2)
    if instrument == 'HFOSC2':
        fix_header(list_files, DIR_SAVE, extnheader=0, extndata=1, prefix='fix_')
        list_files = group_similar_files('', common_text=os.path.join(DIR_SAVE, 'fix_*' + common_text))
        format_header(list_files, instrument)
    else:
        list_files = group_similar_files('', common_text=os.path.join(DIR_SAVE, 'raw_*' + common_text))
        pass

    # Print Header Info of the files in the directory
    read_headerinfo(list_files, outfile=os.path.join(DIR_FILES, output_file), extn=0)
    display_text("Important: Read Instructions in 'RenameInfo.dat' before formatting 'HeaderInfo.dat'.")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Execute the Code
# ------------------------------------------------------------------------------------------------------------------- #
if __name__ == '__main__':
    main()
# ------------------------------------------------------------------------------------------------------------------- #
