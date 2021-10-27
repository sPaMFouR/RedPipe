#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxx---------MODIFIES HEADER IN THE HFOSC2 FITS FILES AND REMOVES UNNECESSARY KEYWORDS---------xxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import glob
from astropy.io import fits
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
BIASSEC_keyword = 'BIASSEC'
CCDSEC_keyword = 'CCDSEC'
TRIMSEC_keyword = 'TRIMSEC'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates AIRMASS And Appends Respective Details In The Header
# ------------------------------------------------------------------------------------------------------------------- #
list_files = glob.glob('*.fits')

for file_name in list_files:
    with fits.open(file_name, mode='update') as hdulist:
        file_header = hdulist[0].header
        for keyword in [TRIMSEC_keyword, BIASSEC_keyword, CCDSEC_keyword]:
            if keyword in file_header.keys():
                file_header.remove(keyword, remove_all=True)
    
# ------------------------------------------------------------------------------------------------------------------- #
