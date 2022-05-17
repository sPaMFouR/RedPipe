#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx------------------------PERFORMS PRE-PROCESSING----------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
import numpy as np
import easygui as eg
import ccdproc as ccdp
from astropy import units as u
from astropy.stats import mad_std
from astropy.nddata import CCDData
from photutils import detect_sources
from astropy.nddata.blocks import block_replicate
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# GLOBAL VARIABLES
# ------------------------------------------------------------------------------------------------------------------- #s
READ_NOISE = 4.87  # Units in electrons/photon
GAIN = 1.22 * u.electron / u.adu
SATURATION = 800000 # Units in electrons/photon
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
# Functions for Performing Pre-Processing
# ------------------------------------------------------------------------------------------------------------------- #

def inv_median(a):
    """
    Used as a scaling function in flat combine task.
    """
    return 1 / np.median(a)


def bias_combine(textlist_bias, master_bias='mbias.fits'):
    """
    Combines bias files to make a master bias.
    Args:
        textlist_bias - A python list object with paths/names to the individual files.
        master_bias   - Output master bias file (default 'mbias.fits')
    Returns:
        None
    """
    new_list = [CCDData.read(filename, unit=u.adu) for filename in textlist_bias]

    combined_bias = ccdp.combine(img_list=new_list, method='average',
                                sigma_clip=True,sigma_clip_low_thresh=5,
                                sigma_clip_high_thresh=5, 
                                sigma_clip_func=np.ma.median,
                                sigma_clip_dev_func=mad_std)

    combined_bias.meta['combined'] = True
    combined_bias.data = combined_bias.data.astype('float32')
    combined_bias.write(master_bias, hdu_mask=None, hdu_uncertainty=None)


def subtract_bias(textlist_files, master_bias='mbias.fits', prefix_str='bs_'):
    """
    Subtract bias from the given files list
    Args:
        textlist_files - A python list object with paths/names to the individual files.
        master_bias    - Master bias used for subtraction (default 'mbias.fits')
        prefix_str     - String to be prefixed for subtracted files
    Returns:
        None
    """
    master = CCDData.read(master_bias, unit=u.adu)

    for filename in textlist_files:
        ccd = CCDData.read(filename, unit=u.adu)
        bias_subtracted = ccdp.subtract_bias(ccd=ccd, master=master)
        bias_subtracted.meta['biassub'] = True
        bias_subtracted.data = bias_subtracted.data.astype('float32')
        bias_subtracted.write(prefix_str+filename, hdu_mask=None, hdu_uncertainty=None)


def flat_combine(textlist_files, band, outfile='mflat'):
    """
    Combines multiple flats for a given input files list

    Args:
        textlist_files - A python list object with paths/names to the individual files.
        band           - Filter name associated with the flat files
        outfile        - Master flat name (default outfile + band + '.fits')
    Returns:
        None
    """
    new_list = [CCDData.read(filename, unit=u.adu) for filename in textlist_files]

    combined_flat = ccdp.combine(img_list=new_list, method='average', scale=inv_median,
                                sigma_clip=True, sigma_clip_low_thresh=5, 
                                sigma_clip_high_thresh=5, sigma_clip_func=np.ma.median,
                                sigma_clip_dev_func=mad_std)

    combined_flat.meta['combined'] = True
    file_name = outfile + band + '.fits'
    combined_flat.data = combined_flat.data.astype('float32')
    combined_flat.write(file_name, hdu_mask=None, hdu_uncertainty=None)


def flat_correction(textlist_files, master_flat, prefix_str='f'):
    """
    To flat field the science images using master flat
    Args:
        textlist_files - A python list object with paths/names to the individual files.
        master_flat    - Master flat used to flat field the science image
        prefix_str     - String to be added to newly created sceince image
    Returns:
        None
    """
    master = CCDData.read(master_flat, unit=u.adu)

    for filename in textlist_files:
        ccd = CCDData.read(filename, unit=u.adu)
        flat_corrected = ccdp.flat_correct(ccd=ccd, flat=master, min_value=0.9)
        flat_corrected.meta['flatcorr'] = True
        flat_corrected.data = flat_corrected.data.astype('float32')

        flat_corrected.write(prefix_str + filename, hdu_mask=None, hdu_uncertainty=None)


def cosmic_ray_corr(textlist_files, prefix_str='c'):
    """
    Gain correction and Cosmic ray correction using LA Cosmic method.
    Args:
        textlist_files - A python list object with paths/names to the individual files.
        prefix_str     - String appended to the name of newly created file
    Returns:
        None
    """
    for filename in textlist_files:

        file_corr = CCDData.read(filename, unit=u.adu)
        file_corr = ccdp.gain_correct(file_corr, gain=GAIN)
        new_ccd = ccdp.cosmicray_lacosmic(file_corr, readnoise=READ_NOISE, sigclip=7,
                                             satlevel=SATURATION, niter=4,
                                             gain_apply=False, verbose=True)
        new_ccd.meta['crcorr'] = True
        new_ccd.data = new_ccd.data.astype('float32')
        new_ccd.write(prefix_str + filename, hdu_mask=None, hdu_uncertainty=None)

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
    DIR_FILES = eg.enterbox('Enter the directory in which preprocessing has to be performed:',
                            title='Enter the Directory Path', default=[os.path.join(os.getcwd(), 'preprocessed')])
    
    telescopename = eg.enterbox('Enter the Name of the Telescope from which the data was observed:',
                                title='Enter the Name of the Telescope', default=['HCT'])

    instrument = eg.enterbox('Enter the Instrument from which the data was observed:',
                                 title='Enter the Short Name of the Instrument', default=['HFOSC2'])

    input_file = eg.enterbox('Enter the name of the output file containing the header info:',
                             title='Enter the Name of the Output File', default=['HeaderInfo.dat'])

    # Group FITS Files whose header are to be Updated + Read Input File

    # Extract the details from 'file_telescopes' in a Pandas DataFrame
    telescope_df = pd.read_csv(file_telescopes, sep='\s+', comment='#').set_index('ShortName')
    list_files = group_similar_files('', common_text=os.path.join(DIR_FILES, prefix + '*.fits'))



# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Execute the Standalone Code
# ------------------------------------------------------------------------------------------------------------------- #
if __name__ == '__main__':
    main()
# ------------------------------------------------------------------------------------------------------------------- #