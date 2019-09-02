# /usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxx----------------------------COMBINE IMAGES----------------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import glob
from pyraf import iraf
from astropy.io import fits
import dateutil.parser as dparser
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Telescope CCD Specifications & Image Header Keywords
# ------------------------------------------------------------------------------------------------------------------- #
read_noise = 4.87
ccd_gain = 1.22
EXPTIME_keyword = 'EXPTIME'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Name Of The Object
# ------------------------------------------------------------------------------------------------------------------- #
OBJECT_name = '2018hna'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For File Handling
# ------------------------------------------------------------------------------------------------------------------- #

def remove_file(file_name):
    """
    Removes the file "file_name" in the constituent directory.
    Args:
         file_name  : Name of the file to be removed from the current directory
    Returns:
        None
    """
    try:
        os.remove(file_name)
    except OSError:
        pass


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
        print ("Error : File '{0}' Not Found".format(text_list))


def calculate_exptime(textlist_images):
    """
    Calculates total exposure for the images in the list 'list_images'.
    Args:
        textlist_images : Text list of subject images which needs to be combined
    Returns:
        total_exptime   : Total exposure time of the combined image
    """
    list_images = text_list_to_python_list(textlist_images)

    total_exptime = 0
    if len(list_images) != 0:
        for image in list_images:
            file_header = fits.getheader(image)
            exptime = file_header[EXPTIME_keyword]
            total_exptime += int(exptime)
    else:
        print ("No Images In The Text List {0}".format(textlist_images))

    return total_exptime


def edit_exptime(comb_image, total_exptime):
    """
    Calculates total exposure for the images in the list 'list_images'.
    Args:
        comb_image     : Image whose header has to be edited
        total_exptime  : Total exposure time of the combined image
    Returns:
        None
    """
    hdulist = fits.open(comb_image, mode='update')
    hdulist[0].header[EXPTIME_keyword] = int(total_exptime)
    hdulist.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def imcombine(textlist_images, type_combine='sum'):
    """
    Combines images in the text list 'text_list_images' using algorithms based on the combine
    type mentioned in the variable 'type_combine'.
    Args:
        textlist_images  : Text list of subject images which needs to be combined
        type_combine     : Type of combining operation to be performed on pixels
    Returns:
        None
    """
    task = iraf.images.immatch.imcombine
    task.unlearn()

    task.combine = type_combine                     # Type Of Combining Operation Performed On Pixels
    task.reject = 'none'                            # Type Of Rejection Operation Performed On Pixels
    task.project = 'no'                             # Combine Across The Highest Dimension Of The Image?
    task.rdnoise = float(read_noise)                # CCD Readout Noise (In e-)
    task.gain = float(ccd_gain)                     # CCD Gain (In e-/ADU)

    task(input='@' + textlist_images, output=textlist_images[5:] + ".fits")
    total_exptime = calculate_exptime(textlist_images)
    edit_exptime(textlist_images[5:] + '.fits', total_exptime)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Images To Be Combined - [ca_jul30_fbs_object-v1.fits]
# ------------------------------------------------------------------------------------------------------------------- #
list_dates = []
for file_name in group_similar_files('', common_text='*.fits'):
    temp_name = file_name.split(OBJECT_name)[0]
    date = dparser.parse(temp_name, fuzzy=True)
    date = date.strftime('%Y-%m-%d')
    list_dates.append(date)

list_dates = set(list_dates)
list_filters = ['U', 'B', 'V', 'R', 'I']
list_patterns = ['ca_' + date + '_cfbs_' + OBJECT_name + '-' + band.lower() for date in list_dates for band in list_filters]

list_list_comb = []
for pattern in list_patterns:
    if len(group_similar_files('', pattern + '*')) != 0:
        pattern_files = group_similar_files('list_' + pattern, common_text=pattern + '*')
        if len(pattern_files) > 1:
            list_list_comb.append('list_' + pattern)
print list_list_comb
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Combines Images In The List 'list_comb'
# ------------------------------------------------------------------------------------------------------------------- #
for list_comb in list_list_comb:
    imcombine(textlist_images=list_comb, type_combine='sum')

# for file_name in group_similar_files('', 'list_ca*'):
#     remove_file(file_name)
# ------------------------------------------------------------------------------------------------------------------- #
