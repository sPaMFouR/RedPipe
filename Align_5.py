#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxx----------------------------IMAGE ALIGNMENT---------------------------xxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import sys
import glob
import math
import shutil
import numpy as np
import pandas as pd
import easygui as eg
from pyraf import iraf
import imreg_dft as ird
from astropy.io import fits
from scipy.ndimage.interpolation import shift
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Telescope CCD Specifications
# ------------------------------------------------------------------------------------------------------------------- #
read_noise = 4.87
ccd_gain = 1.22
ccd_xsize = 2048
ccd_ysize = 2048
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Required IRAF Packages
# ------------------------------------------------------------------------------------------------------------------- #
iraf.noao(_doprint=0)
iraf.images(_doprint=0)
iraf.immatch(_doprint=0)
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


def remove_similar_files(common_text):
    """
    Removes similar files based on the string "common_text".
    Args:
        common_text : String containing partial name of the files to be deleted
    Returns:
        None
    """
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)


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
        with open(text_list, 'r+') as f:
            python_list = f.read().split()
            return python_list
    else:
        print ("ERROR: File '{0}' Not Found".format(text_list))
        sys.exit(1)


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
# Functions For Tasks In IRAF
# ------------------------------------------------------------------------------------------------------------------- #

def imexam_coords(sub_image, ref_coords, log_file):
    """
    Examines the image 'sub_image' at the coordinates mentioned in the file 'stars.coo' and logs the
    output onto the file 'log_imexam' using the IMEXAMINE task.
    Args:
        sub_image   : Subject image which needs to be aligned
        ref_coords  : Text file listing the reference coordinates of the chosen stars in the field
        log_file    : Name of the text list to record log of IMEXAM
    Returns:
        None
    """
    task = iraf.images.tv.imexam
    task.unlearn()

    task.logfile = log_file                         # Log File To Record Output Of IMEXAMINE Task
    task.keeplog = 'yes'                            # Log Output Results?
    task.defkey = 'a'                               # Default Key For Cursor x-y Input List
    task.imagecur = ref_coords                      # Image Display Cursor Input
    task.use_display = 'no'                         # Use The Image Display?

    task(input=sub_image, frame=1)


def imalign(textlist_images, ref_coords, prefix_str='a'):
    """
    Aligns the image 'sub_image' to the template image 'template_image' based on computed relative object shifts.
    Args:
        textlist_images : Text list of subject images which needs to be aligned with the reference image
        ref_coords      : Text file listing the reference coordinates of selected stars in the field
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    task = iraf.images.immatch.imalign
    task.unlearn()

    task.boxsize = 7                                # Box Size Used For Final Centering
    task.bigbox = 11                                # Box Size Used For Coarse Centering
    task.background = 'INDEF'                       # Absolute Reference Level For Centroid Calculation
    task.niterate = 3                               # Number Of Centering Iterations To Perform
    task.maxshift = 10                              # Maximum Permitted Shifts
    task.boundary_type = 'nearest'                  # Boundary Extension Type
    task.shiftimages = 'yes'                        # If 'yes', IMSHIFT Task Will Be Used To Align Images
    task.trimimages = 'no'                          # Whether Output Images Will Be Trimmed
    task.verbose = 'no'                             # Print Task Details?

    list_images = text_list_to_python_list(textlist_images)
    template_image = list_images[0]
    list_images = list_images[1:]

    for sub_image in list_images:
        output_filename = prefix_str + sub_image
        remove_file(output_filename)
        task(input=sub_image, reference=template_image, coords=ref_coords, output=output_filename)

    display_text("Alignment Using IMALIGN Completed")


def geomap(align_coords, align_db):
    """
    Commputes spatial transformation functions requred to map the reference coordinate system
    to the subject coordinate system.
    Args:
        align_coords  : Text file listing the reference and subject coordinates of selected stars in the field
        align_db      : Text file containing the coordinate transformations(produced by GEOMAP)
    Returns:
        None
    """
    task = iraf.images.immatch.geomap
    task.unlearn()

    task.fitgeometry = 'general'                        # Fitting Geometry To Be Used
    task.function = 'polynomial'                        # Type Of Analytic Surface To Be Fit
    task.xxorder = 2                                    # Order Of The Polynomials In 'x' And 'y'
    task.xyorder = 2
    task.yxorder = 2
    task.yyorder = 2
    task.reject = 3                                     # The Reject Limit In Units Of Sigma
    task.verbose = 'no'                                 # Print Task Messages?
    task.interactive = 'no'                             # Switch On Interactive Mode?

    task(input=align_coords, database=align_db, xmin=1, xmax=int(ccd_xsize), ymin=1, ymax=int(ccd_ysize))


def geotran(sub_image, align_db, align_coords, prefix_str):
    """
    Geometrically corrects the image 'sub_image' for distortions based on the coordinate transformations
    determined by GEOMAP.
    Args:
        sub_image      : Subject image which needs to be aligned
        align_db       : Text file containing the coordinate transformations(produced by GEOMAP)
        align_coords   : Text file listing the reference and subject coordinates of selected stars in the field
        prefix_str     : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    task = iraf.images.immatch.geotran
    task.unlearn()

    task.geometry = 'geometric'                         # Type Of Geometric Transformation
    task.xmin = 1                                       # Minimum X-Reference Value
    task.xmax = int(ccd_xsize)                          # Maximum X-Reference Value
    task.ymin = 1                                       # Minimum Y-Reference Value
    task.ymax = int(ccd_ysize)                          # Maximum Y-Reference Value
    task.interpolant = 'linear'                         # Interpolant Used For Rebinning The Image
    task.boundary = 'nearest'                           # Boundary Extension Type
    task.verbose = 'no'                                 # Print Task Messages?

    output_filename = prefix_str + sub_image
    remove_file(output_filename)
    task(input=sub_image, output=output_filename, database=align_db, transforms=align_coords)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Image Alignment
# ------------------------------------------------------------------------------------------------------------------- #

def modify_array(list_raw_array, scale_exponent=0.4):
    """
    Modifies arrays in the list 'list_raw_array' to bring them to a similar scale.
    Args:
        list_raw_array      : List of arrays which needs to be made usable for alignment procedure
        scale_exponent      : Exponent to which the array will be scaled
    Returns:
        list_usable_array   : List of array which has been made usable for alignment procedure
        """
    list_sliced_array = [array[300:1700, 300:1700] for array in list_raw_array]
    template_array = list_sliced_array[0]
    template_array[template_array < 0] = 0
    list_usable_array = []

    for array in list_sliced_array:
        array[array < 0] = 0
        exponent = (float(scale_exponent) * math.log10(np.median(template_array))) / math.log10(np.median(array))
        array_new = np.power(array, exponent)
        array_norm = array_new / np.mean(array_new)

        list_usable_array.append(array_norm)

    return list_usable_array


def translation(template_array, sub_array, scale_exponent=0.4):
    """
    Aligns array 'sub_array' to the template 'template_array' using translational transformation.
    Args:
        template_array  : Template array which will be used as a reference for translational aligning
        sub_array       : Subject array which will be aligned with the template_array using translation
        scale_exponent  : Exponent to which the array will be scaled
    Returns:
        mod_array       : Modified subject array aligned using translational transformation
    """
    list_array = modify_array(np.array([template_array, sub_array]), scale_exponent)
    template_arrayslice = list_array[0]
    sub_arrayslice = list_array[1]

    dict_shift = ird.translation(template_arrayslice, sub_arrayslice)
    list_shift = list(dict_shift['tvec'])
    mod_array = np.empty_like(sub_array)
    shift(sub_array, list_shift, output=mod_array, mode='nearest')

    print ("Translation In X & Y = {0}".format(list_shift))

    return mod_array


def similarity(template_array, sub_array, scale_exponent=0.4):
    """
    Aligns array 'sub_array' to the template 'template_array' using rotational transformation.
    Args:
        template_array  : Template array which will be used as a reference for rotational aligning
        sub_array       : Subject array (already corrected for translation)
        scale_exponent  : Exponent to which the array will be scaled
    Returns:
        mod_array       : Modified subject array aligned using rotational & translational transformation
    """
    list_array = modify_array([template_array, sub_array], scale_exponent)
    template_arrayslice = list_array[0]
    sub_arrayslice = list_array[1]

    dict_constraint = {'scale': [1.0, 0], 'tx': [0, 0], 'ty': [0, 0]}

    dict_rot = ird.similarity(template_arrayslice, sub_arrayslice, numiter=3, order=2, constraints=dict_constraint)
    mod_array = ird.transform_img(sub_array, scale=1.0, angle=dict_rot['angle'], tvec=dict_rot['tvec'], mode='nearest')

    print ("Rotation = {0}, Translation = {1}".format(dict_rot['angle'], dict_rot['tvec']))

    return mod_array


def write_fits(sub_image, mod_array, prefix_str='a'):
    """
    Writes the aligned array data 'mod_array' and header data from the original subject
    image 'sub_image' into a new aligned FITS file.
    Args:
        sub_image       : Subject image for obtaining the header details of the image
        mod_array       : Modified subject array to be written into the aligned FITS file
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    sub_header = fits.getheader(sub_image)
    output_filename = prefix_str + sub_image
    remove_file(output_filename)
    fits.writeto(filename=output_filename, data=mod_array, header=sub_header)


def generate_align_file(ref_coords, log_file, align_file):
    """
    Generates the text list 'align_coords' required by GEOMAP to calculate geometrical transformations
    based on data from reference coordinates file 'ref_coords' and log file 'log_imexam' from IMEXAMINE task.
    Args:
        ref_coords      : Text file listing the reference coordinates of the chosen stars in the field
        log_file        : Text file containing log from the IMEXAMINE task
        align_file      : Name of the text file for listing the reference and subject coordinates
    Returns:
        None
    """
    list_refcoords = text_list_to_python_list(ref_coords)
    no_of_stars = len(list_refcoords) / 2
    list_xref = list_refcoords[::2]
    list_yref = list_refcoords[1::2]

    data = pd.read_csv(log_file, sep='\s+', comment='#', header=None)
    list_xsub = data[0].tolist()
    list_ysub = data[1].tolist()

    remove_file(align_file)
    with open(align_file, 'w') as f:
        for index in range(0, no_of_stars):
            f.write(str(list_xref[index]) + ' ' + str(list_yref[index]) + ' ' +
                    str(list_xsub[index]) + ' ' + str(list_ysub[index]) + '\n')


def align_imgproc(textlist_images, scale_exponent=0.5, prefix_str='a'):
    """
    Executes the alignment procedure for the data in the image 'sub_image' with the template 'template_array'.
    Args:
        textlist_images : Text list of subject images which needs to be aligned with the reference image
        scale_exponent  : Exponent to which the array will be scaled
        prefix_str      : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_images = text_list_to_python_list(textlist_images)
    template_image = list_images[0]
    list_images = list_images[1:]

    shutil.copy(template_image, prefix_str + template_image)
    template_array = fits.getdata(filename=template_image)

    if len(list_images) > 0:
        for image in list_images:
            sub_array = fits.getdata(filename=image)
            mod_array = translation(template_array, sub_array, scale_exponent)

#             mod_array = translation(template_array, mod_array, scale_exponent)
            mod_array = similarity(template_array, mod_array, scale_exponent)
#            mod_array = similarity(template_array, mod_array, scale_exponent)

            write_fits(image, mod_array, prefix_str)
    else:
        print ("ERROR: Too Few Images To Align")

    display_text("Alignment Using Image Registration Completed")


def align_iraf(textlist_images, log_imexam, ref_coords, align_coords, align_db, itern=2, prefix_str='a'):
    """
    Executes the alignment procedure tasks for the images in the list 'textlist_images' given the reference
    coordinates in the file 'align_coords'. This is done using the GEOMAP & GEOTRAN tasks in DAOPHOT Package.
    Args:
        textlist_images : Text list of subject images which needs to be aligned with the reference image
        ref_coords       : Text list containing the reference coordinates of the chosen stars in the field
        log_imexam       : Name of the text file to record log of IMEXAM
        align_coords     : Name of the text file to list the reference & subject coordinates of the chosen stars
        align_db         : Name of the text file to record the coordinate transformations(produced by GEOMAP)
        itern            : No. of iterations for which the task is to be run
        prefix_str       : Prefix to distinguish the aligned FITS file from the original FITS file
    Returns:
        None
    """
    list_images = text_list_to_python_list(textlist_images)
    template_image = list_images[0]
    list_images = list_images[1:]
    shutil.copy(template_image, prefix_str + template_image)

    if len(list_images) > 0 and len(text_list_to_python_list(ref_coords)) != 0:
        for image in list_images:
            print ("File Name: {0}: ".format(image))
            output_filename = prefix_str + image
            list_temp = []
            for value in range(0, int(itern)):
                temp_prefix = str(value)
                remove_file(log_imexam)
                imexam_coords(image, ref_coords=ref_coords, log_file=log_imexam)
                generate_align_file(ref_coords=ref_coords, log_file=log_imexam, align_file=align_coords)

                geomap(align_coords=align_coords, align_db=align_db)
                geotran(image, align_db=align_db, align_coords=align_coords, prefix_str=temp_prefix)
                image = temp_prefix + image

                if value != int(itern) - 1:
                    list_temp.append(image)
                else:
                    shutil.move(image, output_filename)

                remove_file(align_db)
                remove_file(log_imexam)
                remove_file(align_coords)

            for temp_image in list_temp:
                remove_file(temp_image)

    elif len(list_images) > 0 and len(text_list_to_python_list(ref_coords)) == 0:
        print ("ERROR: Reference Coordinates Not Specified In The File {0}".format(ref_coords))
        sys.exit(1)

    else:
        print ("ERROR: Too Few Images To Align")
        sys.exit(1)

    display_text("Alignment Using GEOMAP & GEOTRAN Completed")


def check_aligned(textlist_images, log_align, ref_coords='stars.coo', log_imexam='log_imexam'):
    """
    Check whether the files are aligned by checking coordinate shift in the subject from
    reference coords of the chosen stars.
    Args:
        textlist_images  : Text file listing the names of the files for which alignment has to be checked
        log_align        : Name of the text file for listing the shifts after aligning procedure
        ref_coords       : Text file listing the reference coordinates of the chosen stars in the field
        log_imexam       : Name of the text file to record log of IMEXAMINE task
    Returns:
        None
    """
    list_images = text_list_to_python_list(textlist_images)
    list_refcoords = text_list_to_python_list(ref_coords)

    no_of_stars = len(list_refcoords) / 2
    list_xref = list_refcoords[::2]
    list_yref = list_refcoords[1::2]

    remove_file(log_imexam)
    for image in list_images:
        print ("File Name: {0}: ".format(image))
        imexam_coords(image, ref_coords, log_imexam)

    data_sub = pd.read_csv(log_imexam, sep='\s+', comment='#', header=None)

    remove_file(log_align)
    with open(log_align, 'w') as fout:
        fout.write("{0:>7s}{1:>9s}{2:>11s}{3:>11s}{4:>11s}{5:>11s}{6:>9s}\n\n".format
                   ("Star_ID", "X-Ref", "Y-Ref", "X-Img", "Y-Img", "X-Err", "Y-Err"))
        for i in range(0, len(list_images)):
            fout.write("# IMAGE #{0} - {1}\n".format(i + 1, list_images[i]))
            for j in range(0, no_of_stars):
                fout.write("{0:>4}{1:>13.2f}{2:>11.2f}{3:>11.2f}{4:>11.2f}{5:>9.2f}{6:>9.2f}\n".format
                           (j + 1, float(list_xref[j]), float(list_yref[j]),
                            float(data_sub.loc[i * no_of_stars + j, 0]), float(data_sub.loc[i * no_of_stars + j, 1]),
                            (float(data_sub.loc[i * no_of_stars + j, 0]) - float(list_xref[j])),
                            (float(data_sub.loc[i * no_of_stars + j, 1]) - float(list_yref[j]))))
            fout.write('\n')

    display_text("Log Of Alignment Has Been Generated")

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# ------------------------------------------------------------------------------------------------------------------- #
# rmv_files = eg.boolbox('Remove Residual Files From Previous Run?', title='Remove Residual Files', choices=['Yes', 'No'])
# ctext = eg.enterbox('Enter The Common Text Of Files To Be Aligned?', title='Files To Be Aligned', default='*.fits')
# bool_coarse = eg.boolbox('Perform Aligning Using Image Processing?', title='Coarse Aligning', choices=['Yes', 'No'])
# bool_precise = eg.boolbox('Perform Aligning Using IRAF?', title='Precise Aligning', choices=['Yes', 'No'])
# ref_file = eg.enterbox('Enter The File With Reference Coordinates:', title='Reference Coordinates', default='stars.coo')
#
# choice_precise = False
# if bool_precise:
#     choice_precise = eg.boolbox('Align Using Which IRAF Task?', title='Aligning Task To Be Used',
#                                 choices=['GEOMAP & GEOTRAN', 'IMALIGN'])
#
rmv_files = True
ctext = '*.fits'
bool_coarse = True
bool_precise = True
choice_precise = True
ref_file = 'stars.coo'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Removes Residual Files From Previous Run
# ------------------------------------------------------------------------------------------------------------------- #
if rmv_files:
    remove_similar_files('list_*')
    remove_similar_files('a_*.fits')
    remove_similar_files('ca_*.fits')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups & Performs Aligning Of Images
# ------------------------------------------------------------------------------------------------------------------- #
prefix = 'ca_'
text_list = 'list_align'
group_similar_files(text_list, common_text=ctext)

if bool_coarse and not bool_precise:
    align_imgproc(textlist_images=text_list, scale_exponent=0.4, prefix_str=prefix)
elif bool_coarse and bool_precise:
    align_imgproc(textlist_images=text_list, scale_exponent=0.4, prefix_str='a_')
    prefix = 'c'
    text_list = 'list_align_coarse'
    group_similar_files(text_list, 'a_*.fits')

if bool_precise:
    if choice_precise:
        align_iraf(textlist_images=str(text_list), log_imexam='log_imexam', ref_coords=str(ref_file),
                   align_coords='align.coo', align_db='align.db', prefix_str=prefix)
    else:
        imalign(textlist_images=str(text_list), ref_coords=str(ref_file), prefix_str=prefix)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Logs The Precision Of Images Aligned In The File 'log_alignment'
# ------------------------------------------------------------------------------------------------------------------- #
group_similar_files('list_aligned', prefix + ctext)
check_aligned('list_aligned', 'log_after_alignment')

# try:
#     task(input=str(sub_image), frame=1)
# except iraf.IrafError, e:
#     print e
#     pass
# ------------------------------------------------------------------------------------------------------------------- #
