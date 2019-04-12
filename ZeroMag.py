###########################################################################
#      Calculation Of Instrumental Magntiudes & Zero Point Correction     #
###########################################################################

#####---Read The 'output_mag2' File Generated After 'txdump' Task---#####

f=open('output_mag2')
data=f.read()
f.close()


#####---Make The Data Usable By Splitting The Data Into A List Of Strings---#####

data1=data.split()


#####---Initialize Lists To Be Printed In The Output File"---#####

star_id = []
image_name = []
xcenter = []
ycenter = []
ifilter = []
sky_counts = []
sigma_counts = []
aperture_1 = []
phot_mag_1 = []
phot_mag_err_1 = []
aperture_2 = []
phot_mag_2 = []
phot_mag_err_2 = []
aperture_correction = []
aperture_correction_mean = []
airmass = []
instr_mag = []
instr_mag_err = []


#####---Extinction Coefficients In Mag For Hanle In Different Bands---#####

# FILTER  EXTINCTION_MEAN   EXTINCTION_ERROR      
#   U        0.36              0.07 
#   B        0.21              0.04
#   V        0.12              0.04
#   R        0.09              0.04
#   I        0.05              0.03


#####---Initialize Global Variables---#####

count = 0
eeta = {'7BesU' : 0.36, '6BesB' : 0.21, '5BesV' : 0.12, '4BesR' : 0.09, '3BesI' : 0.05}
eeta_err = {'7BesU' : 0.07, '6BesB' : 0.04, '5BesV' : 0.04, '4BesR' : 0.04, '3BesI' : 0.03}
length_data = len(data1)
rows = length_data/14


for i in range(0, rows):

    #####---Calculate The Aperture Correction For Each Star In An Image---#####

    if str(data1[9 +i*14]) != 'INDEF' and str(data1[10 +i*14]) != 'INDEF':
        mag_phot_1 = float(data1[9 +i*14])
        mag_phot_2 = float(data1[10 +i*14])
        mag_err_phot_1 = float(data1[11 +i*14])
        mag_err_phot_2 = float(data1[12 +i*14])
        correction_aperture = mag_phot_1 - mag_phot_2
    else:
        mag_phot_1 = 'INDEF'
        mag_err_phot_1 = 'INDEF'
        mag_phot_2 = 'INDEF'
        mag_err_phot_2 = 'INDEF'
        mag_instr = 'INDEF'
        mag_err_instr = 'INDEF'
        correction_aperture = 'INDEF'

    #####---Append The Data To The Respective Lists---#####

    star_id.append(float(data1[0 + i*14]))
    image_name.append(str(data1[1 + i*14]))
    xcenter.append(float(data1[2 + i*14]))
    ycenter.append(float(data1[3 + i*14]))
    ifilter.append(str(data1[4 + i*14]))
    sky_counts.append(float(data1[5 + i*14]))
    sigma_counts.append(float(data1[6 + i*14]))
    aperture_1.append(float(data1[7 + i*14]))
    aperture_2.append(float(data1[8 + i*14]))
    airmass.append(float(data1[13 + i*14]))

    phot_mag_1.append(mag_phot_1)
    phot_mag_err_1.append(mag_err_phot_1)
    phot_mag_2.append(mag_phot_2)
    phot_mag_err_2.append(mag_err_phot_2)
    aperture_correction.append(correction_aperture)


#####---Calculate Mean Aperture Correction For Each Image---#####

star_count = int(star_id[-1])

for j in range(0, rows):

    x = j // star_count
    aperture_correction_sum = 0

    if aperture_correction[j] != 'INDEF':
        for i in range(x*star_count, x*star_count + star_count):
            aperture_correction_sum += aperture_correction[i]

        aperture_correction_mean.append(float(aperture_correction_sum/star_count))

    else:
        aperture_correction_mean.append('INDEF')


#####---Caclulate The Instrumental Magnitude With Error---#####

for i in range(0, rows):

    mag_instr = phot_mag_1[i] - airmass[i]*eeta[str(data1[4 + i*14])] - aperture_correction_mean[i]
    mag_err_instr = phot_mag_err_1[i]
    instr_mag.append(mag_instr)
    instr_mag_err.append(mag_err_instr)


#####---Print The Calculated Data Along With Original Data In The Output File---#####

f2=open("output_instr_mag2", "w")

f2.write(str('STAR_ID') + "        " + str('IMAGE_NAME') + "         " + str('FILTER') + "  " + str('SKY_COUNTS') + "  " + str('AIRMASS') + "  " + str('APER_1') + "   " + str('APER_2') + "    " + str('MAG_1') + "    " + str('MAG_2') + "     " + str('APER_CORR') + "   " + str('AP_CO_MEAN') + "   " + str('INSTR_MAG') + "  " + str('INSTR_MAG_ERR') + "\n" + "\n")

for i in range(0, rows):

    if type(phot_mag_1[i]) == float and type(phot_mag_2[i]) == float:
        f2.write(str("%3.0f" % star_id[i]) + "     " + str("%23s" % image_name[i]) + "  " + str("%6s" % ifilter[i]) + "   " +
        str("%8.4f" % sky_counts[i]) + "   " + str("%8.6f" % airmass[i]) + "   " + str("%4.2f" % aperture_1[i]) + "    " +
        str("%4.2f" % aperture_2[i]) + "   " + str("%7.3f" % phot_mag_1[i]) + "   "  + str("%7.3f" % phot_mag_2[i]) + "    " +
        str("%8.6f" % aperture_correction[i]) + "     " + str("%8.6f" % aperture_correction_mean[i]) + "   " +
        str("%8.3f" % instr_mag[i]) + "    " + str("%7.3f" % instr_mag_err[i]) + "\n")

    else:
        f2.write(str("%3.0f" % star_id[i]) + "     " + str("%23s" % image_name[i]) + "   " + str("%6s" % ifilter[i]) + "   " +
        str("%8.4f" % sky_counts[i]) + "   " + str("%8.6f" % airmass[i]) + "   " + str("%4.2f" % aperture_1[i]) + "    " +
        str("%4.2f" % aperture_2[i]) + "   " + str("%7s" % phot_mag_1[i]) + "   "  + str("%7s" % phot_mag_2[i]) + "    " +
        str("%10s" % aperture_correction[i]) + "     " + str("%8s" % aperture_correction_mean[i]) + "   " + str("%8s" % instr_mag[i]) +
        "    " + str("%7s" % instr_mag_err[i]) + "\n")

    if star_id[i] == int(star_count):
        f2.write("\n")

f2.close()


# list_instr_mag = []
# list_instr_mag_err = []
# for ap_index in range(0, len(list_mag)):
#     for index in range(0, rows):
#         mag_instr = list_mag[ap_index][0][index] - airmass[index] * eeta[ifilter[index]]
#         mag_err_instr = list_mag[ap_index][1][index]
#         list_instr_mag.append(mag_instr)
#         list_instr_mag_err.append(mag_err_instr)
#
# with open(str(file_name), 'w') as fout:
#     fout.write(str('STAR_ID') + " " * 8 + str('IMAGE_NAME') + " " * 9 + str('FILTER') + " " * 2 +
#                str('SKY_COUNTS') + " " * 2 + str('AIRMASS') + " " * 2 + str('APER') + " " * 3 +
#                str('MAG') + " " * 4 + str('MAG_ERR') + " " * 5 + str('INSTR_MAG') + " " * 2 +
#                str('INSTR_MAG_ERR') + "\n" + "\n")
#
#     for index in range(0, rows):
#         if type(list_mag[index]) == float:
#             fout.write(str("%3.0f" % star_id[index]) + " " * 5 + str("%23s" % img_name[index]) + " " * 2 +
#                        str("%6s" % ifilter[index]) + " " * 3 + str("%8.4f" % sky_counts[index]) + " " * 3 +
#                        str("%6.4f" % airmass[index]) + " " * 3 + str("%5.2f" % list_aper[index]) + " " * 4 +
#                        str("%7.4f" % list_mag[index]) + " " * 3 + str("%7.4f" % list_instr_mag[index]) + " " * 4 +
#                        str("%6.4f" % list_instr_mag_err[index]) + "\n")
#         else:
#             fout.write(str("%3.0f" % star_id[index]) + " " * 5 + str("%23s" % img_name[index]) + " " * 2 +
#                        str("%6s" % ifilter[index]) + " " * 3 + str("%8.4f" % sky_counts[index]) + " " * 3 +
#                        str("%6.4f" % airmass[index]) + " " * 3 + str("%5.2f" % list_aper[index]) + " " * 4 +
#                        str("%7s" % list_mag[index]) + " " * 3 + str("%7s" % list_instr_mag[index]) + " " * 4 +
#                        str("%6s" % list_instr_mag_err[index]) + "\n")
#
#         if star_id[index] == int(star_count):
#             fout.write("\n")

import re
import sys


def get_file_dmnsn(file_name, title_rows=0):
    """
    Finds out the number of rows and columns in a text file.
    Args:
        file_name   : Text file whose dimensions are to be obtained
        title_rows  : No. of rows used as title description
    Returns:
        rows        : Number of rows in the text file
        columns     : Number of columns in the text file
    """
    with open(str(file_name), 'r') as f:
        columns = len(f.readline().rstrip().split())
    if columns == 0:
        print ("Error: '" + str(file_name) + "' Is An Empty File ")
        sys.exit(1)

    with open(str(file_name), 'r') as f:
        for i in range(0, int(title_rows)):
            f.readline()
        length_data = len(f.read().split())
        rows = length_data / columns

    return rows, columns


def zip_list(list_lists):
    """
    Combines all the lists in a list to a single list element-wise.
    Args:
        list_lists  : List of all the lists which needs to be zipped
    Returns:
        new_list    : List with the combined elements
    """
    new_list = []
    for index in range(0, len(list_lists[0])):
        value = ''
        for val in range(0, len(list_lists)):
            value += list_lists[val][index]
        new_list.append(str(value))

    return new_list


def read_mag(file_name, col_nos, fmt='{:>8}', title_rows=0):
    """
    Extracts the mag file data as a list of columns specified by 'col_nos' from a text file and formats
    the list according to the format specified in 'fmt'.
    Args:
        file_name       : Text file from which file data has to be extracted
        col_nos         : Indexes of columns to be read from the file ('7:9' or '7,8,9')
        fmt             : Format of storing data in the list
        title_rows      : No. of rows used for title description
    Returns:
        list_filedata   : List of all columns extracted from the text file
    """
    rows, columns = get_file_dmnsn(str(file_name), title_rows=title_rows)

    if re.search(':', col_nos):
        col_nos = range(int(col_nos.split(':')[0]), int(col_nos.split(':')[-1]), 1)
    elif re.search(',', col_nos):
        col_nos = col_nos.split(',')
    else:
        print ("Invalid Format Of Entering Column Numbers")
        sys.exit(1)

    with open(str(file_name), 'r') as f:
        for i in range(0, int(title_rows)):
            f.readline().rstrip()
        data_file = f.read().split()

    list_filedata = []
    for col_index in col_nos:
        list_col = []
        for index in range(0, rows):
            try:
                float(data_file[int(col_index) + index * columns])
            except TypeError:
                fmt = fmt[0:3] + "s}"
                list_col.append(fmt.format(str(data_file[int(col_index) + index * columns])))
            else:
                list_col.append(fmt.format(float(data_file[int(col_index) + index * columns])))
        list_filedata.append(list_col)

    return list_filedata


def read_file(file_name, title_rows=0):
    """
    Extracts the file data as a list of columns from a text file.
    Args:
        file_name       : Text file from which file data has to be extracted
        title_rows      : No. of rows used for title description
    Returns:
        list_filedata   : List of all columns extracted from the text file
    """
    rows, columns = get_file_dmnsn(str(file_name), title_rows=title_rows)

    with open(str(file_name), 'r') as f:
        for i in range(0, int(title_rows)):
            f.readline().rstrip()
        data_file = f.read().split()

    list_filedata = []
    for col_index in range(0, columns):
        list_col = []
        for index in range(0, rows):
            list_col.append(data_file[col_index + index * columns])
        list_filedata.append(list_col)

    return list_filedata


def tabular_mag(file_name, output_file):
    """
    Takes the output from a MAG or ALS file and computes a tabular magnitude file.
    Args:
        file_name   : Input MAG or ALS file
        output_file : Name of the output file onto which the tabular magnitudes are to be written
    Returns:
        None
    """
    file_data = read_file(str(file_name))
    [star_id, img_name, ifilter, xcenter, ycenter, sky_counts, airmass] = file_data[0:7]

    col_data = 7
    star_count = int(max(star_id))
    columns = len(file_data)
    rows = len(file_data[0])
    apertures = (columns - col_data) / 3

    list_col = []
    for index in range(0, 3):
        list_col.append(str(col_data + index * apertures) + ':' + str(col_data + (index + 1) * apertures))

    list_aper = zip_list(read_mag(str(file_name), col_nos=list_col[0], fmt="{:8.2f}"))
    list_mag = zip_list(read_mag(str(file_name), col_nos=list_col[1], fmt="{:10.4f}"))
    list_err = zip_list(read_mag(str(file_name), col_nos=list_col[2], fmt="{:8.3f}"))

    aper_names = ""
    mag_names = ""
    err_names = ""
    for value in range(1, apertures + 1):
        aper_names += "{:8s}".format("AP_" + str(value))
        mag_names += "{:10s}".format("MAG_" + str(value))
        err_names += "{:8s}".format("ERR_" + str(value))

    with open(str(output_file), 'w') as fout:
        fout.write("{0:>4s}{1:>10s}{2:>11s}{3:>13s}  {4}{5}{6}\n\n".format
                   ("ID", "XCENTER", "YCENTER", "SKY_COUNTS", aper_names, mag_names, err_names))

        for index in range(0, rows):
            if index % star_count == 0:
                fout.write("# IMAGE #{0} - {1}, FILTER - {2}, AIRMASS - {3}\n".format
                           ((index / star_count) + 1, img_name[index], ifilter[index], airmass[index]))

            fout.write("{0:>3.0f}{1:>11.3f}{2:>11.3f}{3:>11.4f}{4}{5} {6}\n".format
                       (float(star_id[index]), float(xcenter[index]), float(ycenter[index]), float(sky_counts[index]),
                        list_aper[index], list_mag[index], list_err[index]))

            if index % star_count == star_count - 1:
                fout.write("\n")


