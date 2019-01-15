##############################################################################################################
#   To Fit The Curve Of Instrumental Color Magnitudes Vs Actual Color Magnitudes Plot (Obtain 'm' and 'c')   #
##############################################################################################################


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import os
import re
import math
import glob
import numpy as np
from pyraf import iraf
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt
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
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The 'output_instr_mag2' File Generated After running 'zero_mag.py'
# ------------------------------------------------------------------------------------------------------------------- #
with open('output_instr_nov17') as f:  # Nov17
    data_string = f.read()
    data_list = data_string.split()

with open('output_instr_nov17') as f:  # Nov17
    line_string = f.readline()

with open("stars.coo") as f:
    data_stars = f.read().split()

columns = len(line_string.split())
length_data = len(data_list)
rows = length_data / columns
no_of_stars = len(data_stars) / 2
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Initialise Lists
# ------------------------------------------------------------------------------------------------------------------- #
star_id = []
list_u = []
list_b = []
list_v = []
list_r = []
list_i = []

list_u_err = []
list_b_err = []
list_v_err = []
list_r_err = []
list_i_err = []

list_b_v = []
list_u_b = []
list_v_r = []
list_r_i = []
list_v_i = []

list_b_v_err = []
list_u_b_err = []
list_v_r_err = []
list_r_i_err = []
list_v_i_err = []
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Put Respective Band Magnitudes In A List
# ------------------------------------------------------------------------------------------------------------------- #
remove_index = []
for i in range(1, rows):
    if str(data_list[7 + i * columns]) == 'INDEF' or str(data_list[8 + i * columns]) == 'INDEF':
        remove_index.append(data_list[0 + i * columns])

new_remove_index = list(set(remove_index))
new_remove_index.sort()

usable_stars_nov17 = no_of_stars - len(new_remove_index)

for i in range(1, rows):

    if not new_remove_index.__contains__(data_list[0 + i * columns]):
        if data_list[2 + i * columns] == '7BesU':
            list_u.append(float(data_list[11 + i * columns]))
            list_u_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '6BesB':
            star_id.append(int(data_list[0 + i * columns]))
            list_b.append(float(data_list[11 + i * columns]))
            list_b_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '5BesV':
            list_v.append(float(data_list[11 + i * columns]))
            list_v_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '4BesR':
            list_r.append(float(data_list[11 + i * columns]))
            list_r_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '3BesI':
            list_i.append(float(data_list[11 + i * columns]))
            list_i_err.append(float(data_list[12 + i * columns]))

count_nov17 = len(list_i)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Color Magnitudes From Individual Band Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
for i in range(0, count_nov17):
    list_b_v.append(float("%6.4f" % (list_b[i] - list_v[i])))

    if len(list_u) != 0:
        list_u_b.append(float("%6.4f" % (list_u[i] - list_b[i])))
        list_u_b_err.append(float("%6.4f" % (list_u_err[i] + list_b_err[i])))
    else:
        list_u_b.append(0)
        list_u_b_err.append(0)

    list_v_r.append(float("%6.4f" % (list_v[i] - list_r[i])))
    list_r_i.append(float("%6.4f" % (list_r[i] - list_i[i])))
    list_v_i.append(float("%6.4f" % (list_v[i] - list_i[i])))

    list_b_v_err.append(float("%6.4f" % (list_b_err[i] + list_v_err[i])))
    list_v_r_err.append(float("%6.4f" % (list_v_err[i] + list_r_err[i])))
    list_r_i_err.append(float("%6.4f" % (list_r_err[i] + list_i_err[i])))
    list_v_i_err.append(float("%6.4f" % (list_v_err[i] + list_i_err[i])))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Determining The No. Of Standard Stars In The Field And The No. Of Days Of Observation Of The Standard Field
# ------------------------------------------------------------------------------------------------------------------- #
star_count = int(star_id[-1])
no_of_days = len(star_id) / star_count
object_whole_nov17 = []
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# List Of 'star_count' Elements With 12 Elements Each
# ------------------------------------------------------------------------------------------------------------------- #
for i in range(0, usable_stars_nov17):
    object_whole_nov17.append([list_v[i],
                              list_b_v[i],
                              list_u_b[i],
                              list_v_r[i],
                              list_r_i[i],
                              list_v_i[i],
                              list_v_err[i],
                              list_b_v_err[i],
                              list_u_b_err[i],
                              list_v_r_err[i],
                              list_r_i_err[i],
                              list_v_i_err[i]])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# List Of 2 Elements, Each With 'star_counts' Lists Of 12 Elements Each
# ------------------------------------------------------------------------------------------------------------------- #
"""
for i in range(1, no_of_days):
    object_day = []   # Local Definition Of 'standard_day' is important
    for j in range(0, star_count):
        object_day.append([list_v[i][j],
                    list_b_v[i][j],
                    list_u_b[i][j],
                    list_v_r[i][j],
                    list_r_i[i][j],
                    list_v_i[i][j],
                    list_v_err[i][j],
                    list_b_v_err[i][j],
                    list_u_b_err[i][j],
                    list_v_r_err[i][j],
                    list_r_i_err[i][j],
                    list_v_i_err[i][j]])

    object_whole_nov17.append(object_day)
"""
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# V-Band Magnitude and Color Terms For Landolt Standards
# ------------------------------------------------------------------------------------------------------------------- #

# Standard - PG0918+029
# Name	  V       B-V     U-B     V-R     R-I     V-I      V_Err   B-V_Err  U-B_Err  V-R_Err  R-I_Err  V-I_Err
# A 	14.490   0.536	-0.032   0.325   0.336   0.661	  0.0033   0.0058   0.0095   0.0039   0.0076   0.0085
# B		13.963	 0.765   0.366   0.417   0.370   0.787	  0.0034   0.0072   0.0159   0.0025   0.0045   0.0056
# C		13.537	 0.631   0.087   0.367   0.357   0.722	  0.0020   0.0028   0.0048   0.0015   0.0022   0.0028
# D		12.272   1.044   0.821   0.575   0.535   1.108	  0.0021   0.0030   0.0071   0.0016   0.0018   0.0018
# E 	13.327  -0.271  -1.081  -0.129  -0.159  -0.288    0.0024   0.0024   0.0030   0.0019   0.0055   0.0063


# Standard - PG0231+051
# Name	  V       B-V     U-B     V-R     R-I     V-I     V_Err    B-V_Err  U-B_Err  V-R_Err  R-I_Err  V-I_Err
# A 	12.772   0.710   0.270   0.405   0.394   0.799    0.0008   0.0015   0.0030   0.0011   0.0030   0.0030
# B     14.735   1.448   1.342   0.954   0.998   1.951    0.0030   0.0072   0.0178   0.0034   0.0026   0.0057
# C     13.702   0.671   0.114   0.399   0.385   0.783    0.0014   0.0078   0.0148   0.0028   0.0064   0.0085
# D     14.027   1.088   1.046   0.675   0.586   1.256    0.0029   0.0075   0.0312   0.0081   0.0064   0.0110
# E     13.804   0.677   0.201   0.390   0.369   0.757    0.0046   0.0040   0.0075   0.0035   0.0017   0.0023
# F     16.105  -0.329  -1.192  -0.162  -0.371  -0.534    0.0068   0.0083   0.0045   0.0276   0.1066   0.1221


# Standard - PG0942+029
# Name	  V       B-V     U-B     V-R     R-I     V-I     V_Err    B-V_Err  U-B_Err  V-R_Err  R-I_Err  V-I_Err
# A 	14.731   0.783 	 0.339   0.610   0.477   1.081 	 0.0025    0.0028   0.0075   0.0039   0.0022,  0.0042
# B		14.108   0.525 	 0.085 	 0.368 	 0.333   0.697   0.0025    0.0028   0.0075   0.0039   0.0022,  0.0042
# C		14.989   0.727   0.369   0.539   0.376   0.909   0.0025    0.0028   0.0075   0.0039   0.0022,  0.0042
# D     13.707   0.564   0.129   0.348   0.343   0.687   0.0025    0.0028   0.0075   0.0039   0.0022,  0.0042
# E 	14.004  -0.294  -1.175  -0.130  -0.149  -0.280   0.0045    0.0056   0.0069   0.0069   0.0120   0.0144


PG0918 = [[14.490,  0.536, -0.032,  0.325,  0.336,  0.661,  0.0033,  0.0058,  0.0095,  0.0039,  0.0076,  0.0085],
          [13.963,  0.765,  0.366,  0.417,  0.370,  0.787,  0.0034,  0.0072,  0.0159,  0.0025,  0.0045,  0.0056],
          [13.537,  0.631,  0.087,  0.367,  0.357,  0.722,  0.0020,  0.0028,  0.0048,  0.0015,  0.0022,  0.0028],
          [12.272,  1.044,  0.821,  0.575,  0.535,  1.108,  0.0021,  0.0030,  0.0071,  0.0016,  0.0018,  0.0018],
          [13.327, -0.271, -1.081, -0.129, -0.159, -0.288,  0.0024,  0.0024,  0.0030,  0.0019,  0.0055,  0.0063]]

PG0231 = [[12.772,  0.710,  0.270,  0.405,  0.394,  0.799,  0.0008,  0.0015,  0.0030,  0.0011,  0.0030,  0.0030],
          [14.735,  1.448,  1.342,  0.954,  0.998,  1.951,  0.0030,  0.0072,  0.0178,  0.0034,  0.0026,  0.0057],
          [13.702,  0.671,  0.114,  0.399,  0.385,  0.783,  0.0014,  0.0078,  0.0148,  0.0028,  0.0064,  0.0085],
          [14.027,  1.088,  1.046,  0.675,  0.586,  1.256,  0.0029,  0.0075,  0.0312,  0.0081,  0.0064,  0.0110],
          [13.804,  0.677,  0.201,  0.390,  0.369,  0.757,  0.0046,  0.0040,  0.0075,  0.0035,  0.0017,  0.0023],
          [16.105, -0.329, -1.192, -0.162, -0.371, -0.534,  0.0068,  0.0083,  0.0045,  0.0276,  0.1066,  0.1221]]

PG0942 = [[14.731,  0.783,	0.339,  0.610,	0.477,  1.081,	0.0025,  0.0028,  0.0075,  0.0039,  0.0022,  0.0042],
          [14.108,	0.525,	0.085,	0.368,	0.333,  0.697,  0.0025,  0.0028,  0.0075,  0.0039,  0.0022,  0.0042],
          [14.989,  0.727,  0.369,  0.539,  0.376,  0.909,  0.0025,  0.0028,  0.0075,  0.0039,  0.0022,  0.0042],
          [13.707,  0.564,  0.129,  0.348,  0.343,  0.687,  0.0025,  0.0028,  0.0075,  0.0039,  0.0022,  0.0042],
          [14.004, -0.294, -1.175, -0.130, -0.149, -0.280,  0.0045,  0.0056,  0.0069,  0.0069,  0.0120,  0.0144]]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# List Of Slopes 'alpha' And Intercepts 'beta' For 7 Different Plots
# ------------------------------------------------------------------------------------------------------------------- #

# PG0918
list_alpha0 = [[0.91884, 1.24987, 1.02481, 1.02838, 1.02564, 0.07234, 0.13602],
               [0.92806, 1.24361, 1.00528, 1.06763, 1.03397, 0.05589, 0.10163]]
list_beta0 = [[-0.43371, -2.75126, -0.11184, 0.3059, 0.19287, -0.96033, -0.9635],
              [-0.43142, -2.7694, -0.09413, 0.29826, 0.19144, -0.98643, -0.98781]]

# PG0918-[r2,i2]
list_alpha1 = [[0.91884, 1.24987, 1.02481, 1.02838, 1.02564, 0.07234, 0.13602],
               [0.92806, 1.24361, 1.00984, 1.04614, 1.02741, 0.05589, 0.10163],
               [0.92975, 1.25123, 1.07216, 0.98767, 1.03152, 0.08915, 0.16598]]
list_beta1 = [[-0.43371, -2.75126, -0.11184, 0.3059, 0.19287, -0.96033, -0.9635],
              [-0.43142, -2.7694, -0.0994, 0.29052, 0.18313, -0.98643, -0.98781],
              [-0.54569, -2.80458, -0.19206, 0.29728, 0.12391, -1.10283, -1.10622]]

# PG0231
list_alpha2 = [[0.90692, 1.24877, 1.05414, 1.22637, 1.14145, 0.04808, 0.08085],
               [0.90735, 1.25413, 1.05396, 1.23278, 1.14428, 0.04805, 0.08278],
               [0.90502, 1.24325, 1.06287, 1.24341, 1.15353, 0.05074, 0.08818]]
list_beta2 = [[-0.4207, -2.70666, -0.10869, 0.26739, 0.12003, -0.90718, -0.90886],
              [-0.41851, -2.73126, -0.11633, 0.24002, 0.08507, -0.95799, -0.96055],
              [-0.36493, -2.71224, -0.08315, 0.23061, 0.11265, -0.7919, -0.79494]]

# [[0.89561, 1.23645, 1.03478, 1.22206, 1.11949, 0.05322, 0.08953],
# [0.89596, 1.24176, 1.03444, 1.22849, 1.12223, 0.05313, 0.0918],
# [0.89372, 1.23084, 1.04351, 1.23915, 1.13136, 0.05586, 0.09751]]
# [[-0.41727, -2.7006, -0.09839, 0.265, 0.13378, -0.91112, -0.91306],
# [-0.41593, -2.72249, -0.10579, 0.23772, 0.09953, -0.96189, -0.96489],
# [-0.36217, -2.70454, -0.07339, 0.22834, 0.12653, -0.79583, -0.79943]]


# PG0942
list_alpha3 = [[0.83646, 1.09091, 1.1164, 0.9655, 1.03599, 0.0809, 0.12108]]
list_beta3 = [[-0.40574, -2.45694, -0.09777, 0.25939, 0.18824, -0.9687, -0.97341]]

# PG0918 + PG0231
list_alpha4 = [[0.91101, 1.25445, 1.05072, 1.1574, 1.10423, 0.06301, 0.10936],
               [0.91285, 1.25775, 1.03876, 1.16652, 1.1018, 0.05472, 0.09559]]
list_beta4 = [[-0.42562, -2.73903, -0.11401, 0.28979, 0.15175, -0.93483, -0.93682],
              [-0.42087, -2.76785, -0.10797, 0.27143, 0.13649, -0.97322, -0.97518]]

# shubham's - PG0918
list_alpha6 = [[0.92385, 1.25277, 0.98378, 1.0749, 1.02644, 0.05842, 0.10818],
               [0.92873, 1.2344, 1.0147, 1.04385, 1.0274, 0.05482, 0.10031],
               [0.92091, 1.2482, 1.01834, 1.04084, 1.02822, 0.07838, 0.1458]]
list_beta6 = [[-0.45232, -2.82862, -0.10042, 0.26114, 0.14329, -0.95581, -0.95784],
              [-0.45089, -2.85013, -0.12796, 0.27381, 0.13978, -0.98526, -0.9868],
              [-0.47421, -2.82842, -0.09878, 0.26115, 0.15792, -0.8984, -0.90134]]

# shubham's - PG0231
list_alpha7 = [[0.90738, 1.24693, 1.0514, 1.23039, 1.14191, 0.04784, 0.08058],
               [0.91079, 1.23985, 1.05333, 1.23489, 1.14471, 0.04775, 0.08223]]
list_beta7 = [[-0.42369, -2.76708, -0.10077, 0.2276, 0.09023, -0.84134, -0.84307],
              [-0.42841, -2.76269, -0.12705, 0.23667, 0.0699, -0.89444, -0.89697]]


# Standard alpha - PG0231
list_alpha8 = [[0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]] * 3
list_beta8 = [[-0.43904, -2.73217, -0.08239, 0.32753, 0.24899, -0.89123, -0.89401],
              [-0.43563, -2.74358, -0.09072, 0.29861, 0.2179, -0.94617, -0.94824],
              [-0.37687, -2.73693, -0.05211, 0.29379, 0.25121, -0.76848, -0.77126]]
list_beta8_err = [[0.00784, 0.01790, 0.00966, 0.02248, 0.02654, 0.00417, 0.00421],
                  [0.00756, 0.01715, 0.00959, 0.02240, 0.02650, 0.00417, 0.00421],
                  [0.00746, 0.01641, 0.00952, 0.02267, 0.02657, 0.00406, 0.00411]]

# Standard alpha - PG0918
list_alpha9 = [[0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]] * 3
list_beta9 = [[-0.43854, -2.82398, -0.11042, 0.3078, 0.20301, -0.93788, -0.93649],
              [-0.41438, -2.84075, -0.09585, 0.2956, 0.20142, -0.97317, -0.97247],
              [-0.54889, -2.86547, -0.17562, 0.2958, 0.14422, -1.06888, -1.06749]]
list_beta9_err = [[0.00564, 0.01121, 0.00468, 0.00613, 0.00661, 0.00334, 0.00334],
                  [0.00555, 0.01036, 0.00452, 0.00601, 0.00661, 0.00334, 0.00334],
                  [0.00555, 0.01054, 0.00486, 0.00689, 0.00716, 0.00334, 0.00334]]

# PG0231
list_alpha10 = [[0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]] * 3
list_beta10 = [[-0.4276, -2.8318, -0.0668, 0.2819, 0.2249, -0.8375, -0.8393],
               [-0.4284, -2.8314, -0.0959, 0.2939, 0.2067, -0.8907, -0.8924]]

# PG0918
list_alpha11 = [[0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]] * 3
list_beta11 = [[-0.4238, -2.8625, -0.0979, 0.2820, 0.1860, -0.8851, -0.8844],
               [-0.4099, -2.8890, -0.1060, 0.2800, 0.1760, -0.9193, -0.9186]]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots To Determine Coefficients
# ------------------------------------------------------------------------------------------------------------------- #

# 	PLOT		    		X:Y		 	X_ERR:Y_ERR
# (B-V)obs Vs (B-V)		[x,1]:[y,1]	 	[x,7]:[y,7]
# (U-B)obs Vs (U-B)		[x,2]:[y,2]	 	[x,8]:[y,8]
# (V-R)obs Vs (V-R)		[x,3]:[y,3]	 	[x,9]:[y,9]
# (R-I)obs Vs (R-I)		[x,4]:[y,4]	 	[x,10]:[y,10]
# (V-I)obs Vs (V-I)		[x,5]:[y,5]	 	[x,11]:[y,11]
# (B-V) Vs V-Vobs   [y,1]:[y,0]-[x,0]  [y,7]:[y,6]-[x,6]
# (V-R) Vs V-Vobs   [y,3]:[y,0]-[x,0]  [y,9]:[y,6]-[x,6]
# X = standard_whole
# Y = standard(standard_field)[PG0918, PG0231, PG2213 etc.]


list_x = [1, 2, 3, 4, 5, 1, 3]
list_x_err = [7, 8, 9, 10, 11, 7, 9]
list_y = [1, 2, 3, 4, 5, 0, 0]
list_y_err = [7, 8, 9, 10, 11, 6, 6]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Obtaining Secondary-Standard Color Terms And V-Band Magnitudes And Appending To A List 'object_standard'
# ------------------------------------------------------------------------------------------------------------------- #

def object_std(object_whole, list_alpha, list_beta, list_beta_err, day, usable_stars):
    object_stand = []
    object_stand_err = []
    for j in range(0, usable_stars):
        x = []
        x_err = []
        for k in range(0, len(list_y)):
            if list_y[k] != 0 and object_whole[j][list_x[k]] != 0:
                x.append(float("%6.4f" % (list_alpha[day][k] * object_whole[j][list_x[k]] + list_beta[day][k])))
                x_err.append(float("%6.4f" % (list_alpha[day][k] * object_whole[j][list_x[k] + 6] +
                                              list_beta_err[day][k])))
            elif object_whole[j][list_x[k]] == 0:
                x.append(0)
                x_err.append(0)
            else:
                x.append(float("%6.4f" % (object_whole[j][list_y[k]] + list_alpha[day][k] *
                                          (list_alpha[day][int(list_x[k] - 1)] * object_whole[j][list_x[k]] +
                                          list_beta[day][int(list_x[k]-1)]) + list_beta[day][k])))
                x_err.append(float("%6.4f" % (object_whole[j][list_y[k] + 6] + list_alpha[day][k] *
                                              (list_alpha[day][int(list_x[k] - 1)] * object_whole[j][list_x[k] + 6] +
                                              list_beta_err[day][int(list_x[k]-1)]) + list_beta_err[day][k])))

        object_stand.append(x)
        object_stand_err.append(x_err)
    return object_stand, object_stand_err


object_0918_nov17, object_err_0918_nov17 = object_std(object_whole_nov17, list_alpha9, list_beta9,
                                                      list_beta9_err, 0, usable_stars_nov17)
object_0231_nov17, object_err_0231_nov17 = object_std(object_whole_nov17, list_alpha8, list_beta8,
                                                      list_beta8_err, 0, usable_stars_nov17)

object_shubham_0918_nov17, object_shubham_err_0918_nov17 = object_std(object_whole_nov17, list_alpha11, list_beta11,
                                                                      list_beta9_err, 0, usable_stars_nov17)
object_shubham_0231_nov17, object_shubham_err_0231_nov17 = object_std(object_whole_nov17, list_alpha10, list_beta10,
                                                                      list_beta9_err, 0, usable_stars_nov17)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Obtaining Secondary-Standard U, B, V, R & I Magnitudes And Appending To A List 'mag_standard'
# ------------------------------------------------------------------------------------------------------------------- #

def mag_std1(object_standard, v):
    mag_stand = []
    choice_v = v
    for j in range(0, usable_stars_nov17):
        if object_standard == object_0918_nov17 or object_standard == object_0231_nov17 or \
                        object_standard == object_shubham_0918_nov17 or object_standard == object_shubham_0231_nov17:
            u_mag = 'INDEF'
        else:
            u_mag = float("%6.4f" % (object_standard[j][1] + object_standard[j][0] + object_standard[j][choice_v]))
        b_mag = float("%6.4f" % (object_standard[j][0] + object_standard[j][choice_v]))
        v_mag = float("%6.4f" % (object_standard[j][choice_v]))
        r_mag = float("%6.4f" % (object_standard[j][choice_v] - object_standard[j][2]))
        i_mag = float("%6.4f" % (object_standard[j][choice_v] - object_standard[j][4]))
        mag_stand.append([u_mag, b_mag, v_mag, r_mag, i_mag])
    return mag_stand


mag_0918_nov17 = mag_std1(object_0918_nov17, 5)
mag_0231_nov17 = mag_std1(object_0231_nov17, 5)

mag_shubham_0918_nov17 = mag_std1(object_shubham_0918_nov17, 5)
mag_shubham_0231_nov17 = mag_std1(object_shubham_0231_nov17, 5)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read The 'output_instr_mag2' File Generated After running 'zero_mag.py'
# ------------------------------------------------------------------------------------------------------------------- #
with open('output_instr_nov23') as f:  # Nov23
    data_string = f.read()
    data_list = data_string.split()

with open('output_instr_nov23') as f:  # Nov23
    line_string = f.readline()

columns = len(line_string.split())
length_data = len(data_list)
rows = length_data / columns
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Initialize Lists To Be Printed In The Output File
# ------------------------------------------------------------------------------------------------------------------- #
star_id = []
list_u = []
list_b = []
list_v = []
list_r = []
list_i = []

list_u_err = []
list_b_err = []
list_v_err = []
list_r_err = []
list_i_err = []

list_b_v = []
list_u_b = []
list_v_r = []
list_r_i = []
list_v_i = []

list_b_v_err = []
list_u_b_err = []
list_v_r_err = []
list_r_i_err = []
list_v_i_err = []

remove_index = []
for i in range(1, rows):
    if str(data_list[7 + i * columns]) == 'INDEF' or str(data_list[8 + i * columns]) == 'INDEF':
        remove_index.append(data_list[0 + i * columns])


new_remove_index = list(set(remove_index))
new_remove_index.sort()
usable_stars_nov23 = no_of_stars - len(new_remove_index)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Put Respective Band Magnitudes In A List
# ------------------------------------------------------------------------------------------------------------------- #

for i in range(1, rows):
    if not new_remove_index.__contains__(data_list[0 + i * columns]):
        if data_list[2 + i * columns] == '7BesU':
            list_u.append(float(data_list[11 + i * columns]))
            list_u_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '6BesB':
            star_id.append(int(data_list[0 + i * columns]))
            list_b.append(float(data_list[11 + i * columns]))
            list_b_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '5BesV':
            list_v.append(float(data_list[11 + i * columns]))
            list_v_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '4BesR':
            list_r.append(float(data_list[11 + i * columns]))
            list_r_err.append(float(data_list[12 + i * columns]))

        elif data_list[2 + i * columns] == '3BesI':
                list_i.append(float(data_list[11 + i * columns]))
                list_i_err.append(float(data_list[12 + i * columns]))

count_nov23 = len(list_i)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Color Magnitudes From Individual Band Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

for i in range(0, count_nov23):
    list_b_v.append(float("%6.4f" % (list_b[i] - list_v[i])))
    list_u_b.append(float("%6.4f" % (list_u[i] - list_b[i])))
    list_v_r.append(float("%6.4f" % (list_v[i] - list_r[i])))
    list_r_i.append(float("%6.4f" % (list_r[i] - list_i[i])))
    list_v_i.append(float("%6.4f" % (list_v[i] - list_i[i])))

    list_b_v_err.append(float("%6.4f" % (list_b_err[i] + list_v_err[i])))
    list_u_b_err.append(float("%6.4f" % (list_u_err[i] + list_b_err[i])))
    list_v_r_err.append(float("%6.4f" % (list_v_err[i] + list_r_err[i])))
    list_r_i_err.append(float("%6.4f" % (list_r_err[i] + list_i_err[i])))
    list_v_i_err.append(float("%6.4f" % (list_v_err[i] + list_i_err[i])))

# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Color Magnitudes From Individual Band Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# List Of 'star_count' Elements With 12 Elements Each
# ------------------------------------------------------------------------------------------------------------------- #
object_whole_nov23 = []
for i in range(0, usable_stars_nov23):
    object_whole_nov23.append([list_v[i], list_b_v[i], list_u_b[i], list_v_r[i], list_r_i[i], list_v_i[i],
                               list_v_err[i], list_b_v_err[i], list_u_b_err[i], list_v_r_err[i],
                               list_r_i_err[i], list_v_i_err[i]])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Obtaining Secondary-Standard U, B, V, R & I Magnitudes And Appending To A List 'mag_standard'
# ------------------------------------------------------------------------------------------------------------------- #

def mag_std2(object_standard, v):
    mag_stand = []
    choice_v = v
    for j in range(0, usable_stars_nov23):
        u_mag = float("%6.4f" % (object_standard[j][0] + object_standard[j][1] + object_standard[j][choice_v]))
        b_mag = float("%6.4f" % (object_standard[j][0] + object_standard[j][choice_v]))
        v_mag = float("%6.4f" % (object_standard[j][choice_v]))
        r_mag = float("%6.4f" % (object_standard[j][choice_v] - object_standard[j][2]))
        i_mag = float("%6.4f" % (object_standard[j][choice_v] - object_standard[j][4]))
        mag_stand.append([u_mag, b_mag, v_mag, r_mag, i_mag])
    return mag_stand

object_0918_nov23, object_err_0918_nov23 = object_std(object_whole_nov23, list_alpha9, list_beta9,
                                                      list_beta9_err, 1, usable_stars_nov23)
object_0231_nov23, object_err_0231_nov23 = object_std(object_whole_nov23, list_alpha8, list_beta8,
                                                      list_beta8_err, 1, usable_stars_nov23)

object_shubham_0918_nov23, object_shubham_err_0918_nov23 = object_std(object_whole_nov23, list_alpha11, list_beta11,
                                                                      list_beta9_err, 0, usable_stars_nov23)
object_shubham_0231_nov23, object_shubham_err_0231_nov23 = object_std(object_whole_nov23, list_alpha10, list_beta10,
                                                                      list_beta8_err, 0, usable_stars_nov23)

mag_0918_nov23 = mag_std2(object_0918_nov23, 5)
mag_0231_nov23 = mag_std2(object_0231_nov23, 5)

mag_shubham_0918_nov23 = mag_std2(object_shubham_0918_nov23, 5)
mag_shubham_0231_nov23 = mag_std2(object_shubham_0231_nov23, 5)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Mean And Standard Deviation Of Star Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #
plus_minus = "+/-"


def mean_stdev(mag_list):
    mag_value = []
    for i in range(0, usable_stars_nov23):
        mag_u = []
        mag_b = []
        mag_v = []
        mag_r = []
        mag_i = []
        for j in range(0, len(mag_list)):
            if type(mag_list[j][i][0]) == float:
                mag_u.append(mag_list[j][i][0])
            mag_b.append(mag_list[j][i][1])
            mag_v.append(mag_list[j][i][2])
            mag_r.append(mag_list[j][i][3])
            mag_i.append(mag_list[j][i][4])
        if len(mag_u) != 0:
            mag_value.append([str("%6.4f" % np.mean(mag_u)) + plus_minus + str("%6.4f" % np.std(mag_u)),
                              str("%6.4f" % np.mean(mag_b)) + plus_minus + str("%6.4f" % np.std(mag_b)),
                              str("%6.4f" % np.mean(mag_v)) + plus_minus + str("%6.4f" % np.std(mag_v)),
                              str("%6.4f" % np.mean(mag_r)) + plus_minus + str("%6.4f" % np.std(mag_r)),
                              str("%6.4f" % np.mean(mag_i)) + plus_minus + str("%6.4f" % np.std(mag_i))])
        else:
            mag_value.append([str("%16s" % "INDEF      "),
                              str("%6.4f" % np.mean(mag_b)) + plus_minus + str("%6.4f" % np.std(mag_b)),
                              str("%6.4f" % np.mean(mag_v)) + plus_minus + str("%6.4f" % np.std(mag_v)),
                              str("%6.4f" % np.mean(mag_r)) + plus_minus + str("%6.4f" % np.std(mag_r)),
                              str("%6.4f" % np.mean(mag_i)) + plus_minus + str("%6.4f" % np.std(mag_i))])
    return mag_value

mag_list1 = [mag_0231_nov17, mag_0231_nov23, mag_0918_nov17, mag_0918_nov23]
mag_list2 = [mag_0231_nov23, mag_0918_nov23]
mag_list3 = [mag_0231_nov17, mag_0918_nov17]
mag_list4 = [mag_0231_nov17, mag_0231_nov23, mag_0918_nov23]

mag_shubham_list1 = [mag_shubham_0231_nov17, mag_shubham_0231_nov23, mag_shubham_0918_nov17, mag_shubham_0918_nov23]
mag_shubham_list2 = [mag_shubham_0231_nov23, mag_shubham_0918_nov23]
mag_shubham_list3 = [mag_shubham_0231_nov17, mag_shubham_0918_nov17]
mag_shubham_list4 = [mag_shubham_0231_nov17, mag_shubham_0231_nov23, mag_shubham_0918_nov23]

m1 = mean_stdev(mag_list1)
m2 = mean_stdev(mag_list2)
m3 = mean_stdev(mag_list3)
m4 = mean_stdev(mag_list4)

m5 = mean_stdev(mag_shubham_list1)
m6 = mean_stdev(mag_shubham_list2)
m7 = mean_stdev(mag_shubham_list3)
m8 = mean_stdev(mag_shubham_list4)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print The Secondary Standard Magnitudes In The Output File 'output_nov1723'
# ------------------------------------------------------------------------------------------------------------------- #

f = open("output_nov1723", "w")
f.write("Secondary Standards Mean Magnitudes" + "\n" + "\n")
f.write(str('STAR_ID') + "    " + str("%8s" % 'MAG_U') + "                 " + str("%8s" % 'MAG_B') + "             " +
        str("%8s" % 'MAG_V') + "               " + str("%8s" % 'MAG_R') + "            " + str("%8s" % 'MAG_I') +
        "\n" + "\n" + "Mean +/- Standard Deviation" + "\n" + "\n")


def write(message, m1, usable_stars):
    f.write(message + "\n" + "\n")
    for j in range(0, usable_stars):
        f.write(str("%3.0f" % star_id[j]) + "        " +
                str(m1[j][0]) + "     " +
                str(m1[j][1]) + "     " +
                str(m1[j][2]) + "     " +
                str(m1[j][3]) + "     " +
                str(m1[j][4]) + "\n")
    f.write("\n")

write("Mean Magnitudes", m1, usable_stars_nov23)
write("Mean Magnitudes - Shubham", m5, usable_stars_nov23)
write("Best Match", m4, usable_stars_nov17)
write("Best Match - Shubham", m8, usable_stars_nov17)
write("Nov23_Mag", m2, usable_stars_nov23)
write("Nov23_Mag - Shubham", m6, usable_stars_nov23)
write("Nov17_Mag", m3, usable_stars_nov17)
write("Nov17_Mag - Shubham", m7, usable_stars_nov17)
f.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print The Secondary Standard Magnitudes In The Output File 'output_mag_nov1723'
# ------------------------------------------------------------------------------------------------------------------- #
remove_file('output_mag_nov1723')
f = open("output_mag_nov1723", "w")
f.write("Secondary Standards Magnitudes" + "\n" + "\n")


def write_data(message, mag_standard, usable_stars):
    f.write(str(message) + "\n" + "\n")
    for j in range(0, usable_stars):
        if type(mag_standard[j][0]) != str and type(mag_standard[j][1]) != str:
            f.write(str("%3.0f" % star_id[j]) + "              " +
                    str("%6.4f" % mag_standard[j][0]) + "              " +
                    str("%6.4f" % mag_standard[j][1]) + "              " +
                    str("%6.4f" % mag_standard[j][2]) + "               " +
                    str("%6.4f" % mag_standard[j][3]) + "            " +
                    str("%6.4f" % mag_standard[j][4]) + "\n")
        else:
            f.write(str("%3.0f" % star_id[j]) + "              " +
                    str("%8s" % mag_standard[j][0]) + "              " +
                    str("%8s" % mag_standard[j][1]) + "              " +
                    str("%6.4f" % mag_standard[j][2]) + "               " +
                    str("%6.4f" % mag_standard[j][3]) + "            " +
                    str("%6.4f" % mag_standard[j][4]) + "\n")
    f.write("\n")


write_data("PG0918-nov23", mag_0918_nov23, usable_stars_nov23)
write_data("PG0918-nov23 - Shubham", mag_shubham_0918_nov23, usable_stars_nov23)

write_data("PG0918-nov17", mag_0918_nov17, usable_stars_nov17)
write_data("PG0918-nov17 - Shubham", mag_shubham_0918_nov17, usable_stars_nov17)

write_data("PG0231-nov23", mag_0231_nov23, usable_stars_nov23)
write_data("PG0231-nov23 - Shubham", mag_shubham_0231_nov23, usable_stars_nov23)

write_data("PG0231-nov17", mag_0231_nov17, usable_stars_nov17)
write_data("PG0231-nov17 - Shubham", mag_shubham_0231_nov17, usable_stars_nov17)
f.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print The Color Terms And Their Errors The Output File 'output_color_nov1723'
# ------------------------------------------------------------------------------------------------------------------- #
remove_file('output_color_nov1723')
f = open("output_color_nov1723", "w")
f.write("Secondary Standards Color Terms With Errors" + "\n" + "\n")


def write_color_mag(message, object_color, object_color_err, usable_stars):
    f.write(str(message) + "\n" + "\n")
    f.write(str('STAR_ID') + " " * 5 + str('B-V') + " " * 5 + str('B-V Err') + " " * 5 + str('U-B') + " " * 5 +
            str('U-B Err') + " " * 5 + str('V-R') + " " * 5 + str('V-R Err') + " " * 5 + str('R-I') + " " * 5 +
            str('R-I Err') + " " * 5 + str('V-I') + " " * 5 + str('V-I Err') + " " * 5 + str('V-1') + " " * 5 +
            str('V-1 Err') + " " * 5 + str('V-2') + " " * 5 + str('V-2 Err') + "\n" + "\n")

    for i in range(usable_stars):
        f.write(str("%3.0f" % star_id[i]) + " " * 5 + str("%8s" % object_color[i][0]) + " " * 2 +
                str("%8s" % object_color_err[i][0]) + " " * 2 + str("%8s" % object_color[i][1]) + " " * 2 +
                str("%8s" % object_color_err[i][1]) + " " * 2 + str("%8s" % object_color[i][2]) + " " * 2 +
                str("%8s" % object_color_err[i][2]) + " " * 2 + str("%8s" % object_color[i][3]) + " " * 2 +
                str("%8s" % object_color_err[i][3]) + " " * 2 + str("%8s" % object_color[i][4]) + " " * 2 +
                str("%8s" % object_color_err[i][4]) + " " * 2 + str("%8s" % object_color[i][5]) + " " * 2 +
                str("%8s" % object_color_err[i][5]) + " " * 2 + str("%8s" % object_color[i][6]) + " " * 2 +
                str("%8s" % object_color_err[i][6]) + "\n")


write_color_mag("PG0231-nov17", object_0231_nov17, object_err_0231_nov17, usable_stars_nov17)
write_color_mag("PG0231-nov17 - Shubham", object_shubham_0231_nov17, object_shubham_err_0231_nov17, usable_stars_nov17)

write_color_mag("PG0918-nov17", object_0918_nov17, object_err_0918_nov17, usable_stars_nov17)
write_color_mag("PG0918-nov17 - Shubham", object_shubham_0918_nov17, object_shubham_err_0918_nov17, usable_stars_nov17)

write_color_mag("PG0231-nov23", object_0231_nov23, object_err_0231_nov23, usable_stars_nov23)
write_color_mag("PG0231-nov23 - Shubham", object_shubham_0231_nov23, object_shubham_err_0231_nov23, usable_stars_nov23)

write_color_mag("PG0918-nov23", object_0918_nov23, object_err_0918_nov23, usable_stars_nov23)
write_color_mag("PG0918-nov23 - Shubham", object_shubham_0918_nov23, object_shubham_err_0918_nov23, usable_stars_nov23)
f.close()

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Location Of Directory Where PHotometric Analysis Is To Be Performed
# ------------------------------------------------------------------------------------------------------------------- #
DIR_PATH = "/home/avinash/Supernovae_Data/ASASSN14dq/Photometry/"
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Mean Mag Matrix From Individual Magnitude Files
# ------------------------------------------------------------------------------------------------------------------- #
mag_mean = []
for i in range(0, len(mag_0231_nov17)):
    temp_mag = []
    for j in range(0, len(mag_0231_nov17[0])):
        if mag_0231_nov17[i][j] != 'INDEF':
            mean_value = "%6.4f" % np.mean([float(mag_0231_nov17[i][j]), float(mag_0231_nov23[i][j]),
                                            float(mag_0918_nov23[i][j])])
        else:
            mean_value = "%6.4f" % np.mean([float(mag_0231_nov23[i][j]), float(mag_0918_nov23[i][j])])
        temp_mag.append(mean_value)
    mag_mean.append(temp_mag)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print Mean Of Secondary Standard Magnitude In Output File 'output_mean_nov1723'
# ------------------------------------------------------------------------------------------------------------------- #
remove_file('output_mean_nov1723')
with open("output_mean_nov1723", "w") as f:
    f.write("Secondary Standards Mean Magnitudes" + "\n" + "\n")
    for j in range(0, len(mag_mean)):
        f.write(str("%3.0f" % int(j + 1)) + "              " +
                str(mag_mean[j][0]) + "              " +
                str(mag_mean[j][1]) + "              " +
                str(mag_mean[j][2]) + "               " +
                str(mag_mean[j][3]) + "            " +
                str(mag_mean[j][4]) + "\n")
    f.write("\n")
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


def txdump(list_files, output_file, choice):
    """
    Args:
        list_files  : List of all MAG, ALS files from which magnitudes are to be extracted
        output_file : Output file where data from the list of input files is to be written
        choice      : Category of files from which magnitudes are to be extracted
    Returns:
        None
    """
    fields = ""
    if choice == "als":
        fields = "IMAGE, IFILTER, ID, MSKY, XAIRMASS, MAG, MERR"
    elif choice == "mag":
        fields = "ID, IMAGE, XCENTER, YCENTER, IFILTER, MSKY, STDEV, RAPERT, MAG, MERR, XAIRMASS"
    task = iraf.noao.digiphot.ptools.txdump
    task.unlearn()

    task(textfile="@" + str(list_files), fields=fields, expr="yes", Stdout=output_file)


def list_statistics(input_list):
    """
    Returns the statistics of the list of elements in the input list
    Args:
        input_list  : Input list
    Returns:
        value_mean  : Mean of the list of elements
        value_median: Median of the list of elements
        value_std   : Standard Deviation of the list of elements

    """
    value_mean = np.mean(input_list)
    value_median = np.median(input_list)
    value_std = np.std(input_list)

    return value_mean, value_median, value_std


def remove_similar_files(common_text):
    """
    Removes similar files based on the string "common_text"
    Args:
        common_text : String containing partial filename
    Returns:
        None
    """
    os.chdir(DIR_PATH)
    for residual_files in glob.glob(common_text):
        os.remove(residual_files)


def reject(input_list):
    """
    Rejects outliers from the input list
    Args:
        input_list  : Input list
    Returns:
        reject_list : Output list after rejecting outliers from the input list
    """
    reject_list = []
    pop = False
    for index in range(0, len(input_list)):
        reject_list.append(float(input_list[index]))

    reject_list.sort()
    value_mean, value_median, value_std = list_statistics(reject_list)
    if abs(reject_list[0] - value_median) < abs(reject_list[-1] - value_median):
        remove_index = -1
    else:
        remove_index = 0
    if abs(reject_list[remove_index] - value_median) > value_std:
        reject_list.pop(remove_index)
        pop = True

    if pop:
        value_mean, value_median, value_std = list_statistics(reject_list)
        if abs(reject_list[0] - value_median) < abs(reject_list[-1] - value_median):
            remove_index = -1
        else:
            remove_index = 0
        if abs(reject_list[remove_index] - value_median) > 2 * value_std:
            reject_list.pop(remove_index)

    return reject_list

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Makes A List Of FITS Files
# ------------------------------------------------------------------------------------------------------------------- #
os.chdir(DIR_PATH)
remove_similar_files("tmp*.fits")
list_fits = group_similar_files('', common_text="*.fits", exceptions=".fits.")
list_fits.sort()
no_of_fits_files = len(list_fits)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Makes A List Of Dates On Which Observation Was Done
# ------------------------------------------------------------------------------------------------------------------- #
list_dates = []
for index in range(0, no_of_fits_files):
    if not list_fits[index][0:5] in list_dates:
        list_dates.append(list_fits[index][0:5])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Groups Mag Files From PHOT Task Into Separate Lists
# ------------------------------------------------------------------------------------------------------------------- #
list_mag4 = group_similar_files("list_mag4", "*.mag.4")
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Groups Files From The Date Into A File "list_date_mag4", Applies TXDUMP Task To Obtain "output_date_mag4"
# ------------------------------------------------------------------------------------------------------------------- #
for date in list_dates:
    group_similar_files("list_" + str(date) + "_mag4", str(date) + "*.mag.4")
    txdump("list_" + str(date) + "_mag4", "output_" + str(date) + "_mag4", "mag")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Finds Out No. Of Stars Used As Secondary Standards In The Supernova Field
# ------------------------------------------------------------------------------------------------------------------- #
with open("stars.coo") as f:
    data_stars = f.read().split()
no_of_stars = len(data_stars) / 2
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Calculating Instrumental Magnitudes And Color Terms
# ------------------------------------------------------------------------------------------------------------------- #

def get_instr_mag(date):
    with open("output_" + str(date) + "_mag4") as f:
        data_list = f.read().split()

    with open("output_" + str(date) + "_mag4") as f:
        data_line = f.readline().split()

    length_data = len(data_list)
    columns = len(data_line)
    rows = length_data / columns

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

    # FILTER  EXTINCTION_MEAN   EXTINCTION_ERROR
    #   U        0.36              0.07
    #   B        0.21              0.04
    #   V        0.12              0.04
    #   R        0.09              0.04
    #   I        0.05              0.03

    for i in range(0, rows):

        if str(data_list[9 + i * columns]) != 'INDEF' and str(data_list[10 + i * columns]) != 'INDEF':
            mag_phot_1 = "%6.4f" % float(data_list[9 + i * columns])
            mag_phot_2 = "%6.4f" % float(data_list[10 + i * columns])
            mag_err_phot_1 = "%6.4f" % float(data_list[11 + i * columns])
            mag_err_phot_2 = "%6.4f" % float(data_list[12 + i * columns])
            correction_aperture = "%6.4f" % (float(mag_phot_1) - float(mag_phot_2))

        else:
            mag_phot_1 = 'INDEF'
            mag_err_phot_1 = 'INDEF'
            mag_phot_2 = 'INDEF'
            mag_err_phot_2 = 'INDEF'
            correction_aperture = 'INDEF'

        if data_list[13 + i * columns] != 'INDEF':
            airmass_value = "%5.4f" % float(data_list[13 + i * columns])
        else:
            airmass_value = 'INDEF'

        star_id.append(float(data_list[0 + i * columns]))
        image_name.append(str(data_list[1 + i * columns]))
        xcenter.append(float(data_list[2 + i * columns]))
        ycenter.append(float(data_list[3 + i * columns]))
        ifilter.append(str(data_list[4 + i * columns]))
        sky_counts.append("%8.4f" % float(data_list[5 + i * columns]))
        sigma_counts.append(float(data_list[6 + i * columns]))
        aperture_1.append(float(data_list[7 + i * columns]))
        aperture_2.append(float(data_list[8 + i * columns]))
        airmass.append(airmass_value)

        phot_mag_1.append(mag_phot_1)
        phot_mag_err_1.append(mag_err_phot_1)
        phot_mag_2.append(mag_phot_2)
        phot_mag_err_2.append(mag_err_phot_2)
        aperture_correction.append(correction_aperture)

    star_count = int(star_id[-1])

    for value in range(0, rows):
        quotient = value // star_count
        aperture_correction_list = []
        temp_count = 0
        if aperture_correction[value] != 'INDEF':
            for index in range(quotient * star_count, (quotient + 1) * star_count):
                if aperture_correction[index] != 'INDEF':
                    aperture_correction_list.append(float(aperture_correction[index]))
            aperture_correction_mean.append("%6.4f" % np.mean(aperture_correction_list))
        else:
            aperture_correction_mean.append('INDEF')
            temp_count += 1

        if temp_count == star_count:
            aperture_correction_mean.append('INDEF')

    for i in range(0, rows):
        if phot_mag_1[i] != 'INDEF' and aperture_correction_mean[i] != 'INDEF':
            mag_instr = "%6.4f" % (float(phot_mag_1[i]) - float(aperture_correction_mean[i]))
            mag_err_instr = "%6.4f" % (float(phot_mag_err_1[i]))
        else:
            mag_instr = 'INDEF'
            mag_err_instr = 'INDEF'

        instr_mag.append(mag_instr)
        instr_mag_err.append(mag_err_instr)

    with open("output_instr_" + str(date) + "_mag4", "w") as f:
        f.write(str('STAR_ID') + "        " + str('IMAGE_NAME') + "         " + str('FILTER') + "  " +
                str('SKY_COUNTS') + "  " + str('AIRMASS') + "  " + str('APER_1') + "   " + str('APER_2') + "    " +
                str('MAG_1') + "    " + str('MAG_2') + "     " + str('APER_CORR') + "   " + str('AP_CO_MEAN') + "   " +
                str('INSTR_MAG') + "  " + str('INSTR_MAG_ERR') + "\n" + "\n")

        for i in range(0, rows):
            f.write(str("%3.0f" % star_id[i]) + "     " + str("%23s" % image_name[i]) + "   " +
                    str("%6s" % ifilter[i]) + "   " + str("%9s" % float(sky_counts[i])) + "   " +
                    str("%6s" % airmass[i]) + "   " + str("%4.2f" % aperture_1[i]) + "    " +
                    str("%4.2f" % aperture_2[i]) + "   " + str("%7s" % phot_mag_1[i]) + "   " +
                    str("%7s" % phot_mag_2[i]) + "      " + str("%6s" % aperture_correction[i]) + "     " +
                    str("%6s" % aperture_correction_mean[i]) + "      " + str("%6s" % instr_mag[i]) + "      " +
                    str("%6s" % instr_mag_err[i]) + "   " + "\n")

            if star_id[i] == int(star_count):
                f.write("\n")


def calculate_beta_color(file_name, mag_matrix):
    with open(str(file_name), "r") as f:
        data_list = f.read().split()
    with open(str(file_name), "r") as f:
        data_line = f.readline().split()

    columns = len(data_line)
    length_data = len(data_list)
    rows = length_data / columns

    list_u = []
    list_b = []
    list_v = []
    list_r = []
    list_i = []

    list_u_err = []
    list_b_err = []
    list_v_err = []
    list_r_err = []
    list_i_err = []

    list_b_v = []
    list_u_b = []
    list_v_r = []
    list_r_i = []
    list_v_i = []

    list_b_v_err = []
    list_u_b_err = []
    list_v_r_err = []
    list_r_i_err = []
    list_v_i_err = []

    list_filters = []
    for index in range(1, rows):
        if not data_list[2 + index * columns] in list_filters:
            list_filters.append(data_list[2 + index * columns])

    date_filter[str(file_name)] = list_filters

    remove_index = []
    for index in range(1, rows):
        if str(data_list[7 + index * columns]) == 'INDEF' or str(data_list[8 + index * columns]) == 'INDEF':
            remove_index.append(data_list[0 + index * columns])

    new_remove_index = list(set(remove_index))
    new_remove_index.sort()
    usable_stars = no_of_stars - len(new_remove_index)

    for index in range(1, rows):
        if not data_list[0 + index * columns] in new_remove_index:
            if data_list[2 + index * columns] == '7BesU':
                list_u.append(float(data_list[11 + index * columns]))
                list_u_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '6BesB':
                list_b.append(float(data_list[11 + index * columns]))
                list_b_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '5BesV':
                list_v.append(float(data_list[11 + index * columns]))
                list_v_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '4BesR':
                list_r.append(float(data_list[11 + index * columns]))
                list_r_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '3BesI':
                list_i.append(float(data_list[11 + index * columns]))
                list_i_err.append(float(data_list[12 + index * columns]))

    for i in range(0, usable_stars):
        if '6BesB' in list_filters and '5BesV' in list_filters:
            list_b_v.append(float("%6.4f" % (list_b[i] - list_v[i])))
            list_b_v_err.append(float("%6.4f" % (list_b_err[i] + list_v_err[i])))

        if '7BesU' in list_filters and '6BesB' in list_filters:
            list_u_b.append(float("%6.4f" % (list_u[i] - list_b[i])))
            list_u_b_err.append(float("%6.4f" % (list_u_err[i] + list_b_err[i])))

        if '5BesV' in list_filters and '4BesR' in list_filters:
            list_v_r.append(float("%6.4f" % (list_v[i] - list_r[i])))
            list_v_r_err.append(float("%6.4f" % (list_v_err[i] + list_r_err[i])))

        if '4BesR' in list_filters and '3BesI' in list_filters:
            list_r_i.append(float("%6.4f" % (list_r[i] - list_i[i])))
            list_r_i_err.append(float("%6.4f" % (list_r_err[i] + list_i_err[i])))

        if '5BesV' in list_filters and '3BesI' in list_filters:
            list_v_i.append(float("%6.4f" % (list_v[i] - list_i[i])))
            list_v_i_err.append(float("%6.4f" % (list_v_err[i] + list_i_err[i])))

    beta_b_v = []
    beta_u_b = []
    beta_v_r = []
    beta_r_i = []
    beta_v_i = []
    beta_v_1 = []
    beta_v_2 = []

    for i in range(0, usable_stars):
        if '6BesB' in list_filters and '5BesV' in list_filters:
            beta_b_v.append(float(mag_matrix[i][0]) - float(list_alpha[0]) * float(list_b_v[i]))

        if '7BesU' in list_filters and '6BesB' in list_filters:
            beta_u_b.append(float(mag_matrix[i][1]) - float(list_alpha[1]) * float(list_u_b[i]))

        if '5BesV' in list_filters and '4BesR' in list_filters:
            beta_v_r.append(float(mag_matrix[i][2]) - float(list_alpha[2]) * float(list_v_r[i]))

        if '4BesR' in list_filters and '3BesI' in list_filters:
            beta_r_i.append(float(mag_matrix[i][3]) - float(list_alpha[3]) * float(list_r_i[i]))

        if '5BesV' in list_filters and '3BesI' in list_filters:
            beta_v_i.append(float(mag_matrix[i][4]) - float(list_alpha[4]) * float(list_v_i[i]))

        if '5BesV' in list_filters:
            beta_v_1.append(float(mag_matrix[i][5]) - float(list_v[i]) - float(list_b_v[i]) -
                            float(list_alpha[5]) * float(mag_matrix[i][0]))
            # print mag_matrix[i][5]
            # print list_v[i]
            # print list_b_v[i]
            # print list_alpha[5]
            # print mag_matrix[i][0]
            # print "\n"
        if '5BesV' in list_filters:
            beta_v_2.append(float(mag_matrix[i][6]) - float(list_v[i]) - float(list_v_r[i]) -
                            float(list_alpha[6]) * float(mag_matrix[i][2]))

    master_beta = []
    for beta_list in [beta_b_v, beta_u_b, beta_v_r, beta_r_i, beta_v_i, beta_v_1, beta_v_2]:
        if len(beta_list) != 0:
            beta_list = reject(beta_list)
            beta_list = reject(beta_list)
            beta_list = reject(beta_list)
            master_beta.append("%5.4f " % np.mean(beta_list))
        else:
            master_beta.append("INDEF")
    date_filter[str(file_name)] = list_filters

    return master_beta


def calculate_beta_mag(file_name, mag_matrix):
    with open(str(file_name), "r") as f:
        data_list = f.read().split()
    with open(str(file_name), "r") as f:
        data_line = f.readline().split()

    columns = len(data_line)
    length_data = len(data_list)
    rows = length_data / columns

    list_u = []
    list_b = []
    list_v = []
    list_r = []
    list_i = []

    list_u_err = []
    list_b_err = []
    list_v_err = []
    list_r_err = []
    list_i_err = []

    list_filters = []
    for index in range(1, rows):
        if not data_list[2 + index * columns] in list_filters:
            list_filters.append(data_list[2 + index * columns])

    date_filter[str(file_name)] = list_filters

    remove_index = []
    for index in range(1, rows):
        if str(data_list[7 + index * columns]) == 'INDEF' or str(data_list[8 + index * columns]) == 'INDEF':
            remove_index.append(data_list[0 + index * columns])

    new_remove_index = list(set(remove_index))
    new_remove_index.sort()
    usable_stars = no_of_stars - len(new_remove_index)

    for index in range(1, rows):
        if not data_list[0 + index * columns] in new_remove_index:
            if data_list[2 + index * columns] == '7BesU':
                list_u.append(float(data_list[11 + index * columns]))
                list_u_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '6BesB':
                list_b.append(float(data_list[11 + index * columns]))
                list_b_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '5BesV':
                list_v.append(float(data_list[11 + index * columns]))
                list_v_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '4BesR':
                list_r.append(float(data_list[11 + index * columns]))
                list_r_err.append(float(data_list[12 + index * columns]))

            elif data_list[2 + index * columns] == '3BesI':
                list_i.append(float(data_list[11 + index * columns]))
                list_i_err.append(float(data_list[12 + index * columns]))

    beta_u = []
    beta_b = []
    beta_v = []
    beta_r = []
    beta_i = []

    for i in range(0, usable_stars):
        if '7BesU' in list_filters:
            beta_u.append(float(mag_matrix[i][0]) - float(list_u[i]))

        if '6BesB' in list_filters:
            beta_b.append(float(mag_matrix[i][1]) - float(list_b[i]))

        if '5BesV' in list_filters:
            beta_v.append(float(mag_matrix[i][2]) - float(list_v[i]))

        if '4BesR' in list_filters:
            beta_r.append(float(mag_matrix[i][3]) - float(list_r[i]))

        if '3BesI' in list_filters:
            beta_i.append(float(mag_matrix[i][4]) - float(list_i[i]))

    master_beta = []
    for beta_list in [beta_u, beta_b, beta_v, beta_r, beta_i]:
        if len(beta_list) != 0:
            beta_list = reject(beta_list)
            beta_list = reject(beta_list)
            beta_list = reject(beta_list)
            master_beta.append("%6.4f " % np.mean(beta_list))
        else:
            master_beta.append("INDEF")
    date_filter[str(file_name)] = list_filters

    return master_beta
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Gets Instrumental Magnitudes For All The Days From Individual Date-Wise TXDUMP'ed MAG Files
# ------------------------------------------------------------------------------------------------------------------- #
for date in list_dates:
    get_instr_mag(date)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Appends Names Of Instrumental Magnitude Files Onto A List "list_instr_mag4"
# ------------------------------------------------------------------------------------------------------------------- #
list_instr_mag4 = group_similar_files("list_instr_mag4", "output_instr_*_mag4")
list_alpha = [0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]
date_filter = {}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculates Beta Matrix From Beta Values For All The Days On Which Supernova Was Observed
# ------------------------------------------------------------------------------------------------------------------- #
super_master_beta = []
for file_name in list_instr_mag4:
    super_master_beta.append(calculate_beta_mag(file_name, mag_mean))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Performs TXDUMP Task On List Of ALS Files
# ------------------------------------------------------------------------------------------------------------------- #
for date in list_dates:
    group_similar_files("list_" + str(date) + "_psfmag1", str(date) + "*.als.1")
    txdump("list_" + str(date) + "_psfmag1", "output_" + str(date) + "_psfmag1", "als")

list_psfmag1 = group_similar_files("list_files_psfmag1", "output_*_psfmag1")
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Calculating Supernova Magnitudes
# ------------------------------------------------------------------------------------------------------------------- #

# def calculate_supernova_mag_from_color(list_color, v):
#     choice_v = v
#
#     u_mag = float("%6.4f" % (list_color[0] + list_color[1] + list_color[choice_v]))
#     b_mag = float("%6.4f" % (list_color[0] + list_color[choice_v]))
#     v_mag = float("%6.4f" % (list_color[choice_v]))
#     r_mag = float("%6.4f" % (list_color[choice_v] - list_color[2]))
#     i_mag = float("%6.4f" % (list_color[choice_v] - list_color[4]))
#     mag_supernova.append([u_mag, b_mag, v_mag, r_mag, i_mag])
#
#     return mag_supernova


def analyse_psfmag(file_name):
    with open(str(file_name), "r") as f:
        data_list_psf = f.read().split()
    with open(str(file_name), "r") as f:
        data_line_psf = f.readline().split()

    date = str(file_name)[7:12]
    index_date = list_dates.index(str(date))
    columns_psf = len(data_line_psf)
    length_data_psf = len(data_list_psf)
    rows_psf = length_data_psf / columns_psf
    supernova_index = no_of_stars + 1

    value_u = 0
    value_b = 0
    value_v = 0
    value_r = 0
    value_i = 0

    value_u_err = 0
    value_b_err = 0
    value_v_err = 0
    value_r_err = 0
    value_i_err = 0

    list_filters = []
    for index in range(1, rows_psf):
        if not data_list_psf[1 + index * columns_psf] in list_filters:
            list_filters.append(data_list_psf[1 + index * columns_psf])

    for index in range(1, rows_psf):
        if int(data_list_psf[2 + index * columns_psf]) == supernova_index:
            if data_list_psf[1 + index * columns_psf] == '7BesU':
                value_u = (float(data_list_psf[5 + index * columns_psf]))
                value_u_err = (float(data_list_psf[6 + index * columns_psf]))

            elif data_list_psf[1 + index * columns_psf] == '6BesB':
                value_b = (float(data_list_psf[5 + index * columns_psf]))
                value_b_err = (float(data_list_psf[6 + index * columns_psf]))

            elif data_list_psf[1 + index * columns_psf] == '5BesV':
                value_v = (float(data_list_psf[5 + index * columns_psf]))
                value_v_err = (float(data_list_psf[6 + index * columns_psf]))

            elif data_list_psf[1 + index * columns_psf] == '4BesR':
                value_r = (float(data_list_psf[5 + index * columns_psf]))
                value_r_err = (float(data_list_psf[6 + index * columns_psf]))

            elif data_list_psf[1 + index * columns_psf] == '3BesI':
                value_i = (float(data_list_psf[5 + index * columns_psf]))
                value_i_err = (float(data_list_psf[6 + index * columns_psf]))

    mag_u = "INDEF"
    mag_b = "INDEF"
    mag_v = "INDEF"
    mag_r = "INDEF"
    mag_i = "INDEF"

    mag_u_err = "INDEF"
    mag_b_err = "INDEF"
    mag_v_err = "INDEF"
    mag_r_err = "INDEF"
    mag_i_err = "INDEF"

    if '7BesU' in list_filters:
        mag_u = "%6.4f" % (float(value_u) + float(super_master_beta[index_date][0]))
        mag_u_err = "%6.4f" % (float(value_u_err))

    if '6BesB' in list_filters:
        mag_b = "%6.4f" % (float(value_b) + float(super_master_beta[index_date][1]))
        mag_b_err = "%6.4f" % (float(value_b_err))

    if '5BesV' in list_filters:
        mag_v = "%6.4f" % (float(value_v) + float(super_master_beta[index_date][2]))
        mag_v_err = "%6.4f" % (float(value_v_err))

    if '4BesR' in list_filters:
        mag_r = "%6.4f" % (float(value_r) + float(super_master_beta[index_date][3]))
        mag_r_err = "%6.4f" % (float(value_r_err))

    if '3BesI' in list_filters:
        mag_i = "%6.4f" % (float(value_i) + float(super_master_beta[index_date][4]))
        mag_i_err = "%6.4f" % (float(value_i_err))

    list_mag = [mag_u, mag_b, mag_v, mag_r, mag_i]
    list_mag_err = [mag_u_err, mag_b_err, mag_v_err, mag_r_err, mag_i_err]

    return [list_mag, list_mag_err]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Supernova Magnitudes For All The Epochs
# ------------------------------------------------------------------------------------------------------------------- #
mag_supernova = []
mag_supernova_err = []

for file_name in list_psfmag1:
    mag_supernova.append(analyse_psfmag(file_name)[0])
    mag_supernova_err.append(analyse_psfmag(file_name)[1])

mag_supernova_u = []
mag_supernova_b = []
mag_supernova_v = []
mag_supernova_r = []
mag_supernova_i = []

mag_supernova_u_err = []
mag_supernova_b_err = []
mag_supernova_v_err = []
mag_supernova_r_err = []
mag_supernova_i_err = []

for index in range(0, len(mag_supernova)):
    if mag_supernova[index][0] != 'INDEF':
        mag_supernova_u.append(float(mag_supernova[index][0]))
        mag_supernova_u_err.append(float(mag_supernova_err[index][0]))

    if mag_supernova[index][1] != 'INDEF':
        mag_supernova_b.append(float(mag_supernova[index][1]))
        mag_supernova_b_err.append(float(mag_supernova_err[index][1]))

    if mag_supernova[index][2] != 'INDEF':
        mag_supernova_v.append(float(mag_supernova[index][2]))
        mag_supernova_v_err.append(float(mag_supernova_err[index][2]))

    if mag_supernova[index][3] != 'INDEF':
        mag_supernova_r.append(float(mag_supernova[index][3]))
        mag_supernova_r_err.append(float(mag_supernova_err[index][3]))

    if mag_supernova[index][4] != 'INDEF':
        mag_supernova_i.append(float(mag_supernova[index][4]))
        mag_supernova_i_err.append(float(mag_supernova_err[index][4]))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Obtains Julian Day From The FITS Header
# ------------------------------------------------------------------------------------------------------------------- #
file_airmass = "/home/avinash/Dropbox/PyCharm/Reduction_Pipeline/airmas.asc"
year = 2014
month_dict = {'jan': '01', 'feb': '02', 'mar': '03', 'apr': '04', 'may': '05', 'jun': '06', 'jul': '07',
              'aug': '08', 'sep': '09', 'oct': '10', 'nov': '11', 'dec': '12'}

date_dict = {}
for i in range(0, len(list_dates)):
    if list_dates[i][0:3] != 'jan':
        date_dict[list_dates[i]] = str(year) + "-" + str(month_dict[list_dates[i][0:3]]) + \
                                   "-" + str(list_dates[i][3:5])
    else:
        date_dict[list_dates[i]] = str(year + 1) + "-" + str(month_dict[list_dates[i][0:3]]) + \
                                   "-" + str(list_dates[i][3:5])


def edit_file(file_name):
    with open(file_airmass, 'r') as f:
        text = f.readline()
        rest_data = f.read()

    # DATE_OBS = '2014-09-16'
    new_text = text[0:11] + "'" + date_dict[str(file_name)[0:5]] + "'" + "\n"
    with open(file_airmass, 'w') as f:
        f.write(new_text)
        f.write(rest_data)


def calc_julian_day(file_name):
    edit_file(file_name)
    iraf.noao(_doprint=0)
    iraf.astutil(_doprint=0)
    iraf.astutil.asthedit.verbose = "no"
    iraf.astutil.asthedit.update = "yes"
    iraf.astutil.asthedit.oldstyle = "no"
    iraf.astutil.asthedit(images=str(file_name), commands=file_airmass, Stdout="output_asthedit")

for file_name in list_fits:
    calc_julian_day(file_name)


list_fits_u = group_similar_files("", "*-u*.fits", exceptions="u1.fits.,u.fits.")
list_fits_b = group_similar_files("", "*-b*.fits", exceptions="b1.fits.,b.fits.,b2.fits.")
list_fits_v = group_similar_files("", "*-v*.fits", exceptions="v1.fits.,v.fits.,v2.fits.")
list_fits_r = group_similar_files("", "*-r*.fits", exceptions="r1.fits.,r.fits.,r2.fits.")
list_fits_i = group_similar_files("", "*-i*.fits", exceptions="i1.fits.,i.fits.,i2.fits.")

list_JD_u = []
list_JD_b = []
list_JD_v = []
list_JD_r = []
list_JD_i = []

for file_name in list_fits:
    fits_header = fits.getheader(str(file_name))
    filter_value = fits_header['FILTER']
    if filter_value == '7BesU' or filter_value == '7 Bes U':
        list_JD_u.append(fits_header['JD'])
    elif filter_value == '6BesB' or filter_value == '6 Bes B':
        list_JD_b.append(fits_header['JD'])
    elif filter_value == '5BesV' or filter_value == '5 Bes V':
        list_JD_v.append(fits_header['JD'])
    elif filter_value == '4BesR' or filter_value == '4 Bes R':
        list_JD_r.append(fits_header['JD'])
    elif filter_value == '3BesI' or filter_value == '3 Bes I':
        list_JD_i.append(fits_header['JD'])


def sort_sne_mag(list_jd, list_mag, list_mag_err):
    jd_mag_err = zip(list_jd, list_mag, list_mag_err)
    list_sorted_jd = [jd for jd, mag, err in sorted(jd_mag_err, key=lambda x: x[0])]
    list_sorted_mag = [mag for jd, mag, err in sorted(jd_mag_err, key=lambda x: x[0])]
    list_sorted_err = [err for jd, mag, err in sorted(jd_mag_err, key=lambda x: x[0])]
    return list_sorted_jd, list_sorted_mag, list_sorted_err


list_sorted_JD_u, list_sorted_mag_u, list_sorted_err_u = sort_sne_mag(list_JD_u, mag_supernova_u, mag_supernova_u_err)
list_sorted_JD_b, list_sorted_mag_b, list_sorted_err_b = sort_sne_mag(list_JD_b, mag_supernova_b, mag_supernova_b_err)
list_sorted_JD_v, list_sorted_mag_v, list_sorted_err_v = sort_sne_mag(list_JD_v, mag_supernova_v, mag_supernova_v_err)
list_sorted_JD_r, list_sorted_mag_r, list_sorted_err_r = sort_sne_mag(list_JD_r, mag_supernova_r, mag_supernova_r_err)
list_sorted_JD_i, list_sorted_mag_i, list_sorted_err_i = sort_sne_mag(list_JD_i, mag_supernova_i, mag_supernova_i_err)

JD_offset = min(list_JD_v)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print Supernova Magnitudes In A File
# ------------------------------------------------------------------------------------------------------------------- #
def write_sne_mag(list_jd, list_sne_mag, output_file_name):
    with open(str(output_file_name), "w") as f:
        f.write(str('Julian_Day') + " " * 6 + str('JD_offset') + " " * 2 + str('SNe_Mag') + "\n")

        for i in range(0, len(list_sne_mag)):
            f.write(str("%9.3f" % list_jd[i]) + "     " +
                    str("%7.3f" % (float(list_jd[i]) - float(JD_offset))) + " " * 5 +
                    str("%6.3f" % list_sne_mag[i]) + "\n")

write_sne_mag(list_sorted_JD_u, list_sorted_mag_u, output_file_name='OUTPUT_supernova_mag_u')
write_sne_mag(list_sorted_JD_b, list_sorted_mag_b, output_file_name='OUTPUT_supernova_mag_b')
write_sne_mag(list_sorted_JD_v, list_sorted_mag_v, output_file_name='OUTPUT_supernova_mag_v')
write_sne_mag(list_sorted_JD_r, list_sorted_mag_r, output_file_name='OUTPUT_supernova_mag_r')
write_sne_mag(list_sorted_JD_i, list_sorted_mag_i, output_file_name='OUTPUT_supernova_mag_i')

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot Light Curve
# ------------------------------------------------------------------------------------------------------------------- #
list_offset_mag_u = [float(value) + 1.5 for value in list_sorted_mag_u]
list_offset_mag_b = [float(value) + 0.0 for value in list_sorted_mag_b]
list_offset_mag_v = [float(value) - 0.8 for value in list_sorted_mag_v]
list_offset_mag_r = [float(value) - 1.6 for value in list_sorted_mag_r]
list_offset_mag_i = [float(value) - 2.5 for value in list_sorted_mag_i]

list_offset_JD_u = [float(value) - float(JD_offset) for value in list_sorted_JD_u]
list_offset_JD_b = [float(value) - float(JD_offset) for value in list_sorted_JD_b]
list_offset_JD_v = [float(value) - float(JD_offset) for value in list_sorted_JD_v]
list_offset_JD_r = [float(value) - float(JD_offset) for value in list_sorted_JD_r]
list_offset_JD_i = [float(value) - float(JD_offset) for value in list_sorted_JD_i]

# plt.plot(list_offset_JD_u, list_offset_mag_u, 'ko', markersize=8, label='$U+1.5$')
# plt.errorbar(list_offset_JD_u, list_offset_mag_u, yerr=list_sorted_err_u, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_b, list_offset_mag_b, 'bD', markersize=8, label='$B$')
# plt.errorbar(list_offset_JD_b, list_offset_mag_b, yerr=list_sorted_err_b, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_v, list_offset_mag_v, 'gs', markersize=8, label='$V-0.8$')
# plt.errorbar(list_offset_JD_v, list_offset_mag_v, yerr=list_sorted_err_v, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_r, list_offset_mag_r, 'r^', markersize=8, label='$R-1.6$')
# plt.errorbar(list_offset_JD_r, list_offset_mag_r, yerr=list_sorted_err_r, capsize=4, ls='none', color='k')
#
# plt.plot(list_offset_JD_i, list_offset_mag_i, 'mp', markersize=8, label='$I-2.5$')
# plt.errorbar(list_offset_JD_i, list_offset_mag_i, yerr=list_sorted_err_i, capsize=4, ls='none', color='k')
#
# plt.xlabel('$\Delta$t (Days Since $V$-Band Maximum)', fontsize=18)
# plt.ylabel('Apparent Magnitudes', fontsize=18)
# plt.xlim(-20, 200)
# plt.ylim(24.0, 11.0)
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)
# plt.legend(numpoints=1, fontsize=15)
# plt.savefig("light_curve_final.png")
# plt.title("Broadband Photometric Light Curves", fontsize=20)
# plt.grid()
# plt.show()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Fit Spline To The Light Curves
# ------------------------------------------------------------------------------------------------------------------- #
plt.plot(list_offset_JD_u, list_sorted_mag_u, 'ko', markersize=8, label='$U$')
plt.errorbar(list_offset_JD_u, list_sorted_mag_u, yerr=list_sorted_err_u, capsize=4, ls='none', color='k')

plt.plot(list_offset_JD_b, list_sorted_mag_b, 'bD', markersize=8, label='$B$')
plt.errorbar(list_offset_JD_b, list_sorted_mag_b, yerr=list_sorted_err_b, capsize=4, ls='none', color='k')

plt.plot(list_offset_JD_v, list_sorted_mag_v, 'gs', markersize=8, label='$V$')
plt.errorbar(list_offset_JD_v, list_sorted_mag_v, yerr=list_sorted_err_v, capsize=4, ls='none', color='k')

plt.plot(list_offset_JD_r, list_sorted_mag_r, 'r^', markersize=8, label='$R$')
plt.errorbar(list_offset_JD_r, list_sorted_mag_r, yerr=list_sorted_err_r, capsize=4, ls='none', color='k')

plt.plot(list_offset_JD_i, list_sorted_mag_i, 'mp', markersize=8, label='$I$')
plt.errorbar(list_offset_JD_i, list_sorted_mag_i, yerr=list_sorted_err_i, capsize=4, ls='none', color='k')


def fit_spline1d(list_days, list_mag):
    xaxis = np.linspace(min(list_days), max(list_days), 1000)
    f1 = interpolate.interp1d(list_days, list_mag, kind='cubic')
    plt.plot(xaxis, f1(xaxis), 'k-')


def fit_curve_3d(list_days, list_mag):
    xnew = np.linspace(min(list_days) - 2, max(list_days) + 2, 1000)
    spl = interpolate.UnivariateSpline(list_days, list_mag, k=3, s=0.04)
    plt.plot(xnew, spl(xnew), 'r-')


def fit_curve_2d(list_days, list_mag):
    xaxis = np.linspace(min(list_days), max(list_days), 1000)
    spline = interpolate.UnivariateSpline(list_days, list_mag, k=2, s=0.1)
    plt.plot(xaxis, spline(xaxis), 'k-')


def fit_curve_1d(list_days, list_mag):
    xaxis = np.linspace(min(list_days), max(list_days), 1000)
    spline = interpolate.UnivariateSpline(list_days, list_mag, k=1, s=0.01)
    plt.plot(xaxis, spline(xaxis), 'g-')


def delete_element(list_of_lists, list_elements):
    for element in list_elements:
        for i in range(0, len(list_of_lists)):
            list_of_lists[i] = [x for x in list_of_lists[i] if list_of_lists[i].index(x) != int(element)]
    return list_of_lists


[list_offset_JD_u, list_sorted_mag_u, list_sorted_err_u] = delete_element(
    [list_offset_JD_u, list_sorted_mag_u, list_sorted_err_u], list_elements=[len(list_offset_JD_u) - 2])

[list_offset_JD_b, list_sorted_mag_b, list_sorted_err_b] = delete_element(
    [list_offset_JD_b, list_sorted_mag_b, list_sorted_err_b], list_elements=[len(list_offset_JD_b) - 1])

[list_offset_JD_v, list_sorted_mag_v, list_sorted_err_v] = delete_element(
    [list_offset_JD_v, list_sorted_mag_v, list_sorted_err_v], list_elements=[len(list_offset_JD_r) - 4])

[list_offset_JD_r, list_sorted_mag_r, list_sorted_err_r] = delete_element(
    [list_offset_JD_r, list_sorted_mag_r, list_sorted_err_r], list_elements=[len(list_offset_JD_r) - 2, 10])

[list_offset_JD_i, list_sorted_mag_i, list_sorted_err_i] = delete_element(
    [list_offset_JD_i, list_sorted_mag_i, list_sorted_err_i], list_elements=[20])

# list_list_JD = [list_sorted_JD_u, list_sorted_JD_b, list_sorted_JD_v, list_sorted_JD_r, list_sorted_JD_i]
list_list_JD = [list_offset_JD_u, list_offset_JD_b, list_offset_JD_v, list_offset_JD_r, list_offset_JD_i]
list_list_mag = [list_sorted_mag_u, list_sorted_mag_b, list_sorted_mag_v, list_sorted_mag_r, list_sorted_mag_i]

for index in range(0, len(list_list_JD)):
    # fit_spline1d(list_list_JD[index], list_list_mag[index])
    fit_curve_3d(list_list_JD[index], list_list_mag[index])
    fit_curve_2d(list_list_JD[index], list_list_mag[index])
    fit_curve_1d(list_list_JD[index], list_list_mag[index])

plt.xlabel('$\Delta$t (Days Since $V$-Band Maximum)', fontsize=18)
plt.ylabel('Apparent Magnitudes', fontsize=18)
# plt.xlim(-10, 200)
plt.ylim(22.0, 14.0)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.legend(numpoints=1, fontsize=15)
plt.title("Broadband Photometric Light Curves", fontsize=20)
plt.grid()
plt.savefig('fit_light_curve.png')
plt.show()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print Interpolated Supernova Magnitudes In OUTPUT Files For Different Bands
# ------------------------------------------------------------------------------------------------------------------- #
def write_interpolated_sne_mag(list_jd, list_mag, output_file_name):
    with open(str(output_file_name), "w") as f:
        f.write(str('Julian_Day') + " " * 6 + str('JD_offset') + " " * 2 + str('SNe_Mag') + "\n")
        jd_array = np.arange(int(min(list_jd)), int(max(list_jd)), 0.5)
        spl = interpolate.UnivariateSpline(list_jd, list_mag, k=2, s=0.1)
        for i in range(0, len(jd_array)):
            f.write(str("%9.3f" % (float(jd_array[i]) + float(JD_offset))) + "     " +
                    str("%7.3f" % (float(jd_array[i]))) + " " * 5 +
                    str("%6.3f" % float(spl(jd_array[i]))) + "\n")

write_interpolated_sne_mag(list_offset_JD_u, list_sorted_mag_u, output_file_name='OUTPUT_interp_snmag_u')
write_interpolated_sne_mag(list_offset_JD_b, list_sorted_mag_b, output_file_name='OUTPUT_interp_snmag_b')
write_interpolated_sne_mag(list_offset_JD_v, list_sorted_mag_v, output_file_name='OUTPUT_interp_snmag_v')
write_interpolated_sne_mag(list_offset_JD_r, list_sorted_mag_r, output_file_name='OUTPUT_interp_snmag_r')
write_interpolated_sne_mag(list_offset_JD_i, list_sorted_mag_i, output_file_name='OUTPUT_interp_snmag_i')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Print Interpolated Supernova Magnitudes In A Single File
# ------------------------------------------------------------------------------------------------------------------- #
# with open('OUTPUT_interp_snmag', "w") as f:
#     f.write(str('Julian_Day') + " " * 6 + str('JD_offset') + " " * 2 + str('SNe_Mag') + "\n")
#     list_spl = []
#     for list_jd in list_list_JD:
#         for list_mag in list_list_mag:
#             list_spl.append(interpolate.UnivariateSpline(list_jd, list_mag, k=3, s=0.04))
#
#             jd_array = np.arange(int(min(list_jd)), int(max(list_jd)), 0.5)
#     for i in range(0, len(jd_array)):
#         for list_mag
#         f.write(str("%9.3f" % float(jd_array[i])) + "     " +
#                 str("%7.3f" % (float(jd_array[i]) - float(JD_offset))) + " " * 5 +
#                 str("%6.3f" % float(spl(jd_array[i]))) + "\n")
# ------------------------------------------------------------------------------------------------------------------- #


# # ------------------------------------------------------------------------------------------------------------------- #
# # Calculate Flux Values In Individual Bands Day-Wise
# # ------------------------------------------------------------------------------------------------------------------- #
# # mag_lambda = -2.5log(flux_lambda) - 21.100 - zp(f_lambda)
# # Zero Point        U       B       V       R       I
# # zp(f_lambda)   -0.152  -0.602   0.000   0.555   1.271
# zp_u = -0.152
# zp_b = -0.602
# zp_v = 0.000
# zp_r = 0.555
# zp_i = 1.271
# A_lambda_u = 0.303
# A_lambda_b = 0.254
# A_lambda_v = 0.192
# A_lambda_r = 0.152
# A_lambda_i = 0.105
# list_zp = [zp_u, zp_b, zp_v, zp_r, zp_i]
#
#
# def calculate_flux(mag_supernova_day):
#     flux_u = 'INDEF'
#     flux_b = 'INDEF'
#     flux_v = 'INDEF'
#     flux_r = 'INDEF'
#     flux_i = 'INDEF'
#
#     if mag_supernova_day[0] != 'INDEF':
#         flux_u = 10 ** (-(float(mag_supernova_day[0]) - A_lambda_u + zp_u + 21.100) / 2.5)
#     if mag_supernova_day[1] != 'INDEF':
#         flux_b = 10 ** (-(float(mag_supernova_day[1]) - A_lambda_b + zp_b + 21.100) / 2.5)
#     if mag_supernova_day[2] != 'INDEF':
#         flux_v = 10 ** (-(float(mag_supernova_day[2]) - A_lambda_v + zp_v + 21.100) / 2.5)
#     if mag_supernova_day[3] != 'INDEF':
#         flux_r = 10 ** (-(float(mag_supernova_day[3]) - A_lambda_r + zp_r + 21.100) / 2.5)
#     if mag_supernova_day[4] != 'INDEF':
#         flux_i = 10 ** (-(float(mag_supernova_day[4]) - A_lambda_i + zp_i + 21.100) / 2.5)
#     return [flux_u, flux_b, flux_v, flux_r, flux_i]
#
# flux_supernova = []
# for index in range(0, len(mag_supernova)):
#     flux_supernova.append(calculate_flux(mag_supernova[index]))
#
#
# with open('output_flux_values' , 'w') as f:
#     for i in range(0, len(flux_supernova)):
#         for j in range(0, len(flux_supernova[i])):
#             f.write(list_JD[j] + ' ' + str(flux_supernova[i][j]) + '\n')
#         f.write('\n')
#
# # left_limit = 300
# # right_limit = 1100
# # list_band_centers = [370, 420, 530, 600, 800]
# # list_band_centers.sort()
# # test_y = [flux_supernova[0][0], flux_supernova[0][1], flux_supernova[0][2], flux_supernova[0][4]]
# # test_x = [list_band_centers[0], list_band_centers[1], list_band_centers[2], list_band_centers[4]]
# #
# # for i in range(0, len(flux_supernova)):
# #     print flux_supernova[i]
# # plt.plot(list_band_centers, flux_supernova[2],'o')
# # spline = interpolate.UnivariateSpline(list_band_centers, flux_supernova[0], bbox=[300, 800], k=4)
# # spline2 = interpolate.interp1d(list_band_centers, flux_supernova[0], kind='linear')
# # spline3 = interpolate.splrep(test_x, test_y, s=2, k=3)
# # xaxis = np.linspace(200, 1000, 1000)
# # ynew = interpolate.splev(xaxis, spline3, der=0)
# # plt.plot(xaxis, ynew)
# # plt.ylim(-0.1e-15, 2e-15)
# # plt.grid()
# # plt.show()
# # ------------------------------------------------------------------------------------------------------------------- #
#
# for i in range(0, len(list_dates)):
#     with open("output_flux_" + str(list_dates[i]), 'w') as f:
#         f.write(str(310) + " " + str(920) + " " + str(100) + "\n")
#         if flux_supernova[i][0] != 'INDEF':
#             f.write(str(list_band_centers[0]) + " " + str(flux_supernova[i][0]) + "\n")
#         if flux_supernova[i][1] != 'INDEF':
#             f.write(str(list_band_centers[1]) + " " + str(flux_supernova[i][1]) + "\n")
#         if flux_supernova[i][2] != 'INDEF':
#             f.write(str(list_band_centers[2]) + " " + str(flux_supernova[i][2]) + "\n")
#         if flux_supernova[i][3] != 'INDEF':
#             f.write(str(list_band_centers[3]) + " " + str(flux_supernova[i][3]) + "\n")
#         if flux_supernova[i][4] != 'INDEF':
#             f.write(str(list_band_centers[4]) + " " + str(flux_supernova[i][4]) + "\n")
#
# # print list_fits_v[14:17]
# # print list_fits_b[6:10]
# #
# # print mag_supernova_v[14:17]
# # print mag_supernova_b[6:10]
# # trapz(a[:, 1], x=a[:, 0])
# bol_flux_dec19 = 2.672e-13
# bol_flux_dec26 = 2.287e-13
# bol_flux_jan03 = 1.973e-13
# bol_flux_jan25 = 1.794e-13
# bol_flux_jul25 = 5.446e-12
# bol_flux_jul16 = 7.703e-12
# bol_flux_175 = 2.25e-13
# dist = 45 * 3.086e+24
# z = 0.010424
#
# test_lum1 = 4 * np.pi * dist**2 * bol_flux_dec19
# test_lum2 = 4 * np.pi * dist**2 * bol_flux_jan25
# test_mag1 = -2.5 * math.log10(test_lum1 / 3.839e+33) + 4.72
# test_mag2 = -2.5 * math.log10(test_lum2 / 3.839e+33) + 4.72
# bol_lum_175 = bol_flux_175 * 4 * np.pi * (dist ** 2)
# print "bol_lum_dec19 = ", test_lum1
# print "bol_lum_jan25 = ", test_lum2
# print "\n"
# print -0.4 * math.log10(test_lum1)
#
# print "bol_mag_dec19 = ", test_mag1
# print "bol_mag_jan25 = ", test_mag2
# print "\n"
#
# m_ni2 = (bol_lum_175 * 0.07 / 9.92e41) / (math.exp(-175/111.4) - math.exp(-175/8.8))
# nickel_mass = 7.866e-44 * bol_lum_175 * math.exp(((175 * (1 + z) - 6.1) / 111.26))
#
# print "jerkstrand_Ni = ", m_ni2
# print "hamuy_ni =", nickel_mass
#
# y = [bol_flux_dec19, bol_flux_dec26, bol_flux_jan03, bol_flux_jan25]
# x = [list_JD_v[10] + 11, list_JD_v[11] + 11, list_JD_v[12] + 11, list_JD_v[13] + 11]
# # plt.plot(x, y, 'o')
# # y2 = interpolate.interp1d(x, y, kind='linear')
# # plt.plot(x, y2(x), '-')
# # # plt.show()
#
# peak_bolometric_lum = 4 * np.pi * dist**2 * bol_flux_jul16
# peak_bolometric_mag = -2.5 * math.log10(peak_bolometric_lum / 3.839e+33) + 4.72
# print -0.4 * math.log10(peak_bolometric_lum)
# print "peak bolometric lum = ", peak_bolometric_lum
# print "peak bolometric mag =", peak_bolometric_mag
#
# tail_lum = bol_lum_175
# print "tail_lum=", tail_lum
# tail_mag_et = 1e40 * 10**(0.8)
# print tail_lum/tail_mag_et * 0.048
#
# lum_175 = 10**(-((20.0 - 0.26 - 33.2) - 4.72)/2.5) * 3.839e33
# print "tail_lum=", lum_175
#
# # v = 2000 kms-1
# # velociy around v-50
# # spiro et al 2014 mnras
# # Intro - Finding chart
# #
# # Light Curve
# # ColorCurves
# # fitting the bolometric flux
# # Estimation of ni mass
# # Nickel Mass - duration of plateau light curve
# # v -band steepness - slope of 4 points
# # nickel mass - radioactive tail
# # observations - type IIP supernova
# # duartion of the plateau
# #
# # vinko et al
# #
# # absolute v magntitude = 17 - 0.2 - 33.2
# # low luminosity type-IIP
# # 1.25 V-I
# # 0.60 V-R
# # 1.00 B-V
#
# list_steep = list_JD_v[23:27]
# mag_steep = mag_supernova_v[23:27]
#
# plt.plot(list_steep, mag_steep)
#
# slope = (mag_steep[3] - mag_steep[0]) / (list_steep[3] - list_steep[0])
# # print list_steep
# # print mag_steep
# print slope
# M_ni_2 = 10 ** (-6.2295 * slope - 0.8417)
# print M_ni_2

# s1rev = interpolate.InterpolatedUnivariateSpline(list_days, list_mag)
# s3 = interpolate.interp1d(list_days, list_mag, kind='cubic')
# rbf = interpolate.Rbf(list_days, list_mag)
# tck = interpolate.splrep(list_days, list_mag, s=1)
# ynew = interpolate.splev(xaxis, tck, der=0)
# x = np.arange(0, 2*np.pi+np.pi/4, 2*np.pi/8)
# y = np.sin(x)
# s = interpolate.InterpolatedUnivariateSpline(x, y)
# xnew = np.arange(0, 2*np.pi, np.pi/50)
# ynew = s(xnew)
#
# plt.figure()
# plt.plot(x, y, 'x', xnew, ynew, xnew, np.sin(xnew), x, y, 'b')
# plt.legend(['Linear', 'InterpolatedUnivariateSpline', 'True'])
# plt.axis([-0.05, 6.33, -1.05, 1.05])
# plt.title('InterpolatedUnivariateSpline')
# plt.show()

# m_star - m_sun = -2.5*math.log10(f_star/f_sun)
# f_star = f_sun * 10 ** (-0.4 * (m_star - m_sun))
# f_star = 3.839e33 * 10 ** (-0.4 * (m - 4.72))


