# #############################################################################################################
#   To Fit The Curve Of Instrumental Color Magnitudes Vs Actual Color Magnitudes Plot (Obtain 'm' and 'c')   #
# #############################################################################################################

# ####---Import Required Libraries---#####

import math
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ####---Read The 'output_instr_mag2' File Generated After running 'zero_mag.py'---#####

f=open('/home/avinash/Supernovae_Data/ASASSN14dq/Standards/PG0942/output_instr_mag2')
data=f.read()
f.close()


# ####---Make The Data Usable By Splitting The Data Into A List Of Strings---#####

data1=data.split()


# ####---Initialize Global Variables---#####

length_data = len(data1)
rows = length_data/13


# ####---Initialize Lists To Be Printed In The Output File"---#####

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


for i in range(1, rows):

    # ####---Put Respective Band Magnitudes In A List---#####

    if str(data1[7 + i * 13]) != 'INDEF' and str(data1[8 + i * 13]) != 'INDEF':

        if data1[2 + i * 13] == '7BesU':
            star_id.append(int(data1[0+i*13]))
            list_u.append(float(data1[11 +i*13]))
            list_u_err.append(float(data1[12 +i*13]))

        elif data1[2 + i * 13] == '6BesB':
            list_b.append(float(data1[11 +i*13]))
            list_b_err.append(float(data1[12 +i*13]))

        elif data1[2 +i*13] == '5BesV':
            list_v.append(float(data1[11 +i*13]))
            list_v_err.append(float(data1[12 +i*13]))

        elif data1[2 + i*13] == '4BesR':
            list_r.append(float(data1[11 +i*13]))
            list_r_err.append(float(data1[12 +i*13]))

        elif data1[2 + i*13] == '3BesI':
            list_i.append(float(data1[11 +i*13]))
            list_i_err.append(float(data1[12 +i*13]))


count = len(list_u)


for i in range(0, count):

    # ####---Calculate Color Magnitudes From Individual Band Magnitudes---#####

    list_b_v.append(float("%6.4f" % (list_b[i] - list_v[i])))

    list_u_b.append(float("%6.4f" % (list_u[i] - list_b[i])))

    list_v_r.append(float("%6.4f" % (list_v[i] - list_r[i])))

    list_r_i.append(float("%6.4f" % (list_r[i] - list_i[i])))

    list_v_i.append(float("%6.4f" % (list_v[i] - list_i[i])))

    # ####---Calculate Error In Color Magnitude Calculation---#####

    list_b_v_err.append(float("%6.4f" % (list_b_err[i] + list_v_err[i])))

    list_u_b_err.append(float("%6.4f" % (list_u_err[i] + list_b_err[i])))

    list_v_r_err.append(float("%6.4f" % (list_v_err[i] + list_r_err[i])))

    list_r_i_err.append(float("%6.4f" % (list_r_err[i] + list_i_err[i])))

    list_v_i_err.append(float("%6.4f" % (list_v_err[i] + list_i_err[i])))


# ####---Determining The No. Of Standard Stars In The Field And The No. Of Days Of Obsv. Of The Standard Field---#####

star_count = int(star_id[-1])
no_of_days = len(star_id)/ star_count
standard_whole = []


# ####---Creating A Matrix Of Instrumental V-Band And Color Magnitudes For The Standard Field---#####

# ####---List Of 2 Elements, Each With 'star_counts' Lists Of 12 Elements Each---#####

for i in range(0, no_of_days):
    standard_day = []   # Local Definition Of 'standard_day' is important
    for j in range(0, star_count):
        standard_day.append([list_v[star_count*int(i):star_count*int(i+1)][j],
                             list_b_v[star_count*int(i):star_count*int(i+1)][j],
                             list_u_b[star_count*int(i):star_count*int(i+1)][j],
                             list_v_r[star_count*int(i):star_count*int(i+1)][j],
                             list_r_i[star_count*int(i):star_count*int(i+1)][j],
                             list_v_i[star_count*int(i):star_count*int(i+1)][j],
                             list_v_err[star_count*int(i):star_count*int(i+1)][j],
                             list_b_v_err[star_count*int(i):star_count*int(i+1)][j],
                             list_u_b_err[star_count*int(i):star_count*int(i+1)][j],
                             list_v_r_err[star_count*int(i):star_count*int(i+1)][j],
                             list_r_i_err[star_count*int(i):star_count*int(i+1)][j],
                             list_v_i_err[star_count*int(i):star_count*int(i+1)][j]])

    standard_whole.append(standard_day)


# ####---List Of 2 Elements, Each With 12 Lists Of 'star_counts' Elements Each---#####
"""
for i in range(0, no_of_days):

    standard_day = [list_v[star_count*int(i):star_count*int(i+1)],
                list_b_v[star_count*int(i):star_count*int(i+1)],
                list_u_b[star_count*int(i):star_count*int(i+1)],
                list_v_r[star_count*int(i):star_count*int(i+1)],
                list_r_i[star_count*int(i):star_count*int(i+1)],
                list_v_i[star_count*int(i):star_count*int(i+1)],
                list_v_err[star_count*int(i):star_count*int(i+1)],
                list_b_v_err[star_count*int(i):star_count*int(i+1)],
                list_u_b_err[star_count*int(i):star_count*int(i+1)],
                list_v_r_err[star_count*int(i):star_count*int(i+1)],
                list_r_i_err[star_count*int(i):star_count*int(i+1)],
                list_v_i_err[star_count*int(i):star_count*int(i+1)]]

    standard_whole.append(standard_day)
"""

#####---V-Band Magnitude and Color Terms For Landolt Standards---#####

# Standard - PG0918+029
# Name	  V       B-V     U-B     V-R     R-I     V-I      V_Err   B-V_Err  U-B_Err  V-R_Err  R-I_Err  V-I_Err
# A 	14.490   0.536	-0.032   0.325   0.336   0.661	  0.0033   0.0058   0.0095   0.0039   0.0076   0.0085
# B	13.963	 0.765   0.366   0.417   0.370   0.787	  0.0034   0.0072   0.0159   0.0025   0.0045   0.0056
# C	13.537	 0.631   0.087   0.367   0.357   0.722	  0.0020   0.0028   0.0048   0.0015   0.0022   0.0028
# D	12.272   1.044   0.821   0.575   0.535   1.108	  0.0021   0.0030   0.0071   0.0016   0.0018   0.0018
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


PG0918 = [
    [14.490, 0.536, -0.032, 0.325, 0.336, 0.661, 0.0033, 0.0058, 0.0095, 0.0039, 0.0076, 0.0085],
    [13.963, 0.765, 0.366, 0.417, 0.370, 0.787, 0.0034, 0.0072, 0.0159, 0.0025, 0.0045, 0.0056],
    [13.537, 0.631, 0.087, 0.367, 0.357, 0.722, 0.0020, 0.0028, 0.0048, 0.0015, 0.0022, 0.0028],
    [12.272, 1.044, 0.821, 0.575, 0.535, 1.108, 0.0021, 0.0030, 0.0071, 0.0016, 0.0018, 0.0018],
    [13.327, -0.271, -1.081, -0.129, -0.159, -0.288, 0.0024, 0.0024, 0.0030, 0.0019, 0.0055, 0.0063]]

PG0231 = [
    [12.772, 0.710, 0.270, 0.405, 0.394, 0.799, 0.0008, 0.0015, 0.0030, 0.0011, 0.0030, 0.0030],
    [14.735, 1.448, 1.342, 0.954, 0.998, 1.951, 0.0030, 0.0072, 0.0178, 0.0034, 0.0026, 0.0057],
    [13.702, 0.671, 0.114, 0.399, 0.385, 0.783, 0.0014, 0.0078, 0.0148, 0.0028, 0.0064, 0.0085],
    [14.027, 1.088, 1.046, 0.675, 0.586, 1.256, 0.0029, 0.0075, 0.0312, 0.0081, 0.0064, 0.0110],
    [13.804, 0.677, 0.201, 0.390, 0.369, 0.757, 0.0046, 0.0040, 0.0075, 0.0035, 0.0017, 0.0023],
    [16.105, -0.329, -1.192, -0.162, -0.371, -0.534, 0.0068, 0.0083, 0.0045, 0.0276, 0.1066, 0.1221]]

PG0942 = [
    [14.731, 0.783, 0.339, 0.610, 0.477, 1.081, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
    [14.108, 0.525, 0.085, 0.368, 0.333, 0.697, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
    [14.989, 0.727, 0.369, 0.539, 0.376, 0.909, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
    [13.707, 0.564, 0.129, 0.348, 0.343, 0.687, 0.0025, 0.0028, 0.0075, 0.0039, 0.0022, 0.0042],
    [14.004, -0.294, -1.175, -0.130, -0.149, -0.280, 0.0045, 0.0056, 0.0069, 0.0069, 0.0120, 0.0144]]



#####---Function For Finding Out Which Standard Field Is Being Analyzed---#####

standard_field = str(str(data1[14])[11:17])

def standard(standard):
    if standard == 'PG0231' or standard == 'pg0231':
        return PG0231
    elif standard == 'PG0918' or standard == 'pg0918':
        return PG0918

#####---Plots To Determine Coefficients---#####

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

#####---Curve-Fitting And Determining The list of slopes 'alpha' and intercepts 'beta' for 7 different Plots---#####

list_x = [1, 2, 3, 4, 5, 1, 3]
list_y = [1, 2, 3, 4, 5, 0, 0]
list_alpha = [0.9090, 1.2678, 1.0160, 1.0161, 1.0094, 0.0442, 0.0748]
list_xlabel = ['(B-V)obs', '(U-B)obs', '(V-R)obs', '(R-I)obs', '(V-I)obs', '(B-V)', '(V-R)']
list_ylabel = ['(B-V)', '(U-B)', '(V-R)', '(R-I)', '(V-I)', 'V-Vobs', 'V-Vobs']
list_title = ['(B-V) Vs (B-V)obs', '(U-B) Vs (U-B)obs', '(V-R) Vs (V-R)obs', '(R-I) Vs (R-I)obs', '(V-I) Vs (V-I)obs',
              '(V-Vobs) Vs (B-V)', '(V-Vobs) Vs (V-R)']

colors_dot = ['ro', 'go', 'bo', 'ko', 'co', 'mo', 'yo']
colors_line = ['r-', 'g-', 'b-', 'k-', 'c-', 'm-', 'y-']


def func(x, alpha, beta):
    return alpha*x + beta

#def matrix_to_array(a):
#	return np.squeeze(np.asarray(a))

#fig, axs = plt.subplots(2,5, figsize=(15, 6), facecolor='w', edgecolor='k')
#axs = axs.rave1()


def list_statistics(list_input):
    value_mean = float(np.mean(list_input))
    value_median = float(np.median(list_input))
    value_std = float(np.std(list_input))

    return value_mean, value_median, value_std


def reject(list_input):
    diff = 0
    pop = False
    for i in range(0, len(list_input)):
        list_input[i] = float(list_input[i])

    list_input.sort()
    value_mean, value_median, value_std = list_statistics(list_input)

    if abs(list_input[0] - value_median) < abs(list_input[-1] - value_median):
        diff = -1
    else:  # Redundant
        diff = 0

    if abs(list_input[diff] - value_median) > value_std:
        list_input.pop(diff)
        pop = True

    if pop:
        value_mean, value_median, value_std = list_statistics(list_input)

        if abs(list_input[0] - value_median) < abs(list_input[-1] - value_median):
            diff = -1
        else:
            diff = 0

        if abs(list_input[diff] - value_median) > 2 * value_std:
            list_input.pop(diff)

    return list_input

list_beta = []
list_beta_err = []

for i in range(0, no_of_days):
    beta_list = []
    beta_err_list = []
    for k in range(0, len(list_x)):
        beta = []
        beta_err = []
        for j in range(0, star_count):
            x = 0
            y = 0
            x_err = 0
            y_err = 0
            if list_y[k] != 0:
                x = float(standard_whole[i][j][list_x[k]])
                y = float(standard(standard_field)[j][list_y[k]])
                x_err = float(standard_whole[i][j][list_x[k] + 6])
                y_err = float(standard(standard_field)[j][list_y[k] + 6])
            else:
                x = float(standard(standard_field)[j][list_x[k]])
                y = float(standard(standard_field)[j][list_y[k]] - standard_whole[i][j][list_y[k]])
                x_err = float(standard(standard_field)[j][list_x[k] + 6])
                y_err = math.sqrt(float(standard(standard_field)[j][list_y[k] + 6]) ** 2 +
                                  float(standard_whole[i][j][list_y[k] + 6]) ** 2)
            beta.append(y - (x * list_alpha[k]))
            beta_err.append(math.sqrt(y_err ** 2 + (x_err * list_alpha[k]) ** 2))
        beta = reject(beta)
        beta_list.append(float("%8.5f" % (np.mean(beta))))
        beta_err_list.append(float("%7.5f" % (np.mean(beta_err))))
    list_beta.append(beta_list)
    list_beta_err.append(beta_err_list)

print list_beta
print list_beta_err


"""
#####---Print The Calculated Data Along With Original Data In The Output File---#####

f2=open("output_color_mag2", "w")
f2.write(str('STAR_ID') + "     " + str('B-V') + "        " + str('B-V Error') + "      " + str('U-B') + "        " +
         str('U-B Error') + "      " + str('V-R') + "       " + str('V-R Error') + "     " + str('R-I') + "      " +
         str('R-I Error') + "      " + str('V-I') + "      " + str('V-I Error') + "\n" + "\n")

for i in range(count):
    f2.write(str("%3.0f" % star_id[i]) + "    " + str("%10.4f" % list_b_v[i]) + "   " +
             str("%10.4f" % list_b_v_err[i]) + "   " + str("%10.4f" % list_u_b[i]) + "   " +
             str("%10.4f" % list_u_b_err[i]) + "   " + str("%10.4f" % list_v_r[i]) + "  " +
             str("%10.4f" % list_v_r_err[i]) + "  " + str("%10.4f" % list_r_i[i]) + " " +
             str("%10.4f" % list_r_i_err[i]) + "  " + str("%10.4f" % list_v_i[i]) + "  " +
        str("%10.4f" % list_v_i_err[i]) + "\n")

f2.close()

def standard_std(list_alpha, list_beta, day):
    object_stand = []
    for j in range(0, star_count):
        x = []
        for k in range(0, len(list_y)):
            if list_y[k] != 0:
                x.append(float("%8.6f" % (list_alpha[k]*standard_whole[day][j][list_x[k]] + list_beta[day][k])))
            else:
                x.append(
                    float("%8.6f" % (standard_whole[day][j][list_y[k]] +
                    list_alpha[k] * (list_alpha[int(list_x[k] - 1)] * standard_whole[day][j][list_x[k]] +
                    list_beta[day][int(list_x[k] - 1)]) + list_beta[day][k])))

        object_stand.append(x)
    return object_stand


objectnov17 = standard_std(list_alpha, list_beta, 0)
objectnov23 = standard_std(list_alpha, list_beta, 1)
objectsep16 = standard_std(list_alpha, list_beta, 2)


f3=open("/home/avinash/Supernovae_Data/ASASSN14dq/Standards/PG0918/output_standard_mag", "a")

def write_data(std_name, mag_standard):
    f3.write(str(std_name) + "\n" + "\n")
    for j in range(0, star_count):
        if mag_standard == PG0231 or mag_standard == PG0918 or mag_standard == PG0942:
            f3.write(str("%3.0f" % star_id[j]) + "              " +
            str("%8.6f" % mag_standard[j][0]) + "              " +
            str("%8.6f" % mag_standard[j][1]) + "              " +
            str("%8.6f" % mag_standard[j][2]) + "               " +
            str("%8.6f" % mag_standard[j][3]) + "            " +
            str("%8.6f" % mag_standard[j][4]) + "               " +
            str("%8.6f" % mag_standard[j][5]) + "\n")
        else:
            f3.write(str("%3.0f" % star_id[j]) + "              " +
            str("%8.6f" % mag_standard[j][5]) + "              " +
            str("%8.6f" % mag_standard[j][0]) + "              " +
            str("%8.6f" % mag_standard[j][1]) + "              " +
            str("%8.6f" % mag_standard[j][2]) + "               " +
            str("%8.6f" % mag_standard[j][3]) + "           " +
            str("%8.6f" % mag_standard[j][4]) + "\n")

    f3.write("\n")

#write_data('PG0231-slopes - nov23', objectnov23)
#write_data('PG0231-slopes - nov17', objectnov17)
#write_data('PG0231-standard', PG0231)


def write_diff(std_name, standard, mag_standard):
    f3.write(str(std_name) + "\n" + "\n")
    for j in range(0, star_count):
        f3.write(str("%3.0f" % star_id[j]) + "\t\t" +
                 str("%8.6f" % (standard[j][0] - mag_standard[j][5])) + "\t" +
                 str("%8.6f" % (standard[j][0] - mag_standard[j][6])) + "\t" +
                 str("%8.6f" % (standard[j][1] - mag_standard[j][0])) + "\t" +
                 str("%8.6f" % (standard[j][2] - mag_standard[j][1])) + "\t" +
                 str("%8.6f" % (standard[j][3] - mag_standard[j][2])) + "\t" +
                 str("%8.6f" % (standard[j][4] - mag_standard[j][3])) + "\t" +
                 str("%8.6f" % (standard[j][5] - mag_standard[j][4])) + "\n")
    f3.write("\n")

write_diff('PG0231-(Delta-m)-nov23', PG0231, objectnov23)
write_diff('PG0231-(Delta-m)-nov17', PG0231, objectnov17)
write_diff('PG0231-(Delta-m)-sep16', PG0231, objectsep16)


f3.close()
"""
