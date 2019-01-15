#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxx---------------------LEAST SQUARE MINIMIZATION FOR VALENTI MODEL----------------xxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from lmfit import Minimizer, Parameters, report_fit
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Load Observational Data
# ------------------------------------------------------------------------------------------------------------------- #
data1 = np.loadtxt('bol_05cf')
list_time_int = list(data1[:, 0])
list_bolflx = list(data1[:, 1])
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Constants Of Radioactive Decay
# ------------------------------------------------------------------------------------------------------------------- #
eni = 3.90e10
eco = 6.78e9
tauni = 757728.0
tauco = 9616320.0
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Function For Calculating Luminosity From Valenti Model Based On Initial Values
# ------------------------------------------------------------------------------------------------------------------- #

def calc_lph(param, list_time):
    mni = param['mni'].value
    taum = param['taum'].value
    time_int = [float(value / taum) for value in list_time]
    y = float(taum / (2 * tauni))
    s = float(((taum * (tauco - tauni))/(2 * tauco * tauni)))
    lph = []
    lph_err = []

    def A(z, y):
        return 2 * z * np.exp((-2 * z * y) + (z ** 2))

    def B(z, y, s):
        return 2 * z * np.exp((-2 * z * y)+(2 * z * s) + (z ** 2))

    for time in time_int:
        lph.append(1.989e33 * mni * np.exp(-time ** 2)*((eni - eco) * quad(A, 0, time, args=(y,))[0]
                                                        + eco*quad(B, 0, time, args=(y, s))[0]))
        lph_err.append(1.989e33 * mni * np.exp(-time ** 2)*((eni - eco) * quad(A, 0, time, args=(y,))[1]
                                                            + eco*quad(B, 0, time, args=(y, s))[1]))
    return [lph, lph_err]

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Defines Objective Function For Best Fit- Returns Numpy Array To Be Minimized
# ------------------------------------------------------------------------------------------------------------------- #

def func_calc(param, list_time, list_flux):
    lph = calc_lph(param, list_time)[0]
    delta = []
    delta_err = []

    for i in range(0, len(list_flux)):
        delta.append(float(lph[i]) - float(list_flux[i]))
        delta_err.append((float(lph[i]) - float(list_flux[i])) / float(list_flux[i]))

    return np.asarray(delta)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Create A Set Of Parameters
# ------------------------------------------------------------------------------------------------------------------- #
params = Parameters()
params.add('mni', value=0.6)
params.add('taum', value=2e6)
# ------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------------------------------------------------------------------------- #
# Minimize Using 'Least Squares' Technique
# ------------------------------------------------------------------------------------------------------------------- #
minner_object = Minimizer(func_calc, params, fcn_args=(list_time_int, list_bolflx))
kws = {'options': {'maxiter': 10}}
result_fit = minner_object.minimize()
list_lum_fit = list_bolflx + result_fit.residual
report_fit(result_fit)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots 'Observational_Data', 'Model_Init', 'Best_Fit_Valenti' W.R.T Time
# ------------------------------------------------------------------------------------------------------------------- #
plt.scatter(list_time_int, list_bolflx, label='Observation')
plt.plot(list_time_int, calc_lph(params, list_time_int)[0], label='Valenti-Initial Fit')
plt.plot(list_time_int, list_lum_fit, label='Valenti-Best Fit')
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Set Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #
plt.title('Valenti Model - Least Square Fitting', fontsize=20)
plt.xlabel('Time Since Explosion (In Seconds)', fontsize=20)
plt.ylabel('Bolometric Luminosity (In $10^{42}$ ergs)', fontsize=20)
plt.yticks(np.linspace(1e41, 2e43, 10))
plt.xticks(np.linspace(0.5e6, 5e6, 10))
plt.grid()
plt.legend()
plt.show()
# ------------------------------------------------------------------------------------------------------------------- #

# print result.chisqr
# print result.params
# print result1.chisqr
# print result1.params

# params = Parameters()
#
# list_mni = []
# for value in np.linspace(0.1, 0.4, 5):
#     list_mni.append(value)
#
# list_taum = []
# for value in np.linspace(8e4, 2e6, 5):
#     list_taum.append(value)
#
# list_chisqr = []
# for first_value in list_mni:
#     for second_value in list_taum:
#         print first_value
#         print second_value
#         params.add('mni', value=first_value)
#         params.add('taum', value=second_value)
#         result1 = minimize(func_calc, params, args=(list_t, list_bolflx), method='leastsq')
#         list_chisqr.append(result1.chisqr)
#
# print list_chisqr
