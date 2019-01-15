#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx-------------------CALCULATE EXTINCTION COEFFICIENTS (R_Lambda)---------------xxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from scipy.interpolate import UnivariateSpline
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Cubic Spline Anchor Points For R_v = 3.1
# ------------------------------------------------------------------------------------------------------------------- #
R_data = [0, 0.265, 0.829, 2.688, 3.055, 3.806, 4.315, 6.265, 6.591]
wave_data = [float("inf"), 26500, 12200, 6000, 5470, 4670, 4110, 2700, 2600]
waveinv_data = [1./value for value in wave_data]

wave_array = np.linspace(2500, 70000, 100)
waveinv_array = [1./value for value in wave_array]

spline = UnivariateSpline(waveinv_data, R_data, k=3)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plot The Fit To The Extinction Curve
# ------------------------------------------------------------------------------------------------------------------- #
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111)

ax.grid()
ax.yaxis.set_major_locator(FixedLocator(np.arange(-1.0, 8.25, 1)))
ax.yaxis.set_minor_locator(FixedLocator(np.arange(-1.0, 8.25, 0.25)))
ax.xaxis.set_major_locator(FixedLocator(np.arange(0.0000, 0.0005, 0.0001)))
ax.xaxis.set_minor_locator(FixedLocator(np.arange(0.0000, 0.0005, 0.000025)))
ax.tick_params(which='both', direction='in', width=0.2, labelsize=9)
ax.plot(waveinv_data, R_data, 'b*', waveinv_array, spline(waveinv_array), 'r-')

axopp = ax.twiny()
axopp.set_xlim(ax.get_xlim())
axopp.set_xticks(ax.get_xticks())
list_labels = ["{:.1e}".format(1 / value) if value != 0 else float("inf") for value in ax.get_xticks()]
labels = axopp.set_xticklabels(list_labels)
axopp.tick_params(which='both', direction='in', width=0.2, labelsize=9)

plt.savefig("Plot_FitRv.eps", format='eps')
plt.show()
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculate Extinction Coefficients (R_lambda) At Central Wavelength Of Different Filters
# ------------------------------------------------------------------------------------------------------------------- #
name_filter = ['U', 'B', 'V', 'R', 'I', 'u', 'g', 'r', 'i', 'z']
wave_filter = [3700, 4200, 5300, 6000, 8050, 3655, 5105, 6480, 7105, 8640]
waveinv_filter = [1./value for value in wave_filter]

R_output = spline(waveinv_filter)
print(R_output)
# ------------------------------------------------------------------------------------------------------------------- #
