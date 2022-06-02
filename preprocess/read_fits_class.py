#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxx-----------------Read FITS file----------------xxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
from os import read
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Class to Read FITS files of Observed Data
# ------------------------------------------------------------------------------------------------------------------- #

class read_fits:

    def __init__(self, filename):
        self.filename = filename
        self.file = fits.open(filename)

    def file_info(self):
        return self.file.info()

    def file_headers(self):
        return self.file[0].header

    def file_data(self):
        return self.file[0].data

    def file_plot(self, scale=3, colormap='gray'):
        self.multiplier = scale
        self.colormap = colormap
        mu, med, std = sigma_clipped_stats(self.file_data())
        plt.imshow(self.file_data(), cmap=colormap, vmin=mu-scale*std, vmax=mu+scale*std)
        plt.colorbar()
        plt.show()

    def modify_header(self, header_name='COMMENT', header_value='None', ext=0):
        self.header_name = header_name
        self.header_value = header_value
        self.ext = ext
        fits.setval(self.filename, header_name, value=header_value, ext=ext)
        print(header_name+' is added/modifed with value '+str(header_value))


# ------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------------- #

if __name__ == '__main__':
    filename = 'awcs_cfb_2019-05-11-ZTF19aaskkq_g_8.fits'
    file = read_fits(filename)
    file.file_info()
    print(file.file_headers())
    print(file.file_data())
    file.file_plot(scale=2, colormap='gray')

# ------------------------------------------------------------------------------------------------------------------- #
