#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxx----------------------READ DATA FOR SN 2018hna----------------xxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import json
import glob
import shutil
import numpy as np
import pandas as pd
from astropy.table import Table
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2018hna'
epoch = 2400000.5
precision = 4
select_ztfcols = ['utc', 'mjd', 'fid', 'magpsf', 'sigmapsf']
select_asassncols = ['HJD', 'UT Date', 'mag', 'mag_err', 'Filter']
names_ztfcols = ['Date', 'MJD', 'FILTER', 'FMAG', 'FERR']
names_gaiacols = ['Date', 'JD', 'FMAG']
names_asassncols = ['JD', 'Date', 'FMAG', 'FERR', 'FILTER']

file_ztf = 'ZTF.json'
file_gaia = 'Gaia.csv'
file_asassn = 'ASASSN.csv'
DIR_SNe = '/home/avinash/Dropbox/SNData/IIP_Data/LC_Data/'
dict_fmt = {'Date': 's', 'FILTER': 's', 'RA': '0.6f', 'DEC': '0.6f', 'FMAG': '0.3f', 'FERR': '0.3f', 'JD': '0.4f'}
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read JSON File From ZTF
# ------------------------------------------------------------------------------------------------------------------- #
data_raw = json.loads(open(file_ztf).read())
# for key in data_raw:
#     print key
# print data_raw['candidates']

data_ztf = pd.DataFrame(data_raw['candidates'])
data_ztf = data_ztf[select_ztfcols]
data_ztf.columns = names_ztfcols
data_ztf['JD'] = (data_ztf['MJD'] + epoch).round(precision)
data_ztf = data_ztf.dropna(how='any').drop('MJD', axis=1)
data_ztf['Date'] = data_ztf['Date'].apply(lambda x: x.replace(' ', 'T'))
data_ztf['FILTER'] = data_ztf['FILTER'].apply(lambda x: 'ZTFg' if x == 1 else 'ZTFr')

list_reject = [2458440.9512, 2458534.7677, 2458426.9636, 2458427.0097, 2458543.7881, 2458538.7848,
               2458538.8224, 2458543.8655, 2458538.7204, 2458429.9927, 2458472.0528, 2458461.9921,
               2458486.9897, 2458503.0354, 2458511.9548, 2458522.9869, 2458514.9719, 2458432.9883]
data_ztf = data_ztf[~data_ztf['JD'].isin(list_reject)]

data_ztf2 = data_ztf.copy()
data_ztf2['JD'] = data_ztf2['JD'].round(1)
data_ztf2 = data_ztf2.set_index(['JD', 'FILTER'])
data_ztf2 = data_ztf2.groupby(['JD', 'FILTER']).mean().reset_index()

# data_ztftab = Table.from_pandas(data_ztf)
# data_ztftab.write('2018hna_ZTF.dat', format='ascii.fixed_width', delimiter=' ',
#                   formats={x: dict_fmt[x] for x in data_ztf.columns}, overwrite=True)

data_ztftab2 = Table.from_pandas(data_ztf2)
data_ztftab2.write('2018hna_ZTF.dat', format='ascii.fixed_width', delimiter=' ',
                   formats={x: dict_fmt[x] for x in data_ztf2.columns}, overwrite=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read CSV File From Gaia
# ------------------------------------------------------------------------------------------------------------------- #
data_gaia = pd.read_csv(file_gaia, comment='#')
data_gaia.columns = names_gaiacols
data_gaia = data_gaia.replace('null', np.nan).dropna(how='any')
data_gaia['FILTER'] = 'Gaiag'

data_gaiatab = Table.from_pandas(data_gaia)
data_gaiatab.write('2018hna_Gaia.dat', format='ascii.fixed_width', delimiter=' ',
                   formats={x: dict_fmt[x] for x in data_gaia.columns}, overwrite=True)
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read CSV File From ASASSN
# ------------------------------------------------------------------------------------------------------------------- #
data_asassn = pd.read_csv(file_asassn, comment='#')
data_asassn = data_asassn[select_asassncols]
data_asassn.columns = names_asassncols
data_asassn = data_asassn[~data_asassn['FMAG'].str.startswith('>')]
data_asassn['FMAG'] = data_asassn['FMAG'].astype('float64')
data_asassn = data_asassn[~data_asassn['JD'].round(4).isin([2458474.0592])]
data_asassn = data_asassn[data_asassn['FMAG'] != 99.990]

data_asassn2 = data_asassn.copy()
data_asassn2['JD'] = data_asassn2['JD'].round(1)
data_asassn2 = data_asassn2.set_index(['JD', 'FILTER'])
data_asassn2 = data_asassn2.groupby(['JD', 'FILTER']).mean().reset_index()

# data_asassntab = Table.from_pandas(data_asassn)
# data_asassntab.write('2018hna_ASASSN.dat', format='ascii.fixed_width', delimiter=' ',
#                    formats={x: dict_fmt[x] for x in data_asassn.columns}, overwrite=True)

data_asassntab2 = Table.from_pandas(data_asassn2)
data_asassntab2.write('2018hna_ASASSN.dat', format='ascii.fixed_width', delimiter=' ',
                      formats={x: dict_fmt[x] for x in data_asassn2.columns}, overwrite=True)
# ------------------------------------------------------------------------------------------------------------------- #

for file_name in glob.glob('2018hna_*.dat'):
    shutil.copy(file_name, DIR_SNe + file_name)
