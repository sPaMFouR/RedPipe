#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxx----------------FORMAT THE LIGHT CURVE DATA IN ACCORDANCE WITH VOSA----------------xxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #

# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import pandas as pd
from jdcal import jd2gcal
from datetime import date
from astropy.coordinates import Angle
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Global Variables To Be Used In The Code
# ------------------------------------------------------------------------------------------------------------------- #
name_SN = '2016gfy'
epoch = 2400000.5
obj_RA = '07:26:43.67'
obj_DEC = '85:45:51.70'

input_file = 'OUTPUT_FinalSNMag'
DIR_PHOT = '/home/avinash/Supernovae_Data/2016gfy/Photometry/'
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Useful Functions
# ------------------------------------------------------------------------------------------------------------------- #

def jd_to_cald(julian_day):
    """
    Converts julian day into calendar day in string format.
    Args:
        julian_day  : Julian day value to be converted to calendar day
    Returns:
        cal_date    : Calendar date corresponding to input julian day
    """
    time_tuple = jd2gcal(epoch, julian_day - epoch)
    cal_date = date(*time_tuple[0:3]).strftime("%Y-%m-%d")

    return cal_date
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Read Input Data & Convert To VOSA Format
# ------------------------------------------------------------------------------------------------------------------- #

data_raw = pd.read_csv(DIR_PHOT + input_file, sep='\s+')
data_raw['Date'] = data_raw['JD'].apply(jd_to_cald)
master_df = pd.DataFrame()

for date, part_df in data_raw.groupby(['Date']):
    data_df = pd.DataFrame(index=range(0, len(part_df.index.values)))
    data_df['object'] = name_SN + '_' + date
    data_df['RA'] = Angle(obj_RA + ' hours') * 15
    data_df['DEC'] = Angle(obj_DEC + ' degrees')
    data_df['dis'] = 36.0
    data_df['Av'] = 0.0865
    data_df['filter'] = part_df['FILTER'].values
    data_df['filter'] = data_df['filter'].apply(lambda x: 'Generic/Bessell.' + x)
    data_df['flux'] = part_df['FMAG'].values
    data_df['error'] = part_df['FERR'].values
    data_df['pntopts'] = ''
    data_df['objopts'] = ''
    master_df = pd.concat([master_df, data_df], axis=0)

master_df = master_df.reset_index(drop=True)
master_df = master_df.to_csv('VOSA_2016gfy.dat', sep=' ', header=None, index=False)
# ------------------------------------------------------------------------------------------------------------------- #

