#!/usr/bin/env python
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxx-------------------------OBSERVATION PLANNING----------------------xxxxxxxxxxxxxxxxxxxxxxxxx #
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #


# ------------------------------------------------------------------------------------------------------------------- #
# Import Required Libraries
# ------------------------------------------------------------------------------------------------------------------- #
import math
import ephem
import easygui
import datetime
import numpy as np
import astropy.units as u
from astropy.time import Time
from matplotlib import pyplot as plt
from astropy.coordinates import Angle
from matplotlib.dates import DateFormatter, HourLocator
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Observatory & Telescope Details
# ------------------------------------------------------------------------------------------------------------------- #
OBS_NAME = "Indian Astronomical Observatory, Hanle"
OBS_LONG = '78:57:51'
OBS_LAT = '32:46:46'
OBS_ALT = 4486
OBS_TIMEZONE = +5.5
telescope_horizon = 25
telescope_zenith = 85
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Functions For Error Handling
# ------------------------------------------------------------------------------------------------------------------- #

def remove_empty_values(python_list):
    """
    Args:
        python_list : Python list from which empty values are to be removed
    Returns:
        None
    """
    while True:
        try:
            python_list.remove('')
        except ValueError:
            break

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Defaults Used In Plotting Trajectories
# ------------------------------------------------------------------------------------------------------------------- #
time_offset = 0
date_obs = str(datetime.date.today())
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Manual Setup - GUI Code
# -------------------------------------------------------------------------------------------------------------------
box_msg = "Enter 'Name, RA, DEC' of objects to be plotted"
box_title = "Details of objects"
field_names = ['Object 1', 'Object 2', 'Object 3', 'Object 4', 'Object 5', '6', '7']
field_values = ['SN 2018gj 16:32:02.40 +78:12:41:13', '2017gkk 09:13:44.74 +76:28:44.10',
                '2017hpa 04:39:50.75 +07:03:54.90', '2016gfy 07:26:43.67 85:45:51.7',
                '2017ijx 10:51:53.02 +51:00:28.90', 'SN 2018ie 10:54:01.06 -16:01:21.40',
                'SN 2018gv 08:05:34.61 -11:26:16.30']

list_values = easygui.multenterbox(msg=str(box_msg), title=str(box_title), fields=field_names, values=field_values)

if list_values is not None:
    remove_empty_values(list_values)
else:
    exit()

while len(list_values) == 0:
    err_msg = box_msg + '\n\n Error: Aleast 1 Object required for plotting!!'
    list_values = easygui.multenterbox(msg=str(err_msg), title=str(box_title), fields=field_names, values=list_values)
    remove_empty_values(list_values)

choice_utc = easygui.boolbox(msg='Plot Trajectories W.r.t UTC?', title='UTC Or Local Time?', choices=['Yes', 'No'])
setup_manual = easygui.boolbox(msg='Manually Enter Date?', title='Manual Or Current Date?', choices=['Yes', 'No'])

if setup_manual:
    date_obs = easygui.enterbox(msg='Enter The Date Of Observation!', title='Date Of Observation',
                                    default=str(datetime.date.today()))

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Declaring Object 'telescope'
# ------------------------------------------------------------------------------------------------------------------- #
telescope = ephem.Observer()
telescope.lon = OBS_LONG
telescope.lat = OBS_LAT
telescope.elevation = OBS_ALT
telescope.pressure = 0
telescope.epoch = ephem.J2000
telescope.date = (Time(date_obs) + 1 * u.day - abs(OBS_TIMEZONE) * u.hour).utc.datetime
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculation Of Times Of Local Sunset & Sunrise
# ------------------------------------------------------------------------------------------------------------------- #
telescope.horizon = '-0:34'
time_utc_start = telescope.previous_setting(ephem.Sun(), use_center=True)
time_utc_end = telescope.next_rising(ephem.Sun(), use_center=True)

datetime_utc_sunset = Time(datetime.datetime.strptime(str(time_utc_start).split('.')[0], '%Y/%m/%d %H:%M:%S'))
datetime_utc_sunrise = Time(datetime.datetime.strptime(str(time_utc_end).split('.')[0], '%Y/%m/%d %H:%M:%S'))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Calculation Of Times Of Nautical Twilight [Elevation Of Sun = -12 Degrees]
# ------------------------------------------------------------------------------------------------------------------- #
telescope.horizon = '-13'
time_utc_twil_set = telescope.previous_setting(ephem.Sun(), use_center=True)
time_utc_twil_rise = telescope.next_rising(ephem.Sun(), use_center=True)

datetime_utc_night_start = Time(datetime.datetime.strptime(str(time_utc_twil_set).split('.')[0], '%Y/%m/%d %H:%M:%S'))
datetime_utc_night_end = Time(datetime.datetime.strptime(str(time_utc_twil_rise).split('.')[0], '%Y/%m/%d %H:%M:%S'))
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Determining Time Intervals
# ------------------------------------------------------------------------------------------------------------------- #
plot_duration = (datetime_utc_sunrise.utc.datetime - datetime_utc_sunset.utc.datetime).total_seconds() / 3600
datetime_utc_intervals = datetime_utc_sunset + np.linspace(time_offset, time_offset + plot_duration, 100) * u.hour
datetime_local_intervals = datetime_utc_intervals + OBS_TIMEZONE * u.hour
datetime_sep_intervals = datetime_utc_sunset + np.linspace(time_offset, time_offset + plot_duration, 15)[1:-1] * u.hour
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Colors Used In Plotting
# ------------------------------------------------------------------------------------------------------------------- #
colors = ['blue', 'green', 'red', 'saddlebrown', 'fuchsia', 'teal', 'black']
object_count = 0
# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Class 'ObjectToObs' For Declaring Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #

class ObjectToObs:
    def __init__(self, object_name, object_ra, object_dec):
        self.name = str(object_name)
        self.object = ephem.FixedBody()
        self.object._epoch = ephem.J2000
        self.object._ra = str(object_ra)
        self.object._dec = str(object_dec)
        self.list_alt = []

    def get_altitude(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        object_alt = Angle(str(self.object.alt) + ' degrees').degree
        return object_alt

    def get_moon_sep(self, time_obs):
        global telescope
        telescope.date = str(time_obs)
        self.object.compute(telescope)
        moon_pos = ephem.Moon(str(time_obs))
        angle_ephem = ephem.separation(self.object, moon_pos)
        angle_sep = int(Angle(str(angle_ephem) + ' degrees').degree)
        return angle_sep

    def plot_alt(self, utc=True):
        for time_obs in list(datetime_utc_intervals.value):
            self.list_alt.append(self.get_altitude(time_obs))
        if utc:
            self.plot_alt_utc()
        else:
            self.plot_alt_local()

    def plot_alt_utc(self):
        global object_count
        plt.plot(list(datetime_utc_intervals.value), self.list_alt, label=self.name, color=colors[object_count])

        for time_obs in datetime_sep_intervals:
            if int(self.get_altitude(time_obs)) > 0:
                plt.text(time_obs.value, self.get_altitude(time_obs), str(self.get_moon_sep(time_obs)))
        object_count += 1

    def plot_alt_local(self):
        global object_count
        plt.plot(list(datetime_local_intervals.value), self.list_alt, label=self.name, color=colors[object_count])

        for time_obs in datetime_sep_intervals:
            if int(self.get_altitude(time_obs)) > 0:
                local_time = time_obs + OBS_TIMEZONE * u.hour
                plt.text(local_time.value, self.get_altitude(time_obs), str(self.get_moon_sep(time_obs)))
        object_count += 1

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Setting Plot Parameters
# ------------------------------------------------------------------------------------------------------------------- #

def set_plot_params(utc=True):
    """
    Sets plot parameters for plotting the trajectory of objects in sky.
    Args:
        utc : Boolean value to be determine whether UTC or Local Time is to be used for plotting
    Returns:
        None
    """
    def sign(value):
        return (value > 0) - (value < 0)

    lat_deg = "%7.4f" % Angle(OBS_LAT + ' degrees').degree
    long_deg = "%7.4f" % Angle(OBS_LONG + ' degrees').degree

    text_ns = 'N'
    text_ew = 'E'

    if not sign(lat_deg):
        text_ns = 'S'
    if not sign(long_deg):
        text_ew = 'W'

    ax = plt.gca()
    fig = plt.gcf()

    degree_sign = "$^\circ$"
    text_name = str(OBS_NAME) + " [+" + str(OBS_TIMEZONE) + "h]" + "\n"
    text_lat = "Latitude : " + str(lat_deg) + str(degree_sign) + str(text_ns)
    text_long = ", Longitude : " + str(long_deg) + str(degree_sign) + str(text_ew)
    text_alt = ", Altitude : " + str(OBS_ALT) + 'm'
    display_text = str(text_name) + str(text_lat) + str(text_long) + str(text_alt) + "\n"
    ax.set_title(str(display_text), fontsize=20)

    if utc:
        ax.set_xlabel('UTC Time' + '\n\n' + 'DATE : ' + str(date_obs), fontsize=20)
        datetime_sunset = datetime_utc_sunset
        datetime_sunrise = datetime_utc_sunrise
        datetime_night_start = datetime_utc_night_start
        datetime_night_end = datetime_utc_night_end
        datetime_current = Time.now()
    else:
        ax.set_xlabel('Local Time (IST)' + '\n\n' + 'DATE : ' + str(date_obs), fontsize=20)
        datetime_sunset = datetime_utc_sunset + + OBS_TIMEZONE * u.hour
        datetime_sunrise = datetime_utc_sunrise + OBS_TIMEZONE * u.hour
        datetime_night_start = datetime_utc_night_start + OBS_TIMEZONE * u.hour
        datetime_night_end = datetime_utc_night_end + OBS_TIMEZONE * u.hour
        datetime_current = Time.now() + OBS_TIMEZONE * u.hour

    ax.set_ylabel('Elevation (In Degrees)', fontsize=20)
    ax.set_ylim(0, 90, 10)
    ax.set_xlim(datetime_sunset.value, datetime_sunrise.value)
    ax.grid(True)
    ax.legend(shadow=True, loc=1)

    hours = HourLocator()
    hourfmt = DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_locator(hours)
    ax.xaxis.set_major_formatter(hourfmt)
    ax.tick_params(labelsize=12)
    ax.autoscale_view()
    fig.autofmt_xdate()

    if datetime_utc_sunset.value < datetime_current.utc.datetime < datetime_utc_sunrise.value:
        ax.axvline(x=datetime_current.value, linestyle='--', color='k')
        ax.text(datetime_current.value, (ax.get_ybound()[0] + ax.get_ybound()[-1]) / 2, 'Current Time',
                rotation=-90, color='k')

    ax.axvline(x=datetime_night_start.value, linestyle='--', color='k')
    ax.axvline(x=datetime_night_end.value, linestyle='--', color='k')
    ax.text(datetime_night_start.value, ax.get_ybound()[0] + 13, 'Nautical Dusk', rotation=-90, color='k')
    ax.text(datetime_night_end.value, ax.get_ybound()[0] + 13, 'Nautical Dawn', rotation=-90, color='k')

    ax.text(datetime_sunset.value, ax.get_ybound()[1] + 6, 'Sunset', rotation=+50, color='k', fontsize=12)
    ax.text(datetime_sunrise.value,  ax.get_ybound()[1] + 6, 'Sunrise', rotation=+50, color='k', fontsize=12)
    ax.text(datetime_night_start.value, ax.get_ybound()[1] + 7, 'Twilight', rotation=+50, color='k', fontsize=12)
    ax.text(datetime_night_end.value, ax.get_ybound()[1] + 7, 'Twilight', rotation=+50, color='k', fontsize=12)

    ax.set_facecolor('lightgray')
    ax.text(ax.get_xbound()[0], telescope_horizon + 1, 'Telescope Horizon')
    ax.text(ax.get_xbound()[0], telescope_zenith - 2, 'Telescope Zenith')
    ax.fill_between(ax.get_xbound(), telescope_horizon - 0.5, telescope_horizon + 0.5, facecolor='royalblue')
    ax.fill_between(ax.get_xbound(), telescope_zenith - 0.5, telescope_zenith + 0.5, facecolor='royalblue')
    ax.fill_between(ax.get_xbound(), telescope_horizon + 0.5, telescope_zenith - 0.5, facecolor='white')

    ax.fill_between([datetime_sunset.value, datetime_night_start.value], telescope_horizon + 0.5,
                    telescope_zenith - 0.5, facecolor='paleturquoise')
    ax.fill_between([datetime_night_end.value, datetime_sunrise.value], telescope_horizon + 0.5,
                    telescope_zenith - 0.5, facecolor='paleturquoise')

    list_secz = []
    for altitude in ax.get_yticks():
        if (1 / math.cos(math.radians(90 - altitude))) < 10:
            list_secz.append("%5.2f" % (1 / math.cos(math.radians(90 - altitude))))
        else:
            list_secz.append("NaN")

    ax2 = ax.twinx()
    ax2.set_ylim(0, 90)
    ax2.set_yticks(ax.get_yticks())
    ax2.set_yticklabels(list_secz)
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.set_ylabel('Airmass', fontsize=20)

# ------------------------------------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------------------------- #
# Plots The Trajectories Of Objects To Be Observed
# ------------------------------------------------------------------------------------------------------------------- #
for index, value in enumerate(list_values):
    if len(value.split()) >= 3:
        ObjectToObs(object_name=value.split()[-3], object_ra=value.split()[-2],
                    object_dec=value.split()[-1]).plot_alt(utc=choice_utc)
    elif len(value.split()) == 2:
        ObjectToObs(object_name='Object ' + str(int(index) + 1), object_ra=value.split()[-2],
                    object_dec=value.split()[-1]).plot_alt(utc=choice_utc)
    else:
        print ("Error : Both RA & DEC For Object {} Need To Be Specified".format(str(int(index) + 1)))

set_plot_params(utc=choice_utc)
plt.show()
# ------------------------------------------------------------------------------------------------------------------- #
