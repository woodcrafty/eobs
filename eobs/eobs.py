"""
Class to handle the reading in of E-OBS ENSEMBLES temperature and precipitation
data from CSV files based on the ELUM-HWSD 1km UK grid.
"""

import datetime
import os

import netCDF4 as cdf


class EobsError(Exception):
    pass


class Eobs(object):
    """
    Class to read in and access data from the E-OBS ENSEMBLES daily gridded
    observational climate dataset for Europe, which can be found here:

    http://www.ecad.eu/download/ensembles/download.php

    Note: This class only handles the 0.25 degree grid version of the E-OBS
    dataset and only handles 3 variables: mean temperature (tg),
    precipitation (rr) and elevation (elev).
    """

    def __init__(self, directory, version):
        """
        Return a new Eobs object.

        :param directory: string giving the directory that contains the
            E-OBS netCDF files.
        :param version: string giving the E-OBS dataset version number, e.g.
            '12.0'.
        """
        self.directory = directory
        tmean_filename = 'tg_0.25deg_reg_v{0}.nc'.format(version)
        rain_filename = 'rr_0.25deg_reg_v{0}.nc'.format(version)
        elev_filename = 'elev_0.25deg_reg_v{0}.nc'.format(version)

        self.tmean = cdf.Dataset(os.path.join(directory, tmean_filename))
        self.rain = cdf.Dataset(os.path.join(directory, rain_filename))
        self.elev = cdf.Dataset(os.path.join(directory, elev_filename))

        # Get the start and end date of the time series (as datetime objects):
        self.startdate = cdf.num2date(
            self.tmean.variables['time'][0],
            units=self.tmean.variables['time'].units,
            calendar=self.tmean.variables['time'].calendar,
        )
        self.enddate = cdf.num2date(
            self.tmean.variables['time'][-1],
            units=self.tmean.variables['time'].units,
            calendar=self.tmean.variables['time'].calendar,
        )
        # Calc number of days in time sereies
        self.ndays = self.enddate - self.startdate + datetime.timedelta(1).days
        assert(self.ndays == len(tmean.variables['time']))

        # Obtain the minimum latitude and longitude and grid size
        self.minlat = self.tmean.variables['latitude'][0]
        self.minlon = self.tmean.variables['longitude'][0]
        self.gridsize = (self.tmean.variables['latitude'][0]
                         - self.tmean.variables['latitude'][1]

    def get(self, lat, lon, date):
        lat_index = int((lat - self.minlat) / self.gridsize)
        lon_index = int((lon - self.minlon) / self.gridsize)
        t = date - self.startdate
        tmean = self.tmean.variables['tg'][t.days, lat_index, lon_index]
        rain = self.rain.variables['rr'][t.days, lat_index, lon_index]
        if tmean == '--':
            tmean = None
        if rain == '--':
            rain = None
        return tmean, rain

    def get_daterange(self, lat, lon, startdate, enddate):
        lat_index = int((lat - self.minlat) / self.gridsize)
        lon_index = int((lon - self.minlon) / self.gridsize)
        t1 = startdate - self.startdate
        t2 = enddate - self.startdate + datetime.timedelta(1)
        tmean = (
            self.tmean.variables['tg'][t1.days:t2.days, lat_index, lon_index])
        rain = (
            self.rain.variables['rr'][t1.days:t2.days, lat_index, lon_index])
        tmean = [None if t == '--' else t for t in tmean]
        rain = [None if r == '--' else r for r in rain]

        return tmean, rain

    def get_elev(self, lat, lon):
        lat_index = int((lat - self.minlat) / self.gridsize)
        lon_index = int((lon - self.minlon) / self.gridsize)
        elev = self.elev.variables['elevation'][lat_index, lon_index]
        if elev == '--':
            elev = None
        return elev


def test():
    eobs = Eobs(directory=r'c:\m\ghs\data\e-obs')

    print('get():')
    tmean, rain = eobs.get(lon=-3.0, lat=55.50, date=datetime.date(2000, 6, 1))
    print('tmean={0}, rain={1}'.format(tmean, rain))

    print('\nget_elev():')
    elev = eobs.get_elev(lon=-3.0, lat=55.50)
    print('elev={0}'.format(elev))

    print('\nget_daterange():')
    startdate = datetime.date(2000, 6, 1)
    enddate = datetime.date(2000, 6, 30)
    tmeans, rains = eobs.get_daterange(
        lon=-3.0, lat=55.50, startdate=startdate, enddate=enddate)
    for tmean, rain in zip(tmeans, rains):
        print('tmean={0}, rain={1}'.format(tmean, rain))

