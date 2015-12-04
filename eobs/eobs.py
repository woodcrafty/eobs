"""
Class to handle the reading in of E-OBS ENSEMBLES temperature and precipitation
data from CSV files based on the ELUM-HWSD 1km UK grid.
"""

import datetime

import netCDF4 as cdf
from netCDF4 import date2index


class EobsError(Exception):
    pass


class Eobs(object):
    def __init__(self, filename):

        self.dataset = cdf.Dataset(filename)
        self.variables = self.dataset.variables  # For brevity

        # Get the variable of interest (e.g. mean temperature,
        # precipitation etc) from the variables OrderedDict
        var_name = [k for k in self.dataset.variables][3]
        # The main variable of the dataset (e.g. mean temperature)
        self.variable = self.variables[var_name]

        # Get the start and end date of the time series (as datetime objects):
        self.startdate = cdf.num2date(
            self.variables['time'][0],
            units=self.variables['time'].units,
            calendar=self.variables['time'].calendar,
        )
        self.enddate = cdf.num2date(
            self.variables['time'][-1],
            units=self.variables['time'].units,
            calendar=self.variables['time'].calendar,
        )

        # Calc number of days in time series
        self.ndays = (
            self.enddate - self.startdate + datetime.timedelta(1)).days
        assert(self.ndays == len(self.variables['time']))

        # Grid size in degrees
        self.gridsize = abs(
            (self.variables['latitude'][0] - self.variables['latitude'][1]))
        # Calculate the minimum latitude and longitude of the lower left (SW)
        # corner of the grid
        self.minlat = self.variables['latitude'][0] - (self.gridsize / 2.0)
        self.minlon = self.variables['longitude'][0] - (self.gridsize / 2.0)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        self.f.close()

    def __getitem__(self, *args):
        """
        Example usage:
        >>> import netCDF4 as cdf
        >>> # Retrieve the value on a single day for a single grid cell
        >>> latitude = 51.5
        >>> longitude = 1.2
        >>> tg = Eobs('filename.nc')
        >>> t = tg[datetime.date(2001, 1, 1), latitude, longitude]
        >>> # Retrieve the values for a whole year for a single grid cell
        >>> start_date = datetime.date(2001, 1, 1)
        >>> end_date = datetime.date(2001, 12, 31)
        >>> t_2001 = tg[start_date:end_date, latitude, longitude]
        >>> # Retrieve the values for a whole year for a 1 x 1 degree sub grid
        >>> ll_lat = 51.0  # Latitude of lower left corner
        >>> ll_lon = 1.0   # Longitude of lower left corner
        >>> ur_lat = 52.0  # Latitude of upper right corner
        >>> ur_lon = 2.0   # Longitude of upper right corner
        >>> t_2001 = tg[start_date:end_date, ul_lat:urlat, ll_lon:ur_lon]

        :param date: datetime object or a slice with a start date and end date
            datetime object and an optional step value with units of days.
        :param latitude: latitude of the grid point or a slice specifying the
            latitude of the lower left (SW) and upper right (NE) corners of a
            rectangle. Note that step is not supported in the latitude slice.
        :param longitude: longitude of the grid point or a slice specifying the
            longitude of the lower left (SW) and upper right (NE) corners of a
            rectangle. Note that step is not supported in the longitude slice.
        """
        if len(args[0]) > 3:
            raise IndexError('too many indices for array')

        date = args[0][0]
        if isinstance(date, slice):
            d1 = date2index(
                date.start,
                nctime=self.variables['time'],
                calendar=self.variables['time'].calendar,
            )
            d2 = date2index(
                date.stop,
                nctime=self.variables['time'],
                calendar=self.variables['time'].calendar,
            )
            date_idx = slice(d1, d2, date.step)
        else:
            date_idx = date2index(
                date,
                nctime=self.variables['time'],
                calendar=self.variables['time'].calendar,
            )

        try:
            lat = args[0][1]
        except IndexError:
            lat_idx = slice(None)
        else:
            # Note that the slice.step is not supported for latitude
            if isinstance(lat, slice):
                lat_idx = slice(
                    int((lat.start - self.minlat) / self.gridsize),
                    int((lat.stop - self.minlat) / self.gridsize)
                )
            else:
                lat_idx = int((lat - self.minlat) / self.gridsize)

        try:
            lon = args[0][2]
        except IndexError:
            lon_idx = slice(None)
        else:
            # Note that the slice.step is not supported for longitude
            if isinstance(lon, slice):
                lon_idx = slice(
                    int((lon.start - self.minlon) / self.gridsize),
                    int((lon.stop - self.minlon) / self.gridsize)
                )
            else:
                lon_idx = int((lon - self.minlon) / self.gridsize)

        return self.variable[date_idx, lat_idx, lon_idx]


def test():
    eobs = Eobs(r'c:\m\ghs\data\e-obs\tg_0.25deg_reg_v12.0.nc')
    lat = 48.86
    lon = 2.35
    t = eobs[datetime.datetime(2001, 1, 1), lat, lon]
    print(t)
    print('---------------')

    # Retrieve the values for a week year for a single grid cell
    start_date = datetime.datetime(2001, 1, 1)
    end_date = datetime.datetime(2001, 1, 8)
    t_week = eobs[start_date:end_date, lat, lon]
    print(t_week)
    print('---------------')

    # Retrieve the values for a whole year for a 1 x 1 degree sub grid
    ll_lat = 48.0  # Latitude of lower left corner
    ll_lon = 2.0   # Longitude of lower left corner
    ur_lat = 49.0  # Latitude of upper right corner
    ur_lon = 3.0   # Longitude of upper right corner
    t_grid = eobs[start_date, ll_lat:ur_lat, ll_lon:ur_lon]
    print(t_grid)


test()
#
# class Eobs(object):
#     """
#     Class to read in and access data from the E-OBS ENSEMBLES daily gridded
#     observational climate dataset for Europe, which can be found here:
#
#     http://www.ecad.eu/download/ensembles/download.php
#
#     Note: This class only handles the 0.25 degree grid version of the E-OBS
#     dataset and only handles 3 variables: mean temperature (tg),
#     precipitation (rr) and elevation (elev).
#     """
#
#     def __init__(self, directory, version):
#         """
#         Return a new Eobs object.
#
#         :param directory: string giving the directory that contains the
#             E-OBS netCDF files.
#         :param version: string giving the E-OBS dataset version number, e.g.
#             '12.0'.
#         """
#         self.directory = directory
#         tmean_filename = 'tg_0.25deg_reg_v{0}.nc'.format(version)
#         rain_filename = 'rr_0.25deg_reg_v{0}.nc'.format(version)
#         elev_filename = 'elev_0.25deg_reg_v{0}.nc'.format(version)
#
#         self.tmean = cdf.Dataset(os.path.join(directory, tmean_filename))
#         self.rain = cdf.Dataset(os.path.join(directory, rain_filename))
#         self.elev = cdf.Dataset(os.path.join(directory, elev_filename))
#
#         # Get the start and end date of the time series (as datetime objects):
#         self.startdate = cdf.num2date(
#             self.tmean.variables['time'][0],
#             units=self.tmean.variables['time'].units,
#             calendar=self.tmean.variables['time'].calendar,
#         )
#         self.enddate = cdf.num2date(
#             self.tmean.variables['time'][-1],
#             units=self.tmean.variables['time'].units,
#             calendar=self.tmean.variables['time'].calendar,
#         )
#         # Calc number of days in time sereies
#         self.ndays = self.enddate - self.startdate + datetime.timedelta(1).days
#         assert(self.ndays == len(tmean.variables['time']))
#
#         # Obtain the minimum latitude and longitude and grid size
#         self.minlat = self.tmean.variables['latitude'][0]
#         self.minlon = self.tmean.variables['longitude'][0]
#         self.gridsize = (self.tmean.variables['latitude'][0]
#                          - self.tmean.variables['latitude'][1])
#
#     def get(self, lat, lon, date):
#         lat_index = int((lat - self.minlat) / self.gridsize)
#         lon_index = int((lon - self.minlon) / self.gridsize)
#         t = date - self.startdate
#         tmean = self.tmean.variables['tg'][t.days, lat_index, lon_index]
#         rain = self.rain.variables['rr'][t.days, lat_index, lon_index]
#         if tmean == '--':
#             tmean = None
#         if rain == '--':
#             rain = None
#         return tmean, rain
#
#     def get_daterange(self, lat, lon, startdate, enddate):
#         lat_index = int((lat - self.minlat) / self.gridsize)
#         lon_index = int((lon - self.minlon) / self.gridsize)
#         t1 = startdate - self.startdate
#         t2 = enddate - self.startdate + datetime.timedelta(1)
#         tmean = (
#             self.tmean.variables['tg'][t1.days:t2.days, lat_index, lon_index])
#         rain = (
#             self.rain.variables['rr'][t1.days:t2.days, lat_index, lon_index])
#         tmean = [None if t == '--' else t for t in tmean]
#         rain = [None if r == '--' else r for r in rain]
#
#         return tmean, rain
#
#     def get_elev(self, lat, lon):
#         lat_index = int((lat - self.minlat) / self.gridsize)
#         lon_index = int((lon - self.minlon) / self.gridsize)
#         elev = self.elev.variables['elevation'][lat_index, lon_index]
#         if elev == '--':
#             elev = None
#         return elev
#
#
# def test():
#     eobs = Eobs(directory=r'c:\m\ghs\data\e-obs')
#
#     print('get():')
#     tmean, rain = eobs.get(lon=-3.0, lat=55.50, date=datetime.date(2000, 6, 1))
#     print('tmean={0}, rain={1}'.format(tmean, rain))
#
#     print('\nget_elev():')
#     elev = eobs.get_elev(lon=-3.0, lat=55.50)
#     print('elev={0}'.format(elev))
#
#     print('\nget_daterange():')
#     startdate = datetime.date(2000, 6, 1)
#     enddate = datetime.date(2000, 6, 30)
#     tmeans, rains = eobs.get_daterange(
#         lon=-3.0, lat=55.50, startdate=startdate, enddate=enddate)
#     for tmean, rain in zip(tmeans, rains):
#         print('tmean={0}, rain={1}'.format(tmean, rain))
#
