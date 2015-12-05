
import datetime
import unittest

from eobs import Eobs

tg = Eobs(r'c:\m\ghs\data\e-obs\tg_0.25deg_reg_v12.0.nc')

# Paris coordinates and corresponding indices obtained from manual lookup in
# version 12 of the "tg" netCDF:
lat = 48.66
lon = 2.35
lat_idx = 93
lon_idx = 171


class TestEobs(unittest.TestCase):

    def test_lat2index(self):
        lat_index = tg.lat2index(lat)
        self.assertEqual(lat_index, lat_idx)

        # Check that the latitude referred to by this index represents the
        # grid square that the original latitude falls within
        lat2 = tg.variables['latitude'][lat_index]
        self.assertTrue(lat2 - tg.gridsize <= lat <= lat2 + tg.gridsize)

    def test_lon2index(self):
        lon_index = tg.lon2index(lon)
        self.assertEqual(lon_index, lon_idx)

        # Check that the latitude referred to by this index represents the
        # grid square that original latitude falls within
        lon2 = tg.variables['longitude'][lon_index]
        self.assertTrue(lon2 - tg.gridsize <= lon <= lon2 + tg.gridsize)

    def test_getitem_with_date_lat_lon(self):
        t = tg[datetime.datetime(1950, 1, 1), lat, lon]
        # Using lat-lon indices that have been looked up
        self.assertEqual(t, tg.variable[0, lat_idx, lon_idx])

