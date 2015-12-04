
import datetime
import unittest

from eobs import Eobs

tg = Eobs(r'c:\m\ghs\data\e-obs\tg_0.25deg_reg_v12.0.nc')

# Paris coords and corresponding indices:
lat = 48.66
lon = 2.35
lat_idx = 85
lon_idx = 171


class TestEobs(unittest.TestCase):

    def test_lat2index(self):
        lat_idx = tg.lat2index(lat)
        lat2 = tg.variables['latitude'][lat_idx]
        # Check that the latitude referred to by this index represents the
        # grid square that original latitude falls within
        self.assertTrue(lat2 - tg.gridsize <= lat <= lat2 + tg.gridsize)

    def test_lon2index(self):
        lon_idx = tg.lon2index(lon)
        lon2 = tg.variables['longitude'][lon_idx]
        # Check that the latitude referred to by this index represents the
        # grid square that original latitude falls within
        self.assertTrue(lon2 - tg.gridsize <= lon <= lon2 + tg.gridsize)

    def test_getitem_with_date_lat_lon(self):
        t = tg[datetime.datetime(1950, 1, 1), lat, lon]
        # Using lat-lon indices that have been looked up
        self.assertEqual(t, tg.variable[0, lat_idx, lon_idx])