""" astropy_test.py

Created on Mar 15, 2023
@author: Shin Satoh

Description:


Version
1.0.0 (Mar 15, 2023) Coordinate conversion

"""

# %% Importing libraries
import numpy as np
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, solar_system_ephemeris, EarthLocation
from astropy.coordinates import get_body_barycentric, get_body, get_moon

t = Time('2014-09-22 23:22')
loc = EarthLocation.of_site('greenwich')
with solar_system_ephemeris.set('builtin'):
    jup = get_body('jupiter', t, loc)
print(jup)
print(jup.ra.degree)
print(jup.dec.degree)
print(type(jup))
