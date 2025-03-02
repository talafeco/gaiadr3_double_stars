#! /usr/bin/python3

# Usage: gdr3_comp <Image_folder>

import sys
import os
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky, position_angle, angular_separation

first_ra = float(sys.argv[1])
first_dec = float(sys.argv[2])
second_ra = float(sys.argv[3])
second_dec = float(sys.argv[4])

first_object = SkyCoord(ra=first_ra*u.hourangle, dec=first_dec*u.degree)
second_object = SkyCoord(ra=second_ra*u.hourangle, dec=second_dec*u.degree)

object_position_angle = first_object.position_angle(second_object).degree
object_separation = first_object.separation(second_object).arcsecond

print('### Calculate Separation between two celestial coordinates ###\n', 'First object coordinates (Ra/Dec):', first_object.ra.degree, '°,', first_dec, '°', '\n', 'Second object coordinates (Ra/Dec): ', second_object.ra.degree, '°, ', second_dec, '°', '\n', 'Position angle:', object_position_angle, '°', '\n', 'Separation:', object_separation, '"')