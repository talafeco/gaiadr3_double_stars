#! /usr/bin/python3

# Usage: gdr3_comp <Image_folder>

import sys
import os
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky, position_angle, angular_separation

first_ra = sys.argv[1]
first_dec = sys.argv[2]
second_ra = sys.argv[3]
second_dec = sys.argv[4]

first_object = SkyCoord(ra=first_ra*u.degree, dec=first_dec*u.degree)
second_object = SkyCoord(ra=second_ra*u.degree, dec=second_dec*u.degree)

object_position_angle = first_object.position_angle(second_object)
object_separation = first_object.separation(second_object)

print('First object coordinates (Ra/Dec): ', first_ra, '° ', first_dec, '°', '\n', 'Second object coordinates (Ra/Dec): ', second_ra, '° ', second_dec, '°', '\n', 'Position angle: ', object_position_angle, '°', '\n', 'Separation: ', object_separation*u.arcsecond, '"')