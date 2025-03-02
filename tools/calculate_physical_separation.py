#! /usr/bin/python3

# Usage: gdr3_comp <Image_folder>

import sys
import os
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky, position_angle, angular_separation

auToParsec = 0.0000048481368111358

def calculate_3d_distance(parallax1_mas, parallax2_mas, angular_separation_arcsec):
    """
    Calculate the distance between two stars given their parallax values and angular separation.
    
    Parameters:
    parallax1_mas (float): Parallax of the first star in milliarcseconds (mas).
    parallax2_mas (float): Parallax of the second star in milliarcseconds (mas).
    angular_separation_arcsec (float): Angular separation between the two stars in arcseconds.
    
    Returns:
    float: The distance between the two stars in parsecs.
    """
    
    # Convert parallaxes from milliarcseconds to parsecs
    d1 = 1000.0 / parallax1_mas
    d2 = 1000.0 / parallax2_mas
    
    # Convert angular separation from arcseconds to radians
    theta_radians = np.deg2rad(angular_separation_arcsec / 3600.0)
    
    # Apply the law of cosines to find the distance between the stars
    distance_parsecs = np.sqrt(d1**2 + d2**2 - 2 * d1 * d2 * np.cos(theta_radians))
    
    # Convert distance from parsecs to light-years
    distance_lightyears = distance_parsecs * 3.262

    distance_au = distance_parsecs / auToParsec
    
    return distance_parsecs, distance_lightyears, distance_au

first_ra = float(sys.argv[1])
first_dec = float(sys.argv[2])
second_ra = float(sys.argv[3])
second_dec = float(sys.argv[4])
first_par = float(sys.argv[5])
second_par = float(sys.argv[6])

first_object = SkyCoord(ra=first_ra*u.degree, dec=first_dec*u.degree)
second_object = SkyCoord(ra=second_ra*u.degree, dec=second_dec*u.degree)

object_position_angle = first_object.position_angle(second_object).degree
object_separation = first_object.separation(second_object).arcsecond

print('### Calculate Separation between two celestial coordinates ###', '\nFirst object coordinates (Ra/Dec):', first_object.ra.degree, '°,', first_dec, '°', '\nSecond object coordinates (Ra/Dec): ', second_object.ra.degree, '°, ', second_dec, '°', '\nPosition angle:', object_position_angle, '°', '\nSeparation:', object_separation, '"')

physical_distance = calculate_3d_distance(first_par, second_par, object_separation)

print('### Calculate Physical Separation between two celestial coordinates ###', '\nSeparation in parsecs: ', physical_distance[0], '\nSeparation in lightyears: ', physical_distance[1], '\nSeparation in astronomical units: ', physical_distance[2])