#! /usr/bin/python3

from astropy import units as u
import sys
from astropy.coordinates import Angle, SkyCoord, ICRS, FK5
from astropy.time import Time, TimeDelta

# Get coordinates
star_ra = sys.argv[1].replace(",",".")
#star_ra_in = Angle(star_ra, u.degree)
star_dec = sys.argv[2].replace(",",".")
#star_dec_in = Angle(star_dec, u.degree)
date_of_observation = "2023-01-02T12:34:56"


def getPreciseCoord(ra, dec, date):
    coords = [str(ra) + ' ' + str(dec)]
    coord_object = SkyCoord(coords, frame='icrs', unit=(u.degree, u.degree), obstime=date)
    j2000_coord = coord_object.transform_to(FK5(equinox='J2000.0'))
    j2000_coords = j2000_coord.to_string(style='hmsdms', precision=2)[0]
    j2000_coord_formatted = str(j2000_coords).replace("d",":").replace("h",":").replace("m",":").replace("s","")
    return j2000_coord_formatted

def getUTC(date_time):
    date_of_observation_time_cet = Time(date_time, precision=0)
    time_zone_delta = TimeDelta(-3600, format='sec')
    date_of_observation_time_utc = date_of_observation_time_cet + time_zone_delta
    return str(date_of_observation_time_utc.jyear)

print('Precise coordinates (J2000): ' + str(getPreciseCoord(star_ra, star_dec, date_of_observation)))
print('Date of observation: ' + getUTC(date_of_observation))
