from astropy import units as u
import sys
from astropy.coordinates import Angle

# Get coordinates
star_ra = sys.argv[1].replace(",",".")
star_ra = Angle(star_ra, u.degree)
star_dec = sys.argv[2].replace(",",".")
star_dec = Angle(star_dec, u.degree)

print(star_ra.to_string(unit=u.hour) + ' ' + star_dec.to_string(unit=u.degree))