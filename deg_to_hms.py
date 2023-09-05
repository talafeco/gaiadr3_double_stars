from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from astropy.coordinates import SkyCoord

# Get coordinates
star_ra = sys.argv[1].replace(",",".")
star_dec = sys.argv[2].replace(",",".")

print(star_ra + ' ' + star_dec)
print(SkyCoord(star_ra, star_dec, unit='hour, dms'))