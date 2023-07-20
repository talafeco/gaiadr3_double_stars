#Test
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
from astropy.coordinates import SkyCoord

# Get coordinates
star_ra = sys.argv[1].replace(",",".")[0:1] + 'h' + sys.argv[1].replace(",",".")[2:3] + 'm' + sys.argv[1].replace(",",".")[4:9] + 's'
star_dec = sys.argv[2].replace(",",".")

print(star_ra + ' ' + star_dec)
print(SkyCoord(star_ra, star_dec, unit='hour, degree'))