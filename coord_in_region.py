from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion, PixCoord, CirclePixelRegion

#sky_center = SkyCoord(42, 43, unit='deg')
#sky_radius = Angle(25, 'deg')
#sky_region = CircleSkyRegion(sky_center, sky_radius)
#print(sky_region)

#skypoint = SkyCoord(40, 50, unit='deg')
#skypoint in sky_region

#Region: CircleSkyRegion
#center: <SkyCoord (ICRS): (ra, dec) in deg
    #(42., 43.)>
#radius: 25.0 deg

#c1 = SkyCoord('5h23m34.5s', '-69d45m22s', frame='icrs')
#c2 = SkyCoord('0h52m44.8s', '-72d49m43s', frame='fk5')
#sep = c1.separation(c2)
#print(sep.arcsecond)

