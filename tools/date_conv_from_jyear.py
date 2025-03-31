#! /usr/bin/python3
# Convert date object to Julian year format
# Usage: date_conv.py '<datetime>'
# More info: https://docs.astropy.org/en/stable/time/

import sys
from astropy.time import Time, TimeDelta

date_time = sys.argv[1]
date_of_observation_time_cet = Time(date_time, precision=0, format='jyear')

#print('CET:', date_time_cet_jyear)
#print('UTC:', date_time_utc_jyear)
print(str(date_of_observation_time_cet.jyear))
print(str(date_of_observation_time_cet.iso))

