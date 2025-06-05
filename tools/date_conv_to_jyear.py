#! /usr/bin/python3
# Convert date object to Julian year format
# Usage: python3 date_conv_to_jyear.py 2024-12-01T18:01:31
# More info: https://docs.astropy.org/en/stable/time/

import sys
from astropy.time import Time, TimeDelta

date_time = sys.argv[1]
date_of_observation_time_cet = Time(date_time) # precision=0
date_time_cet_jyear = str(date_of_observation_time_cet.jyear)

time_zone_delta = TimeDelta(3600, format='sec')
date_time_utc_jyear = str((date_of_observation_time_cet - time_zone_delta).jyear)
date_time_utc_jyear2 = str((date_of_observation_time_cet - time_zone_delta))

#print('CET:', date_time_cet_jyear)
#print('UTC:', date_time_utc_jyear)
#print(date_of_observation_time_cet)
print(date_time_utc_jyear)
#print(date_time_utc_jyear2)