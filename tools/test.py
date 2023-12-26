import os
import numpy as np
import math
import sys
from astropy.coordinates import SkyCoord, Angle, FK5
from astropy.coordinates import match_coordinates_sky
from astropy.coordinates import search_around_sky
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack, hstack
from astropy.table import Column, MaskedColumn
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
import warnings
from matplotlib import pyplot as plt
from astropy.wcs import utils
from astropy.time import Time, TimeDelta


date_time = sys.argv[1]
date_of_observation_time_cet = Time(date_time, precision=0)
utc_date_time = str(date_of_observation_time_cet.jyear)
print(utc_date_time)