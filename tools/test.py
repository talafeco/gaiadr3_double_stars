#! /usr/bin/python3

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

gravConst = 0.0043009
massa = 1.1623341720418554
massb = 0.8373715693669024
sep_pc = 0.019690628023291214
half_axis = 1357195.0254232308

max_orbit_velolicy = math.sqrt(gravConst * ((massa + massb) * (2 / (sep_pc)) - 1 / (half_axis * 0.00000485)))

print(max_orbit_velolicy)