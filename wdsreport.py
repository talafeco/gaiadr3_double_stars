import csv
import os
import sys
import numpy as np
import datetime
import math
import sys
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import QTable
from astropy.table import Table, vstack
from photutils.detection import DAOStarFinder
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
import warnings
from io import StringIO
from astropy.io import ascii
warnings.filterwarnings("ignore")

""" wds_identifier = np.array([], dtype=str)
discovr = np.array([], dtype=str)
comp = np.array([], dtype=str)
theta = np.array([], dtype=np.float64)
rho = np.array([], dtype=np.float64)
mag_pri = np.array([], dtype=np.float64)
mag_sec = np.array([], dtype=np.float64)
spectra = np.array([], dtype=str)
pm_a_ra = np.array([], dtype=np.float64)
pm_a_dec = np.array([], dtype=np.float64)
pm_b_ra = np.array([], dtype=np.float64)
pm_b_dec = np.array([], dtype=np.float64)
ra = np.array([], dtype=str)
dec = np.array([], dtype=str)
wdsTable = QTable([wds_identifier, discovr, comp, theta, rho, mag_pri, mag_sec, spectra, pm_a_ra, pm_a_dec, pm_b_ra, pm_b_dec, ra, dec], names=('wds_identifier', 'discovr', 'comp', 'theta', 'rho', 'mag_pri', 'mag_sec', 'spectra', 'pm_a_ra', 'pm_a_dec', 'pm_b_ra', 'pm_b_dec', 'ra', 'dec'), meta={'name': 'source table'}) """

converterslist = {'wds_identifier': ('U10'),
              'discovr': ('U10'),
              'comp': ('U10'),
              'theta': np.float64,
              'rho': np.float64,
              'mag_pri': np.float64,
              'mag_sec': np.float64,
              'spectra': ('U10'),
              'pm_a_ra': np.float64,
              'pm_a_dec': np.float64,
              'pm_b_ra': np.float64,
              'pm_b_dec': np.float64,
              'ra': ('U10'),
              'dec': ('U10')}

wdsfilename = sys.argv[1]
#doubleStars = np.empty((0, 13))
doubleStars = Table.read(wdsfilename, converters=converterslist)

print(doubleStars.info)