import os
import subprocess
import glob
import random
import numpy as np
import matplotlib.pyplot as plt
from nustar import *
from helpers import *
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
from astropy.table import QTable, Table, Column
from tqdm import tqdm
import radial_profile
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion


basepath = "./../bifrost_data/**/**/detections/5000.tbl"
table_files = glob.glob(basepath)
print(len(table_files))