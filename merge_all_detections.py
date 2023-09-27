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
from astropy.table import QTable, Table, Column, vstack
from tqdm import tqdm
import radial_profile
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion

cols = ['INDEX', 'COUNTS', 'XPIX', 'YPIX', 'VIGCOR', 'ERR', 'HBOX', 'PROB', 'SNR', 'TSTART', 'TSTOP', 'COUNTSERR', 'SEP']
scols = ['RA', 'DEC', 'SEQID', 'BOUND']
basepath = "/Volumes/data_ssd_1/bifrost_data/**/**/detections/500_poisson.tbl"
table_files = glob.glob(basepath)
print(table_files)
total_table = Table.read(table_files[0], format='ipac')
for file in tqdm(table_files[1:]):
    data_table = Table.read(file, format='ipac')
    for col in cols:
        data_table[col] = data_table[col].astype(float)
        total_table[col] = total_table[col].astype(float)
    for scol in scols:
        data_table[scol] = data_table[scol].astype(str)
        total_table[scol] = total_table[scol].astype(str)
    if len(data_table['RA']) == 0:
        continue
    total_table = vstack([total_table, data_table])
print(len(total_table['RA']))
total_table.write("poisson_1000_september.tbl", format='ipac', overwrite=True)
