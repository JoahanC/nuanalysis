"""
The file runs a routine to take extract transient detections that are outside of the galactic plane.
"""
import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io.fits import getdata, getheader
from tqdm import tqdm

ssd_val = 2
time_val = 5000
detpath = f"ssd{ssd_val}_jan9_{time_val}_nss.tbl"
detections = Table.read(detpath, format='ipac')
galactic_detections = {}
non_galactic_detections = {}

for col in detections.colnames:
    galactic_detections[col] = []
    non_galactic_detections[col] = []

for idx in tqdm(range(len(detections["INDEX"]))):
    flag = True
    ra = detections["RA"][idx]
    dec = detections["DEC"][idx]
    detection_pos = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')
    
    # Calculate the galactic coordinates for this source
    if np.abs(detection_pos.galactic.b.deg) < 15:
        print("Detection discarded. Found in the galactic plane.")
        for col in detections.colnames:
            galactic_detections[col].append(detections[col][idx])
    
    else:
        for col in detections.colnames:
            non_galactic_detections[col].append(detections[col][idx])

detect_table = Table()
for key in galactic_detections:
    detect_table[key] = galactic_detections[key]
finalpath = f"ssd{ssd_val}_jan9_{time_val}_gal.tbl"
detect_table.write(finalpath, format='ipac', overwrite=True)

detect_table2 = Table()
for key in galactic_detections:
    detect_table2[key] = non_galactic_detections[key]
finalpath = f"ssd{ssd_val}_jan9_{time_val}_nongal.tbl"
detect_table2.write(finalpath, format='ipac', overwrite=True)