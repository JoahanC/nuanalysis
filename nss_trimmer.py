"""
The file runs a comparison routine to take final detection positions and discard 
ones that intersect within 1 arcminute of NSS sources.
"""
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.io.fits import getdata, getheader
from tqdm import tqdm 


nssheader1 = getheader("nustarnss.fits", 1)
nssvalues1 = getdata("nustarnss.fits", 1)
nssheader2 = getheader("nustarnss2.fits", 1)
nssvalues2 = getdata("nustarnss2.fits", 1)

ssd_val = 1
time_val = 500
detpath = f"ssd{ssd_val}_jan9_{time_val}_mainmask.tbl"
detections = Table.read(detpath, format='ipac')
filtered_detections = {}

for col in detections.colnames:
    filtered_detections[col] = []

for idx in tqdm(range(len(detections["INDEX"]))):
    flag = True
    ra = detections["RA"][idx]
    dec = detections["DEC"][idx]
    detection_pos = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg), frame='fk5')

    for idx2 in range(len(nssvalues1["name"])):
        nss1ra = nssvalues1["ra"][idx2]
        nss1dec = nssvalues1["dec"][idx2]
        nss_pos1 = SkyCoord(f"{nss1ra} {nss1dec}", unit=(u.deg, u.deg), frame='icrs')
        sep = detection_pos.separation(nss_pos1)
        if sep.arcsecond < 60:
            flag = False
    if not flag:
        print("Detection discarded. Found in NSS Primary.")
        continue
    for idx3 in range(len(nssvalues2["name"])):
        nss2ra = nssvalues2["ra"][idx3]
        nss2dec = nssvalues2["dec"][idx3]
        nss_pos2 = SkyCoord(f"{nss2ra} {nss2dec}", unit=(u.deg, u.deg), frame='icrs')
        sep = detection_pos.separation(nss_pos2)
        if sep.arcsecond < 50:
            flag = False
    if not flag:
        print("Detection discarded. Found in NSS Secondary.")
        
    if flag:
        for col in detections.colnames:
            filtered_detections[col].append(detections[col][idx])

detect_table = Table()
for key in filtered_detections:
    detect_table[key] = filtered_detections[key]
finalpath = f"ssd{ssd_val}_jan9_{time_val}_nss.tbl"
detect_table.write(finalpath, format='ipac', overwrite=True)