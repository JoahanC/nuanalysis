from helpers import *
from nuanalysis import NuAnalysis
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob

fits_path = "./../data/30501002002/detections/559-1934_10000-3/"
det_files = glob.glob(fits_path)
det_file = det_files[0]

detect_info = read_detection_dir(fits_path)
for info in detect_info:
    for key in detect_info[info]:
        print(detect_info[info][key])
   
seqid = "30501002002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
trimmed_detect_info = {}
for time in detect_info:
    for i in detect_info[time]["INDEX"]:
        ra = detect_info[time]["RA"][int(i) - 1]
        dec = detect_info[time]["DEC"][int(i) - 1]
        c = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
        print(c.ra.deg, c.dec.deg)
        print(test._source_position.separation(c).arcsec)
        if test._source_position.separation(c).arcsec > 50:
            if time not in trimmed_detect_info.keys():
                trimmed_detect_info[time] = {}
                for key in detect_info[time]:
                    trimmed_detect_info[time][key] = []
            for key in detect_info[time]:
                trimmed_detect_info[time][key].append(detect_info[time][key][int(i) - 1])

first_time = list(detect_info.keys())[0]
print(len(detect_info[first_time]["INDEX"]))
print(len(trimmed_detect_info[first_time]["INDEX"]))

for i in detect_info[first_time]["INDEX"]:
    ra = detect_info[first_time]["RA"][int(i) - 1]
    dec = detect_info[first_time]["DEC"][int(i) - 1]
    c = SkyCoord(f"{ra} {dec}", unit=(u.hourangle, u.deg))
    print(c.ra.deg, c.dec.deg)
    print(test._source_position.separation(c).arcsec)