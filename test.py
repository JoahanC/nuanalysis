from helpers import *
from nuanalysis import NuAnalysis
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob

seqids = []
data_dirs = (glob.glob("./../data/*"))
for dir in data_dirs:
    seqids.append(dir.replace("./../data/", ''))

for seqid in seqids:
    path = f"../data/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    test = NuAnalysis(10000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False)
    test.test_gen()
    test.extract_detections()
    remove_tmp_folders(path)
