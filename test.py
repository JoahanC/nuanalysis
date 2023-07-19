from helpers import *
from nuanalysis import NuAnalysis
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob


seqid = "30501002002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
test.generate_detections()
test.extract_detections()
