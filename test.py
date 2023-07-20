from helpers import *
from nuanalysis import NuAnalysis
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob


seqid = "30501002002"
path = f"../data/test/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
test.test_gen()
#test.generate_detections()
#test.extract_detections()
#remove_tmp_folders("../data/30501002002/")