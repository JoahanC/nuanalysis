from helpers import *
from nuanalysis import NuAnalysis
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
test.nuproducts(detect_info)