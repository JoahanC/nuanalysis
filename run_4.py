from helpers import *
from nuanalysis import NuAnalysis

seqid = "70901004002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(100, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
test.test_gen()
test.extract_detections()
remove_tmp_folders(path)
