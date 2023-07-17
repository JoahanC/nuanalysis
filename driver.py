from nuanalysis import *

seqid = "30501002002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
#test.generate_detections()
detect_info = test.read_detections()
