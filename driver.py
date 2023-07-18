from nuanalysis import *

seqid = "30501002002"
path = f"../data/{seqid}/"
evdir = f"{path}event_cl/"
out_path = f"{path}products/"
test = NuAnalysis(10000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True)
#print(test._evdir)
#print(test._refoutpath)
#print(test._refpath)
#test.generate_detections()
#for i in [0, 1]:
#    for interval in test._time_bins[i]:
#        print(interval, test._time_bins[i][interval][0], test._time_bins[i][interval][1])
detections = test.detection_dir_processing(test.phi_bounds[1])
#test.nuproducts(detections)
test.ds9_detections()
#print(detections)