from nuanalysis import *
from astropy.io.fits import getdata, getheader
import matplotlib.pyplot
import os

starting_directory = os.getcwd()
seqids = ["80602315006"]
dtimes = [500]#, 1000, 500]
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
for dtime in dtimes:
    object_name = "test"
    for seqid in seqids:
        #path = f"../bifrost_data/{seqid}/"
        path = f"/Volumes/data_ssd_1/bifrost_data/8/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        run_object.display_detections()
        run_object.display_detections_det()
        #run_object.straycat_removal()
        #run_object.collect_det1coords()
        os.chdir(starting_directory)