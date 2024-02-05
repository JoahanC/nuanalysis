"""
This script runs lightcurve diagnostics on a specific observation.
"""
from nuanalysis import *
from astropy.io.fits import getdata, getheader
import os

seqids = ["80902404002", "80902404004", "80902404006", "80902404008"]
starting_directory = os.getcwd()

dtimes = [5000]#, 1000, 500]
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
for dtime in dtimes:
    object_name = "test"
    for seqid in seqids:
        path = f"../bifrost_data/{seqid}/"
        #path = f"/Volumes/data_ssd_2/16/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        print(run_object.source_position)
        print(min(run_object._event_times), max(run_object._event_times), max(run_object._event_times) - min(run_object._event_times))
        #run_object.detection_lightcurve(float(run_object._source_pix_coordinates[0][0]), float(run_object._source_pix_coordinates[0][1]), 34, 1934, min(run_object._event_times), max(run_object._event_times))
        os.chdir(starting_directory)
