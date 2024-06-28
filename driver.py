""" 
This is a basic script for calling nuanalysis on a specific observation.
"""
from nuanalysis import *
import os

starting_directory = os.getcwd()

# Set the sequence ids and characteristic times for processing.
seqids = ["80902404002"]
dtimes = [5000]#, 1000, 500]

low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"

for dtime in dtimes:
    object_name = "test"
    for seqid in seqids:
        path = f"../bifrost_data/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        os.chdir(starting_directory)
