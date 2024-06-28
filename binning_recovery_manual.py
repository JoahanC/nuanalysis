""" 
This script runs a verification check on observations which may have 
had interrupted pipeline calls or otherwise some issue with their 
processing.
"""
from helpers import *
from nuanalysis import NuAnalysis


starting_directory = os.getcwd()
run_order = []
seqids = []
dtimes = []

# Create a text file which tracks the all sequence ids to be inspected
with open("recovery_list_3.txt", 'r') as file:
    for line in file.readlines():
        datum = line.split()
        run_order.append(datum[0])
        seqids.append(datum[1])
        dtimes.append(datum[2])


for idx in range(len(seqids)):
    run_cycle = int(run_order[idx])
    seqid = int(seqids[idx])
    dtime = int(dtimes[idx])
    object_name = 'test'
    low_phi_file = "./ref_files/nustar_pilow.txt"
    high_phi_file = "./ref_files/nustar_pihi.txt"
    starting_directory = os.getcwd()
    
    path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
    evdir = os.path.join(path, "event_cl")
    out_path = os.path.join(path, "products")
    run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)
    
    # Check if observation has a valid completeness flag, if not then 
    # rerun pipeline.
    if run_object.completeness_flag:
        run_object.event_extraction()
        os.chdir(starting_directory)
    os.chdir(starting_directory)
