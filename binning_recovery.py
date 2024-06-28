""" 
This file runs a recovery routine on observations which did had 
interrupted pipeline calls pertaining to their binning and 
detection calls.
"""
from helpers import *
from nuanalysis import NuAnalysis


starting_directory = os.getcwd()
dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 22
    counter = 0
    low_phi_file = "./ref_files/nustar_pilow.txt"
    high_phi_file = "./ref_files/nustar_pihi.txt"
    starting_directory = os.getcwd()
    with open(f"./pathfiles/runlist_{run_cycle}.txt", 'r') as run_file:
        run_data = run_file.readlines()
    for idx, datum in enumerate(run_data):
        run_order[idx] = datum.split()

    for idx in run_order:
        object_name = run_order[idx][0]
        seqid = run_order[idx][1]
        path = os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            run_object.write_net_detections()
            os.chdir(starting_directory)
        os.chdir(starting_directory)


dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 29
    counter = 0
    low_phi_file = "./ref_files/nustar_pilow.txt"
    high_phi_file = "./ref_files/nustar_pihi.txt"
    starting_directory = os.getcwd()
    with open(f"./pathfiles/runlist_{run_cycle}.txt", 'r') as run_file:
        run_data = run_file.readlines()
    for idx, datum in enumerate(run_data):
        run_order[idx] = datum.split()

    for idx in run_order:
        object_name = run_order[idx][0]
        seqid = run_order[idx][1]
        path = os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            run_object.write_net_detections()
            os.chdir(starting_directory)
        os.chdir(starting_directory)