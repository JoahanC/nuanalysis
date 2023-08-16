from helpers import *
from nuanalysis import NuAnalysis
from tqdm import tqdm


dtime = 1000
run_order = {}
counter = 0
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
starting_directory = os.getcwd()
with open("../test/runlist_1.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/1/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/1/{seqid}/event_cl/"):
            path = f"../bifrost_data/1/{seqid}/"
            
            run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)
            run_object.write_net_detections()
            run_object.recalculate_poisson()
            os.chdir(starting_directory)

run_order = {}
counter = 0
with open("../test/runlist_11.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/11/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/11/{seqid}/event_cl/"):
            path = f"../bifrost_data/11/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)
            run_object.write_net_detections()
            run_object.recalculate_poisson()
            os.chdir(starting_directory)


run_order = {}
counter = 0
with open("../test/runlist_21.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/21/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/21/{seqid}/event_cl/"):
            path = f"../bifrost_data/21/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)
            run_object.write_net_detections()
            run_object.recalculate_poisson()
            os.chdir(starting_directory)

