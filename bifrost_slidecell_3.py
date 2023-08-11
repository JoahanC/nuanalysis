from helpers import *
from nuanalysis import NuAnalysis
from tqdm import tqdm


run_order = {}
counter = 0
with open("../test/runlist_3.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"1000_binning_flag.txt" in os.listdir(f"../bifrost_data/3/{seqid}/event_cl/"):
        if f"1000_flag.txt" not in os.listdir(f"../bifrost_data/3/{seqid}/event_cl/"):
            path = f"../bifrost_data/3/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=2)
            run_object.sliding_cell_detection()

run_order = {}
counter = 0
with open("../test/runlist_13.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"1000_binning_flag.txt" in os.listdir(f"../bifrost_data/13/{seqid}/event_cl/"):
        if f"1000_flag.txt" not in os.listdir(f"../bifrost_data/13/{seqid}/event_cl/"):
            path = f"../bifrost_data/13/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=12)
            run_object.sliding_cell_detection()

