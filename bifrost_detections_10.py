from helpers import *
from nuanalysis import NuAnalysis
from tqdm import tqdm


run_order = {}
counter = 0
with open("../test/runlist_10.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/10/{seqid}/event_cl/"):
        if f"5000_flag.txt" in os.listdir(f"../bifrost_data/10/{seqid}/event_cl/"):
            path = f"../bifrost_data/10/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=10)
            run_object.detection_merging()

run_order = {}
counter = 0
with open("../test/runlist_20.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/20/{seqid}/event_cl/"):
        if f"5000_flag.txt" in os.listdir(f"../bifrost_data/20/{seqid}/event_cl/"):
            path = f"../bifrost_data/20/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=20)
            run_object.detection_merging()

run_order = {}
counter = 0
with open("../test/runlist_30.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/30/{seqid}/event_cl/"):
        if f"5000_flag.txt" in os.listdir(f"../bifrost_data/30/{seqid}/event_cl/"):
            path = f"../bifrost_data/30/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=30)
            run_object.detection_merging()

