from helpers import *
from nuanalysis import NuAnalysis
from tqdm import tqdm

dtime = 1000
run_order = {}
counter = 0
with open("../test/runlist_24.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/24/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/24/{seqid}/event_cl/"):
            path = f"../bifrost_data/24/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(dtime, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=22)
            #run_object.detection_merging()
            run_object.write_net_detections()

run_order = {}
counter = 0
with open("../test/runlist_23.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/23/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/23/{seqid}/event_cl/"):
            path = f"../bifrost_data/23/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(dtime, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=22)
            #run_object.detection_merging()
            run_object.write_net_detections()

run_order = {}
counter = 0
with open("../test/runlist_31.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in tqdm(run_order):
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"{dtime}_binning_flag.txt" in os.listdir(f"../bifrost_data/31/{seqid}/event_cl/"):
        if f"{dtime}_flag.txt" in os.listdir(f"../bifrost_data/31/{seqid}/event_cl/"):
            path = f"../bifrost_data/31/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(dtime, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=31)
            #run_object.detection_merging()
            run_object.write_net_detections()

