from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_7.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/7/{seqid}/event_cl/"):
        if f"5000_flag.txt" not in os.listdir(f"../bifrost_data/7/{seqid}/event_cl/"):
            path = f"../bifrost_data/7/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=7)
            run_object.sliding_cell_detection()

run_order = {}
counter = 0
with open("../test/runlist_17.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/17/{seqid}/event_cl/"):
        if f"5000_flag.txt" not in os.listdir(f"../bifrost_data/17/{seqid}/event_cl/"):
            path = f"../bifrost_data/17/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=17)
            run_object.sliding_cell_detection()

run_order = {}
counter = 0
with open("../test/runlist_27.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/27/{seqid}/event_cl/"):
        if f"5000_flag.txt" not in os.listdir(f"../bifrost_data/27/{seqid}/event_cl/"):
            path = f"../bifrost_data/27/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=27)
            run_object.sliding_cell_detection()