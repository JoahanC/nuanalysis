from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_5.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/5/{seqid}/event_cl/"):
        path = f"../bifrost_data/5/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=5)
        run_object.sliding_cell_detection()

run_order = {}
counter = 0
with open("../test/runlist_15.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/15/{seqid}/event_cl/"):
        path = f"../bifrost_data/15/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=15)
        run_object.sliding_cell_detection()

run_order = {}
counter = 0
with open("../test/runlist_25.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"5000_binning_flag.txt" in os.listdir(f"../bifrost_data/25/{seqid}/event_cl/"):
        path = f"../bifrost_data/25/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=25)
        run_object.sliding_cell_detection()