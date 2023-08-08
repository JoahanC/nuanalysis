from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_4.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/4/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/4/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/4/{seqid}/event_cl/"):
                print(seqid)
                path = f"../bifrost_data/4/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=4)
                run_object.event_extraction()

run_order = {}
counter = 0
with open("../test/runlist_14.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/14/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/14/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/14/{seqid}/event_cl/"):
                print(seqid)
                path = f"../bifrost_data/14/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=14)
                run_object.event_extraction()

run_order = {}
counter = 0
with open("../test/runlist_24.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/24/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/24/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/24/{seqid}/event_cl/"):
                print(seqid)
                path = f"../bifrost_data/24/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=24)
                run_object.event_extraction()