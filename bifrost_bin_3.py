from helpers import *
from nuanalysis import NuAnalysis



run_order = {}
counter = 0
with open("../test/runlist_3.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/3/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/3/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/3/{seqid}/event_cl/"):
                path = f"../bifrost_data/3/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=1)
                run_object.event_extraction()

run_order = {}
counter = 0
with open("../test/runlist_13.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/13/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/13/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/13/{seqid}/event_cl/"):
                path = f"../bifrost_data/13/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=11)
                run_object.event_extraction()

run_order = {}
counter = 0
with open("../test/runlist_23.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/23/{seqid}/event_cl/"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/23/{seqid}/event_cl/"):
            if "1000_binning_flag.txt" not in os.listdir(f"../bifrost_data/23/{seqid}/event_cl/"):
                path = f"../bifrost_data/23/{seqid}/"
                evdir = f"{path}event_cl/"
                out_path = f"{path}products/"
                run_object = NuAnalysis(1000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=21)
                run_object.event_extraction()