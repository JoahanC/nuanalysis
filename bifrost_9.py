from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_9.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/9/{seqid}"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/9/{seqid}"):
            path = f"../bifrost_data/9/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=9)

run_order = {}
counter = 0
with open("../test/runlist_19.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/19/{seqid}"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/19/{seqid}"):
            path = f"../bifrost_data/19/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=19)

run_order = {}
counter = 0
with open("../test/runlist_29.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/29/{seqid}"):
        if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/29/{seqid}"):
            path = f"../bifrost_data/29/{seqid}/"
            evdir = f"{path}event_cl/"
            out_path = f"{path}products/"
            run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=29)