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
    if seqid not in os.listdir(f"../bifrost_data/3/"):
        path = f"../bifrost_data/3/{seqid}/"
        nupath = f"./1/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=3)

run_order = {}
counter = 0
with open("../test/runlist_13.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if seqid not in os.listdir(f"../bifrost_data/13/"):
        path = f"../bifrost_data/13/{seqid}/"
        nupath = f"./1/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=13)

run_order = {}
counter = 0
with open("../test/runlist_23.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    if seqid not in os.listdir(f"../bifrost_data/23/"):
        path = f"../bifrost_data/23/{seqid}/"
        nupath = f"./1/{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=23)