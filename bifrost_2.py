from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_2.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    path = f"../bifrost_data/2/{seqid}/"
    nupath = f"./1/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=2)

run_order = {}
counter = 0
with open("../test/runlist_12.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    path = f"../bifrost_data/12/{seqid}/"
    nupath = f"./1/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=12)

run_order = {}
counter = 0
with open("../test/runlist_22.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    object_name = run_order[idx][0]
    seqid = run_order[idx][1]
    path = f"../bifrost_data/22/{seqid}/"
    nupath = f"./1/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True, object_name=object_name, nupath=nupath, sessionid=22)