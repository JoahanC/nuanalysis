from helpers import *
from nuanalysis import NuAnalysis


run_order = {}
counter = 0
with open("../test/runlist_1.txt", 'r') as run_file:
    run_data = run_file.readlines()
for idx, datum in enumerate(run_data):
    run_order[idx] = datum.split()

for idx in run_order:
    path = f"../bifrost_data/1/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    run_object = NuAnalysis(5000, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=False, bifrost=True)

