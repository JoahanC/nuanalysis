import os 
import glob
import numpy as np


ref_files = glob.glob("./pathfiles/runlist_*.txt", recursive=True)
targets = []
obsids = []
for file in ref_files:
    with open(file, 'r') as fl:
        for line in fl.readlines():
            datum = line.split()
            targets.append(datum[0])
            obsids.append(datum[1])
print(len(set(targets)))
print(len(set(obsids)))
