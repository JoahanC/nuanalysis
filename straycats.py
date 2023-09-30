import os
import glob
import numpy as np
import matplotlib.pyplot as plt


path_stub = "/Volumes/data_ssd_1/bifrost_data/straycat/"
strayfiles = glob.glob(f"{path_stub}*")
nonfiles = ["import_regions.py", "aggregate_regions.py", "stray_seqids.txt"]

for filename in strayfiles:
    file_stub = filename.replace(path_stub, '')
    if file_stub in nonfiles:
        continue
    seqid = file_stub[:11]
    module = file_stub[11]
    catalog = file_stub.replace('_', '').replace(seqid, '').replace(module, '').replace(".reg", '')
    
    if "straydir.txt" not in os.listdir("."):
        with open("straydir.txt", 'w') as file:
            file.write(f"{file_stub} {seqid} {module} {catalog}\n")
    with open("straydir.txt", 'a') as file:
        file.write(f"{file_stub} {seqid} {module} {catalog}\n")
