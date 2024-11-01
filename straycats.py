""" 
This script generates a text file which summarizes all of the 
StrayCats regions files and sorts them by module and sequence 
id for running `nuanalysis`.
"""
import os
import glob

# Aggregate all straycats region files and organize them into a list
path_stub = "/Volumes/data_ssd_1/bifrost_data/straycat/"
strayfiles = glob.glob(f"{path_stub}*")
nonfiles = ["import_regions.py", "aggregate_regions.py", "stray_seqids.txt"]

# Loop through all straycats region files and write them to a nice 
# summary file.
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
