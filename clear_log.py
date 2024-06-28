"""
Housekeeping script created to remove unnecesssary log files which 
bloat post-observation folders.
"""
import os
import glob

ssd_val = 1
file_set = 2
if ssd_val == 1:
    log_files = glob.glob(f"/Volumes/data_ssd_{ssd_val}/bifrost_data/{file_set}/**/logs/*.log")
if ssd_val == 2:
    log_files = glob.glob(f"/Volumes/data_ssd_{ssd_val}/{file_set}/**/logs/*.log")

print(log_files)
for file in log_files:
    os.remove(file)
print("Removed all log files!")