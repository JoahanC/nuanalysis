import os
import glob


log_files = glob.glob("/Volumes/data_ssd_1/bifrost_data/12/**/logs/*.log")

print(log_files)
for file in log_files:
    os.remove(file)