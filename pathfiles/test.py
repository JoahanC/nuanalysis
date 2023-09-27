import glob
import subprocess
import os
import re
import itertools
from more_itertools import batched
from astropy.io.fits import getheader


file_list = subprocess.run(["ls", "-lct"], cwd="../../../nustar/fltops", encoding="utf-8", capture_output=True)
with open("files.txt", 'w') as file:
    file.write(file_list.stdout)
with open("files.txt", "r") as file:
    file.readline()
    objects = []
    for line in file:
        words = file.readline().split()
        if words[-1]=="ECDFS":
            break
        objects.append(words[-1])

data_paths = {}
for object in objects[15:-15]:
    data_paths[object] = []
    folders = os.listdir(f"../../../nustar/fltops/{object}/")
    trim_folders = []
    for folder in folders:
        if len(folder) == 11:
            data_paths[object].append(folder)

path_stub = "../../../nustar/fltops/"
data_sets = {}
size = 0
exposure_time = 0
for object in data_paths:
    for folder in data_paths[object]:
        files = os.listdir(path_stub + f"{object}/{folder}")
        if "event_cl" in files:
            event_files = os.listdir(path_stub + f"{object}/{folder}/event_cl")
            if f"nu{folder}A01_cl.evt" in event_files and f"nu{folder}B01_cl.evt in event_files":
                evtfile = path_stub + f"{object}/{folder}/event_cl/nu{folder}A01_cl.evt"
                hdr = getheader(evtfile)
                exposure_time += float(hdr['EXPOSURE'])
                if object not in data_sets:
                    data_sets[object] = [folder]
                else:
                    data_sets[object].append(folder)

total_runs = 0
for object in data_sets:
    for folder in data_sets[object]:
        total_runs += 1        
        data_size = subprocess.run(["du", "-sh", "."], cwd=path_stub+f"{object}/{folder}", capture_output=True, encoding="utf-8")
        size_string = data_size.stdout.split()[0]
        if 'M' in size_string:
            size += float(size_string.replace("M", ''))
            continue
        if 'G' in size_string:
            size += float(size_string.replace("G", '')) * 1000
            continue
        else:
            print("edge case", size_string)


running_list = {}
counter = 0
net_objects = []
for object in data_sets:
    for folder in data_sets[object]:
        net_objects.append(object)
        running_list[counter] = [object, folder]
        counter += 1


partition_list = {}
index_list = list(running_list.keys())
partitions = list(batched(index_list, 75))
print(len(partitions))

for i in range(len(partitions)):
    for idx in partitions[i]:
        object = running_list[idx][0]
    #with open(f"runlist_{i}.txt", "w") as file:
    #    for idx in partitions[i]:
    #        object = running_list[idx][0]
    #        net_objects.append(object)
    #        folder = running_list[idx][1]
    #        file.write(f"{object} {folder} ")
    #        file.write(path_stub + f"{object}/{folder}/event_cl/nu{folder}A01_cl.evt ")
    #        file.write(path_stub + f"{object}/{folder}/event_cl/nu{folder}B01_cl.evt\n")

print(len(set(net_objects)))
print(f"Total exposure time: {exposure_time}")
print(f"Total runs: {total_runs}")
print(f"Total data: {size/1000}GB")

