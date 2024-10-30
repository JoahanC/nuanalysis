import os 
import subprocess

directory = "./test_data/60160619002/event_cl"
files = os.listdir(directory)
for file in files:
    if ".gz" in file:
        subprocess.run(["mv", os.path.join(directory, file), os.path.join(directory, file.replace(".gz", ""))])
directory = "./test_data/60160619002/auxil"
files = os.listdir(directory)
for file in files:
    if ".gz" in file:
        subprocess.run(["mv", os.path.join(directory, file), os.path.join(directory, file.replace(".gz", ""))])