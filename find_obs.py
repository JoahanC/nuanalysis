import glob 
import os


testid = "80102048008"
pathfiles = glob.glob("./pathfiles/runlist_*")
for pathfile in pathfiles:
    with open(pathfile, 'r') as file:
        data = file.readlines()
    for line in data:
        seqid = line.split()[1]
        if testid == seqid:
            print(pathfile)
            print(pathfile.replace("./pathfiles/runlist_", '').replace(".txt", ''))
            break