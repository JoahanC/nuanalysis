import numpy as np
import glob
import os


for i in range(1, 32):
    basepath = f"./../bifrost_data/{i}/"
    seqids = os.listdir(basepath)
    for seqid in seqids:
        seqpath = basepath + f"{seqid}/event_cl/"
        files = os.listdir(seqpath)
        if "5000_flag.txt" in files:
            print(i, seqid)
