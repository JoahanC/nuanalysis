import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from astropy.io.fits import getdata, getheader
from tqdm import tqdm


data_path = "./../bifrost_data/"
valid_data_files = {}
event_files = glob.glob(data_path + "**/**/event_cl/**A01_cl.evt")
event_doublets = []
for file in event_files:
    event_doublets.append((file, file.replace("A01", "B01")))

bright_sources = []
for doublet in tqdm(event_doublets):
    data = getdata(doublet[0])
    header = getheader(doublet[0])    
    test_set = [i[0] for i in data]
    test_indices = np.arange(np.min(test_set), np.max(test_set), 1)
    histy = np.histogram(test_set, bins=test_indices)   
    if np.max(histy[0]) > 100:
        bright_sources.append(doublet[0].replace("nu", '').replace("A01_cl.evt", ''))

with open("bright_sources.txt", 'w') as file:
    for source in bright_sources:
        file.write(f"{source}\n")

