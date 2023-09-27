from nuanalysis import *
from astropy.io.fits import getdata, getheader
import matplotlib.pyplot
import os

starting_directory = os.getcwd()
seqids = glob.glob("/Volumes/data_ssd_2/23/*")
obsids = []
for seqid in seqids:
    obsids.append(seqid.replace("/Volumes/data_ssd_2/23/", ''))

low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
object_name = "test"
bright = 0
for seqid in obsids:
    path = f"/Volumes/data_ssd_2/23/{seqid}/"
    evdir = f"{path}event_cl/"
    out_path = f"{path}products/"
    run_object = NuAnalysis(5000, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
    if run_object.determine_bright():
        bright += 1
    os.chdir(starting_directory)
print(bright)