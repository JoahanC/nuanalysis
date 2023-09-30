from nuanalysis import *
from astropy.io.fits import getdata, getheader
import matplotlib.pyplot
import os

starting_directory = os.getcwd()
seqids = ["80102048008"]#["80902312004"]
dtimes = [5000]#, 1000, 500]
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
for dtime in dtimes:
    object_name = "test"
    for seqid in seqids:
        path = f"../bifrost_data/{seqid}/"
        #path = f"/Volumes/data_ssd_1/bifrost_data//{seqid}/"
        evdir = f"{path}event_cl/"
        out_path = f"{path}products/"
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        #run_object.display_detections()
        run_object.straycat_removal()
        #run_object.collect_det1coords()
        os.chdir(starting_directory)


"""data = getdata(f"../bifrost_data/{seqid}/event_cl/nu{seqid}_detpos.fits")
header = getheader(f"../bifrost_data/{seqid}/event_cl/nu{seqid}_detpos.fits")
times = []
detx = []
dety = []
for datum in data:
    times.append(float(datum[0]))
    detx.append(float(datum[1]))
    dety.append(float(datum[2]))

plt.scatter(detx, dety)
plt.xlim(0, 300)
plt.ylim(0, 300)
plt.show()"""