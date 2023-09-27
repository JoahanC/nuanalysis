from helpers import *
from nuanalysis import NuAnalysis
from astropy.io import fits



starting_directory = os.getcwd()
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
starting_directory = os.getcwd()

evt_files = glob.glob("/Volumes/data_ssd_1/bifrost_data/4/60160684001/detections/34-159_500-3/*.evt")
print(evt_files)
for file in evt_files:
    print(file)
    with fits.open(file) as testing:
        print(file, "Passed")

object_name = "test"
run_cycle = 4
seqid = 60160684001
dtime = 500
path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
evdir = os.path.join(path, "event_cl")
out_path = os.path.join(path, "products")
run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                seqid=seqid, clean=True, bifrost=True, object_name=object_name)
if run_object.completeness_flag:
    vals = run_object.determine_bright()
    detections, flag = run_object.read_final_detections()
    val2 = 0
    NoneType = type(None)
    flag = False
    if type(detections) == NoneType:
        val2 = 0
        flag = True
    if not flag:
        if len(detections["INDEX"]) == 0:
            val2 = 0
    if not flag:
        val2 = len(detections["INDEX"])
    print(vals)
    print(val2)


"""with open(f"./recovery_list_2.txt", 'r') as run_file:
    run_data = run_file.readlines()
for rdatum in run_data:
    object_name = "test"
    datum = rdatum.split()
    run_cycle = int(datum[0])
    seqid = int(datum[1])
    dtime = int(datum[2])
    path = os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/")
    evdir = os.path.join(path, "event_cl")
    out_path = os.path.join(path, "products")
    run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)
    if run_object.completeness_flag:
        vals = run_object.verify_event_extraction()
        os.chdir(starting_directory)
        with open("binning_check_2.txt", 'a') as file:
            file.write(f"{run_cycle} {seqid} {dtime} {vals[0]} {vals[1]} {vals[2]} {vals[3]}\n")
    os.chdir(starting_directory)
30001011007"""