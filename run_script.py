from helpers import *
from nuanalysis import NuAnalysis


# Define script parameters
dtimes = [1000, 500]
ssd_idx = 1
run_cycles = [1]
action = "poissondet1"
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
check_flag = False


# Sanity_check
valid_cycle_1 = list(range(1, 16))
valid_cycle_2 = list(range(16, 33))

if ssd_idx == 1:
    for cycle in run_cycles:
        if cycle not in valid_cycle_1:
            raise ValueError(f"Invalid run cycle found: {cycle}")
if ssd_idx == 2:
    for cycle in run_cycles:
        if cycle not in valid_cycle_2:
            raise ValueError(f"Invalid run cycle found: {cycle}")


# Run script
starting_directory = os.getcwd()
for dtime in dtimes:
    for run_cycle in run_cycles:
        run_order = {}
        counter = 0

        with open(f"./pathfiles/runlist_{run_cycle}.txt", 'r') as run_file:
            run_data = run_file.readlines()
        for idx, datum in enumerate(run_data):
            run_order[idx] = datum.split()

        for idx in run_order:
            object_name = run_order[idx][0]
            seqid = run_order[idx][1]
            if ssd_idx == 1:
                path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
            if ssd_idx == 2:
                path = os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/")
            evdir = os.path.join(path, "event_cl")
            out_path = os.path.join(path, "products")
            run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
            if run_object.completeness_flag:
                
                if action == "poissondet1":
                    run_object.convert_to_det1()
                
            os.chdir(starting_directory)
