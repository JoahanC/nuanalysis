from helpers import *
from nuanalysis import NuAnalysis


# Define script parameters
dtimes = [1000, 500]
ssd_idx = 1
run_cycles = [10]
action = "basicdet1"
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"
check_flag = False


# Sanity_check
valid_cycle_1 = list(range(1, 16))
valid_cycle_2 = list(range(16, 33))


# High detection import 
high_seqids = {}
with open("./pathfiles/high_detection_2.txt") as high_det:
    for line in high_det.readlines():
        line_data = line.replace('\n', '').split()
        seqid = line_data[0]
        dtime = line_data[1]
        if seqid in high_seqids:
            high_seqids[seqid].append(dtime)
        if seqid not in high_seqids:
            high_seqids[seqid] = [dtime]

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
            
            if action == "basicdet1":
                test_flag = os.path.join(evdir, f"{dtime}_det1basic_flag.txt")
                if os.path.isfile(test_flag):
                    #print("Already processed")
                    continue
                
            if action == "mainsource":
                test_flag = os.path.join(evdir, f"{dtime}_det1basic_flag.txt")
                if not os.path.isfile(test_flag):
                    print("Please process first")
                    continue
            
            out_path = os.path.join(path, "products")
            if seqid in high_seqids:
                if dtime in high_seqids[seqid]:
                    #print(f"High detection case: {dtime}! Skipping")
                    continue
                
            print(seqid, dtime)
            
            if seqid not in high_seqids:
                run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                                seqid=seqid, clean=True, bifrost=True, object_name=object_name)
                if run_object.completeness_flag:
                    
                    if action == "poissondet1":
                        run_object.convert_to_det1()
                        
                    if action == "poisson":
                        run_object.recalculate_poisson()
                        
                    if action == "straycut":
                        run_object.straycat_removal()
                        
                    if action == "basicdet1":
                        run_object.convert_basic_to_det1()
                        
                    if action == "mainsource":
                        run_object.remove_main_source()
                    
                os.chdir(starting_directory)
