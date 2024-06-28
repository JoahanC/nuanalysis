"""
This script crudely estimates the intensity of main sources by 
estimating the number of main source counts through the native 
`nuanalysis` estimate routine.
"""
from helpers import *
from nuanalysis import NuAnalysis


starting_directory = os.getcwd()
dtimes = [5000, 1000, 500]
run_ids = list(range(1, 16))
run_ids.remove(11)
print(run_ids)
for dtime in dtimes:
    for run_cycle in run_ids:
        run_order = {}
        counter = 0
        low_phi_file = "./ref_files/nustar_pilow.txt"
        high_phi_file = "./ref_files/nustar_pihi.txt"
        starting_directory = os.getcwd()
        with open(f"./pathfiles/runlist_{run_cycle}.txt", 'r') as run_file:
            run_data = run_file.readlines()
        for idx, datum in enumerate(run_data):
            run_order[idx] = datum.split()

        for idx in run_order:
            object_name = run_order[idx][0]
            seqid = run_order[idx][1]
            path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
            evdir = os.path.join(path, "event_cl")
            out_path = os.path.join(path, "products")
            run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                            seqid=seqid, clean=True, bifrost=True, object_name=object_name)
            if run_object.completeness_flag:
                val = run_object.determine_bright()
                detections, flag = run_object.read_final_detections(detectiontype="basic")
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
                os.chdir(starting_directory)
                with open("bright.txt", 'a') as file:
                    file.write(f"{run_cycle} {seqid} {dtime} {val} {val2}\n")
            os.chdir(starting_directory)
