from helpers import *
from nuanalysis import NuAnalysis


starting_directory = os.getcwd()
dtimes = [5000]
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
                vals = run_object.verify_event_extraction()
                os.chdir(starting_directory)
                with open("binning_check_3.txt", 'a') as file:
                    file.write(f"{run_cycle} {seqid} {dtime} {vals[0]} {vals[1]} {vals[2]} {vals[3]}\n")
            os.chdir(starting_directory)