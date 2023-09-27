from helpers import *
from nuanalysis import NuAnalysis


starting_directory = os.getcwd()
"""dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 1
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 < 20:
                    run_object.recalculate_poisson()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)

dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 7
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 < 20:
                    run_object.recalculate_poisson()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)


dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 14
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 < 20:
                    run_object.recalculate_poisson()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)

"""






starting_directory = os.getcwd()
dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 1
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 >= 20:
                    run_object.convert_detection_coords()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)

dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 7
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 >= 20:
                    run_object.convert_detection_coords()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)


dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    run_cycle = 14
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
        #if f"{dtime}_binning_flag.txt" not in os.listdir(os.path.relpath(f"/Volumes/data_ssd_2/{run_cycle}/{seqid}/event_cl/")):
        path = os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/")
        evdir = os.path.join(path, "event_cl")
        out_path = os.path.join(path, "products")
        run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                        seqid=seqid, clean=True, bifrost=True, object_name=object_name)
        if run_object.completeness_flag:
            if f"{dtime}_det1convert_flag.txt" in os.listdir(os.path.relpath(f"/Volumes/data_ssd_1/bifrost_data/{run_cycle}/{seqid}/event_cl/")):
                detections, flag, filey = run_object.read_final_detections(detectiontype="basic")
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
                print(val2)
                if val2 >= 20:
                    run_object.convert_detection_coords()
                os.chdir(starting_directory)
            else:
                print("NULL or HEAVY DETECTIONS")
                os.chdir(starting_directory)
        os.chdir(starting_directory)


