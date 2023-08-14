from helpers import *
from nuanalysis import NuAnalysis
from tqdm import tqdm


starting_directory = os.getcwd()
dtimes = [5000, 1000, 500]
for dtime in dtimes:
    run_order = {}
    counter = 0
    with open("../test/runlist_15.txt", 'r') as run_file:
        run_data = run_file.readlines()
    for idx, datum in enumerate(run_data):
        run_order[idx] = datum.split()

    for idx in tqdm(run_order):
        object_name = run_order[idx][0]
        seqid = run_order[idx][1]
        if f"nu{seqid}A01_cl.evt" in os.listdir(f"../bifrost_data/15/{seqid}/event_cl/"):
            if f"nu{seqid}B01_cl.evt" in os.listdir(f"../bifrost_data/15/{seqid}/event_cl/"):
                if f"{dtime}_binning_flag.txt" not in os.listdir(f"../bifrost_data/15/{seqid}/event_cl/"):
                    path = f"../bifrost_data/15/{seqid}/"
                    evdir = f"{path}event_cl/"
                    out_path = f"{path}products/"
                    run_object = NuAnalysis(dtime, 3, path=path, evdir=evdir, seqid=seqid, out_path=out_path, clean=True, bifrost=True, object_name=object_name, sessionid=15)
                    run_object.event_extraction()
                    os.chdir(starting_directory)
