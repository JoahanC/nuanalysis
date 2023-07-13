import glob
import subprocess

def recover_events():
    download_files = glob.glob("../../../../Downloads/*.tgz")
    for file in download_files:
        obsid = file.replace("../../../../Downloads/", '')[:11]
        current_files = glob.glob("*/")
        if len(current_files) == 0:
            subprocess.run(["cp", file, "."])
            subprocess.run(["tar", "-xf", file.replace("../../../../Downloads/", '')])
            subprocess.run(["rm", file.replace("../../../../Downloads/", '')])
        current_files = glob.glob("*/")
        present_flag = False
        for cur_file in current_files:
            if obsid in cur_file:
                present_flag = True 
        if not present_flag:
            subprocess.run(["cp", file, "."])
            subprocess.run(["tar", "-xf", file.replace("../../../../Downloads/", '')])
            subprocess.run(["rm", file.replace("../../../../Downloads/", '')])