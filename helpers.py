import glob
import os
import shutil
import subprocess


def generate_directory(path, overwrite=False):
    if os.path.isdir(path):
        if overwrite:
            clear_directory(path)
        else:
            print("Please enable overwriting to remove this directory!")
    else:
        os.mkdir(path)


def clear_directory(path):
    if os.path.isdir(path):
        for file in os.listdir(path):
            os.remove(f"{path}/{file}")
    else:
        NotADirectoryError("This directory doesn't exist!")

        
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