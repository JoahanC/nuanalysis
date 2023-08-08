import glob
import os
import subprocess
import random
from utils import energy_to_chan


def generate_directory(path, overwrite=False):
    if os.path.isdir(path):
        if overwrite:
            clear_directory(path)
        #else:
        #    print("Please enable overwriting to remove this directory!")
    else:
        os.mkdir(path)


def clear_directory(path):
    if os.path.isdir(path):
        print(os.listdir(path))
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


def read_detection_file(det_file):   
    with open(det_file) as detections:
        for i in range(14):
            detections.readline()
        detect_info = {}
        detect_info["INDEX"] = []
        detect_info["COUNTS"] = []
        detect_info["XPIX"] = []
        detect_info["YPIX"] = []
        detect_info["VIGCOR"] = []
        detect_info["RA"] = []
        detect_info["DEC"] = []
        detect_info["ERR"] = []
        detect_info["HBOX"] = []
        detect_info["PROB"] = []
        detect_info["SNR"] = []

        for line in detections:
            line_info = line.split()
            detect_info["INDEX"].append(line_info[0])
            detect_info["COUNTS"].append(line_info[1])
            detect_info["XPIX"].append(line_info[2])
            detect_info["YPIX"].append(line_info[3])
            detect_info["VIGCOR"].append(line_info[4])
            detect_info["RA"].append(f"{line_info[5]} {line_info[6]} {line_info[7]}")
            detect_info["DEC"].append(f"{line_info[8]} {line_info[9]} {line_info[10]}")
            detect_info["ERR"].append(line_info[11])
            detect_info["HBOX"].append(line_info[12])
            detect_info["PROB"].append(line_info[13])
            detect_info["SNR"].append(line_info[14])
    return detect_info


def read_detection_dir(path):
    det_files = glob.glob(path + "/*.det")
    detect_info = {}
    for file in det_files:
        with open(file) as detections:
            for i in range(14):
                detections.readline()
            times = tuple(file.replace("./../data/30501002002/detections/559-1934_10000-3/", '').replace(".det", '').split(sep="-"))
            detect_info[times] = {}
            detect_info[times]["INDEX"] = []
            detect_info[times]["COUNTS"] = []
            detect_info[times]["XPIX"] = []
            detect_info[times]["YPIX"] = []
            detect_info[times]["VIGCOR"] = []
            detect_info[times]["RA"] = []
            detect_info[times]["DEC"] = []
            detect_info[times]["ERR"] = []
            detect_info[times]["HBOX"] = []
            detect_info[times]["PROB"] = []
            detect_info[times]["SNR"] = []
            detect_info[times]["FILE"] = []

            for line in detections:
                line_info = line.split()
                detect_info[times]["FILE"].append(file)
                detect_info[times]["INDEX"].append(line_info[0])
                detect_info[times]["COUNTS"].append(line_info[1])
                detect_info[times]["XPIX"].append(line_info[2])
                detect_info[times]["YPIX"].append(line_info[3])
                detect_info[times]["VIGCOR"].append(line_info[4])
                detect_info[times]["RA"].append(f"{line_info[5]} {line_info[6]} {line_info[7]}")
                detect_info[times]["DEC"].append(f"{line_info[8]} {line_info[9]} {line_info[10]}")
                detect_info[times]["ERR"].append(line_info[11])
                detect_info[times]["HBOX"].append(line_info[12])
                detect_info[times]["PROB"].append(line_info[13])
                detect_info[times]["SNR"].append(line_info[14])

    return detect_info


def generate_gti_files(path, tstart, tstop):
    """
    Creates a Good Time Interval (GTI) file for temporal filtering calls using XSelect and 
    nuproducts.

    Parameters
    ----------
    path : str
        The path at which to execute the XSelect script
    tstart : str or float
        The starting time of the interval in SCC
    tstop : str or float
        The stopping time of the interval in SCC
    """
    with open(path + "xselect.xco", 'w') as script:
        script.write('\n')
        script.write('\n')
        script.write("filter time scc\n")
        script.write(f"{tstart.replace('nu_', '')} , {tstop}\n")
        script.write("x\n")
        script.write("save time keyboard\n")
        script.write(f"{tstart}_{tstop}_gti.fits\n")
        script.write("exit no")
        script.write('\n')
    subprocess.run(["xselect", "@xselect.xco"], cwd=path, capture_output=False)
    subprocess.run(["rm", "xselect.xco"], cwd=path, capture_output=False)


def chan_to_energy(chan):
    return chan * 0.04 + 1.6


def remove_tmp_folders(path):
    tmp_dirs = glob.glob(path + "*tmp_nuproducts")
    for dir in tmp_dirs:
        subprocess.run(["rm", "-rf", dir])


def make_xselect_commands(infile, outfile, dir, elow, ehigh, sessionid, usrgti=False, evt_extract=True):
    '''
    Helper script to generate the xselect commands to make an image in a given NuSTAR range
    '''
    n1 = random.random()
    n2 = random.random()
    n3 = random.random()
    n4 = random.choice([100, 1000, 10000])
    number = n1*n2*n3*n4
    xsel=open(f"{number}xsel.xco","w")
    xsel.write(f'{number}\n')
    xsel.write('\n')
    xsel.write(f"{sessionid}\n")
    xsel.write("read events \n")
    xsel.write(f'{dir}\n')
    xsel.write(f'{infile}\n ')
    xsel.write('yes \n')
    pi_low = energy_to_chan(elow)
    pi_high = energy_to_chan(ehigh)
    if usrgti is not False:
        xsel.write(f'filter time \n')
        xsel.write('file \n')
        xsel.write(f'{usrgti}\n')
        xsel.write('extract events\n')
    xsel.write('filter pha_cutoff {} {} \n'.format(pi_low, pi_high))

    xsel.write('set xybinsize 1\n')
    xsel.write("extract image\n")
    xsel.write("save image\n")
    xsel.write("%s \n" % outfile)
    if evt_extract:
        xsel.write("extract events")
        xsel.write("\n")
        xsel.write(f"save events {outfile.replace('.fits', '.evt')}\n")
        xsel.write('no\n')
    xsel.write('exit\n')           
    xsel.write('n \n')
    xsel.close()
    return 'xsel.xco', number