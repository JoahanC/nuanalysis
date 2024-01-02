import glob
from nustar import *
from helpers import *
from astropy.table import Table, vstack
from tqdm import tqdm


"""cols = ['INDEX', 'COUNTS', 'XPIX', 'YPIX', 'VIGCOR', 'ERR', 'HBOX', 'PROB', 'SNR', 'TSTART', 'TSTOP', 'COUNTSERR', 'SEP']#, 'DET1X', 'DET1Y']
scols = ['RA', 'DEC', 'SEQID', 'BOUND']
basepath = "/Volumes/data_ssd_2/**/**/detections/*500_detections.tbl"
table_files = glob.glob(basepath)
total_table = Table.read(table_files[0], format='ipac')
for file in tqdm(table_files[1:]):
    data_table = Table.read(file, format='ipac')
    for col in cols:
        data_table[col] = data_table[col].astype(float)
        total_table[col] = total_table[col].astype(float)
    for scol in scols:
        data_table[scol] = data_table[scol].astype(str)
        total_table[scol] = total_table[scol].astype(str)
    if len(data_table['RA']) == 0:
        continue
    total_table = vstack([total_table, data_table])
print(len(total_table['RA']))
total_table.write("ssd2_basic_500_final.tbl", format='ipac', overwrite=True)"""

basepath = "/Volumes/data_ssd_1/bifrost_data/**/**/detections/*_detections.tbl"
table_files = glob.glob(basepath)
high_detection_seqids = []
for file in tqdm(table_files[1:]):
    data_table = Table.read(file, format='ipac')
    n_det = len(data_table['RA'])
    if n_det > 400:
        seqid = file.replace("/Volumes/data_ssd_1/bifrost_data/", '')
        run_cycle = int(seqid[:2].replace('/', ''))
        if run_cycle < 10:
            seqid = seqid[seqid.find('/') + 1:13]
        if run_cycle >= 10:
            seqid = seqid[seqid.find('/')+ 1:14]
        #print(file)
        print(seqid)
        #dtime = file[42:]
        dtime = file[54:]
        dtime = dtime[dtime.find('/') + 1:].replace("_detections.tbl", '')
        high_detection_seqids.append([seqid, dtime, n_det])

print(high_detection_seqids)
with open("./pathfiles/high_detection_1.txt", 'w') as high_det:
    for seq in high_detection_seqids:
        high_det.write(f"{seq[0]} {seq[1]}\n")

high_seqids = []
with open("./pathfiles/high_detection_1.txt") as high_det:
    for line in high_det.readlines():
        high_seqids.append(line.replace('\n', ''))