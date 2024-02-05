import glob
from nustar import *
from helpers import *
from astropy.table import Table, vstack
from tqdm import tqdm


data_cols = ['INDEX', 'COUNTS', 'XPIX', 'YPIX', 'VIGCOR', 'ERR', 'HBOX', 'PROB', 'SNR', 'TSTART', 'TSTOP', 'COUNTSERR', 'SEP', 'DET1X', 'DET1Y']
string_cols = ['RA', 'DEC', 'SEQID', 'BOUND']
basepath = "/Volumes/data_ssd_1/bifrost_data/**/**/detections/*500_basic_det1.tbl"
table_files = glob.glob(basepath)
total_table = Table.read(table_files[0], format='ipac')
for file in tqdm(table_files[1:]):
    data_table = Table.read(file, format='ipac')
    for col in data_cols:
        data_table[col] = data_table[col].astype(float)
        total_table[col] = total_table[col].astype(float)
    for scol in string_cols:
        data_table[scol] = data_table[scol].astype(str)
        total_table[scol] = total_table[scol].astype(str)
    if len(data_table['RA']) == 0:
        continue
    total_table = vstack([total_table, data_table])

total_table.write("ssd1_jan9_500_final.tbl", format='ipac', overwrite=True)