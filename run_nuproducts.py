"""
Script dedicated to calling different variations of the ``nuproducts`` routine.
"""
from nuanalysis import *
import os

# Setup appropriate pathing and diagnostic information
seqid = ""
path = ""
outpath = ""

# Select analysis parameters
instrument = ""
dt = 10
# Should be provided in degrees, for the moment to be done manually
src_position = (1 , 1)
bkg_position = (1, 1)

os.chdir("")

script_string = f"nuproducts indir=./event_cl instrument={instrument} stemimputs=nu{seqid} outdir={outpath} binsize={dt}"


#if gti:
#    script_string += " gtifile={}"
#os.chdir("")