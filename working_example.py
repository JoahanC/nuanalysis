""" 
This file demonstrates how to utilize the main class structure in this repository 
and provides some commentary on how the different methods work and caveats with 
using them.
"""
from nuanalysis import *
import os

# Initialize starting working directory
starting_directory = os.getcwd()

# Set the sequence ids and characteristic times for processing.
seqid = "60160619002"
dtime = 5000

# Initialize reference files, these can be set to an absolute path elsewhere
low_phi_file = "./ref_files/nustar_pilow.txt"
high_phi_file = "./ref_files/nustar_pihi.txt"

# Initialize path parameters
object_name = "test" # Not necessary but nice for bookkeeping
path = f"./test_data/{seqid}/" # Sequence id folder
evdir = f"{path}event_cl/" # Directory containing all clean events files
out_path = f"{path}products/" # Directory to contain all product files

# Uncomment different methods to test functionality. They have been provided in 
# the order they should be run

# INITIALIZE AND CLEAN SEQUENCE ID
# What this line will do is it will clean out a sequence id directory and 
# run the NuSTAR pipeline on the data contained within the directory. It 
# makes use of aggressive settings to minimize background counts
run_object = NuAnalysis(dtime, 3, path=path, low_phi_file=low_phi_file, high_phi_file=high_phi_file,
                                    seqid=seqid, clean=True, bifrost=True, object_name=object_name)

##### RUN DETECTION STEPS #####

# OPTIONAL STEP
# If the astrometry of the target source is poor, i.e., the source is actually 
# not at the header coordinates, run this command to recalculate the target location 
# so that the main source is properly masked
# The test source provided for this script actually suffers from this issue, but commenting 
# this out allows for an easy demonstration of the further analysis methods.

#run_object.recalculate_center()

# RUNNING THE SLIDING CELL SEARCH
# This method calls three submethods and comprises the bulk of the transient 
# searching step that is built into the class. It begins by first extracting 
# every event in this NuSTAR observation and rewritting it into seperate 
# FITS files which correspond to the time bins set by the `dtime` parameter
# It then performs the sliding cell detection step whereby it utilizes the 
# `ximage` sliding cell search algorithm to identify local peaks in each 
# binned fits image. It concludes by merging the detections together to 
# attempt to eliminate duplicates and leave uniquely identified transient 
# candidates

#run_object.run_detection_pipeline()

# AGGREGATING THE PRELIMINARY CANDIDATE LOCATIONS
# This method then aggregates all of the detections and writes them to a file for 
# future use and analysis. The file written to is located in the detections folder 
# in the observation directory: `5000_detections.tbl`

#run_object.write_net_detections()

# POST PROCESSING 
# This method recalculates the poisson statistic onto each candidate location to 
# address the background value issues that `ximage` may have. This step also takes 
# each potential candidate and converts its position into DET1 coordinates. As a result, 
# this method can take some time to run since the calculation of DET1 coordinates is 
# relatively slow

# The newly written table should be called `5000_poisson.tbl`

#run_object.recalculate_poisson()

# BRIGHTNESS FILTERING
# This method provides a crude metric for filtering bright sources. This calculates the 
# highest photon count present in any of the given time bins within 30 arcseconds of the 
# center target. 

#max_count = run_object.determine_bright() 
#print(max_count)
os.chdir(starting_directory)