""" 
This is a testing script for combining different region files for use 
in straycats.py. This is done using the shapely geometry package since 
the python `regions` package has not yet supported operations to 
combine regions.
"""
from regions import Regions
from shapely.geometry import Point
from shapely.plotting import plot_polygon
import matplotlib.pyplot as plt


# Read in a regions file and record its contents
filepath = "/Volumes/data_ssd_1/bifrost_data/straycat/40014024001A_StrayCatsI_260.reg"
region_list = Regions.read(filepath, format='ds9')
with open(filepath, 'r') as region_file:
    reg_info = region_file.readlines()


# Initialize sequence of operations
construct = []
flag = False
for line in reg_info:
    if flag:
        if line[0] == '-':
            construct.append('-')
        else:
            construct.append('+')
    if line == "image\n":
        flag = True


# If only one shape is in region file, create a simple circular 
# mask to represent the region file.
if len(region_list) == 1:
    pos = region_list[0].center.xy
    radius = region_list[0].radius
    src_x = pos[0]
    src_y = pos[1]
    mask = Point(src_x, src_y).buffer(radius, resolution=1000)


# Routine for combining multiple regions using the shapely package
if len(region_list) > 1:
    pos = region_list[0].center.xy
    radius = region_list[0].radius
    src_x = pos[0]
    src_y = pos[1]
    mask = Point(src_x, src_y).buffer(radius, resolution=1000)
    
    for idx, region in enumerate(region_list):
        if idx == 0:
            continue
        
        # Combine two regions via a union 
        if construct[idx] == '+':
            pos = region.center.xy
            radius = region.radius
            src_x = pos[0]
            src_y = pos[1]
            mask_add = Point(src_x, src_y).buffer(radius, resolution=1000)
            mask = mask.union(mask_add)
        
        # Subtract two regions via a difference
        if construct[idx] == '-':
            pos = region.center.xy
            radius = region.radius
            src_x = pos[0]
            src_y = pos[1]
            mask_subtract = Point(src_x, src_y).buffer(radius, resolution=1000)
            mask = mask.difference(mask_subtract)

# Create a test plot to show the polygon in DET1 coordinates.
fig, ax = plt.subplots()
plot_polygon(mask, ax=ax, add_points=False)
plt.xlim(0, 300)
plt.ylim(0, 300)
plt.show()

# Test line for seeing if points are located within the plotted 
# StrayCat region.
test_point = Point(22, ).buffer(1, resolution=1000)
print(mask.contains(test_point))