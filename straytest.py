from regions import Regions
from shapely.geometry import Point, Polygon
from shapely.plotting import plot_polygon

import matplotlib.pyplot as plt


region_list = Regions.read("/Volumes/data_ssd_1/bifrost_data/straycat/40014024001A_StrayCatsI_260.reg", format='ds9')
with open("/Volumes/data_ssd_1/bifrost_data/straycat/40014024001A_StrayCatsI_260.reg", 'r') as filey:
    reg_info = filey.readlines()
    
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

if len(region_list) == 1:
    pos = region_list[0].center.xy
    radius = region_list[0].radius
    src_x = pos[0]
    src_y = pos[1]
    mask = Point(src_x, src_y).buffer(radius, resolution=1000)
    
if len(region_list) > 1:
    pos = region_list[0].center.xy
    radius = region_list[0].radius
    src_x = pos[0]
    src_y = pos[1]
    mask = Point(src_x, src_y).buffer(radius, resolution=1000)
    
    for idx, region in enumerate(region_list):
        if idx == 0:
            continue 
        if construct[idx] == '+':
            pos = region.center.xy
            radius = region.radius
            src_x = pos[0]
            src_y = pos[1]
            mask_add = Point(src_x, src_y).buffer(radius, resolution=1000)
            mask = mask.union(mask_add)
            
        if construct[idx] == '-':
            pos = region.center.xy
            radius = region.radius
            src_x = pos[0]
            src_y = pos[1]
            mask_subtract = Point(src_x, src_y).buffer(radius, resolution=1000)
            mask = mask.difference(mask_subtract)

fig, ax = plt.subplots()
plot_polygon(mask, ax=ax, add_points=False)
plt.xlim(0, 300)
plt.ylim(0, 300)
plt.show()

test_point = Point(22, ).buffer(1, resolution=1000)
print(mask.contains(test_point))