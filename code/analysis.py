'''
working dir: folder containing ROI files
'''

from shapely.geometry import Polygon
import glob
from read_roi import read_roi_file 
import pandas as pd
import numpy as np



############################################################################
# compute region area

region_roi_files = [f for f in glob.glob('*.roi') if "epi" in f or "dermis" in f]
region_roi_files.sort()
print(region_roi_files)

# whole tissue boundary ROI file
# tissue_boundary_roi = [f for f in glob.glob('*.roi') if "whole" in f]
# tmp = read_roi_file(tissue_boundary_roi)    
# val = list(tmp.values()) [0]
# whole_tissue = Polygon(list(zip(val["x"], val["y"])))

# compute the area of each region
areas=[]
for region_roi_file in region_roi_files:
    region = read_roi_file(region_roi_file)    
    val = list(region.values()) [0]
    p1 = Polygon(list(zip(val["x"], val["y"])))
    areas.append(p1.area)
    # areas.append(p1.intersection(whole_tissue).area)


############################################################################
# compute the length of the lines as measure of lesion or normal tissue size

length_roi_files = [f for f in glob.glob('*.roi') if "length" in f]
length_roi_files.sort()
print(length_roi_files)

lengths=[]
for length_roi_file in length_roi_files:
    tmp = read_roi_file(length_roi_file) 
    val = list(tmp.values()) [0]
    
    if val["type"]=="line":
        points = np.array(list(zip([val["x1"], val["x2"]], [val["y1"], val["y2"]])))
    else:
        points = np.array(list(zip(val["x"], val["y"])))
    
    # Calculate the length of each segment
    segment_lengths = np.linalg.norm(points[:-1] - points[1:], axis=1)
    
    # Sum the lengths of all segments to get the length of the polyline
    polyline_length = np.sum(segment_lengths)
    lengths.append(polyline_length)



###############################################################################
# save

data = pd.DataFrame({"name": region_roi_files+length_roi_files, 'size': areas+lengths})
data
data.to_csv("../region_area_length.csv", header=True, index=None, sep=',')




# # lesion
# tmp = read_roi_file(s+"_Les_lesion_line.roi")    
# val = list(tmp.values()) [0]

# if val["type"]=="line":
#     points = np.array(list(zip([val["x1"], val["x2"]], [val["y1"], val["y2"]])))
# else:
#     points = np.array(list(zip(val["x"], val["y"])))

# # Calculate the length of each segment
# segment_lengths = np.linalg.norm(points[:-1] - points[1:], axis=1)

# # Sum the lengths of all segments to get the length of the polyline
# polyline_length = np.sum(segment_lengths)
# polyline_length


# # normal
# segment_lengths=[]
# for i in range(1,2): # for DD
# #for i in range(1,3): # for CF
#     tmp = read_roi_file(s+f"_Les_normal_line_{i}.roi")    
#     val = list(tmp.values()) [0]
    
#     if val["type"]=="line":
#         points = np.array(list(zip([val["x1"], val["x2"]], [val["y1"], val["y2"]])))
#     else:
#         points = np.array(list(zip(val["x"], val["y"])))
    
#     # Calculate the length of each segment
#     segment_lengths.append(np.linalg.norm(points[:-1] - points[1:], axis=1))

# # Sum the lengths of all segments to get the length of the polyline
# np.sum(segment_lengths)



