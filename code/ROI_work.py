'''
working dir: folder containing ROI files
'''

from shapely.geometry import Polygon
import glob,argparse,os
from read_roi import read_roi_file 
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(description='parameters')
parser.add_argument('--XX', type=str, help='tissue label')    
args = parser.parse_args()

# set tissue label here if not set in command line
if args.XX is None:
    args.XX=os.environ.get('XX')
print(args.XX)



############################################################################
# compute size of regions

region_roi_files = [f for f in glob.glob('*.roi') if "epi_region" in f or "dermis_region" in f]
region_roi_files.sort()
print(region_roi_files)

# whole tissue boundary ROI file
tissue_boundary_roi = [f for f in glob.glob('*.roi') if "whole" in f]
print(tissue_boundary_roi)
if len(tissue_boundary_roi)!=1:
    print("STOP: there should only be one whole tissue roi file !!!!!!!!!!!!!!!!!!!!!!!")
try:
    tmp = read_roi_file(tissue_boundary_roi[0])    
    val = list(tmp.values()) [0]
    whole_tissue = Polygon(list(zip(val["x"], val["y"])))
except Exception: # as e
    # Catch any other type of exception
    print("Something went wrong reading tissue boundary roi")
    # print("Error:", str(e))
    whole_tissue = None


# compute the area of each region
areas=[]
for region_roi_file in region_roi_files:
    region = read_roi_file(region_roi_file)    
    val = list(region.values()) [0]
    p1 = Polygon(list(zip(val["x"], val["y"])))
    if whole_tissue is None:
        areas.append(p1.area)
    else:
        areas.append(p1.intersection(whole_tissue).area)


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

# remove .roi from names
region_roi_files=[f.replace(".roi","") for f in region_roi_files]

if args.XX == 'CJ':
    # remove _v3 from name
    region_roi_files=[f.replace("_v3","") for f in region_roi_files]
    # add the two line segments length to get tissue length 
    data = pd.DataFrame({"name": region_roi_files+["CJ_Les_length_1_normal"], 'size': areas+[lengths[0]+lengths[1]]})
else: 
    data = pd.DataFrame({"name": region_roi_files+length_roi_files, 'size': areas+lengths})
print(data)

data.to_csv("../"+args.XX+"_Les_region_area_length.csv", header=True, index=None, sep=',')
