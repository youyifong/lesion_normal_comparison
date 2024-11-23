"""
working dir: folder containing ROI files
"""

from shapely.geometry import Polygon, Point

import glob,argparse,os
from read_roi import read_roi_file
import pandas as pd
import numpy as np

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--XX', type=str, default='CF_Les', help='tissue label')

args = parser.parse_known_args()[0]
print(args, '\n')


# set tissue label here if not set in command line
if args.XX is None:
    args.XX=os.environ.get('XX')
print(args.XX)


# read whole tissue boundary ROI file
tissue_boundary_roi = [f for f in glob.glob('*.roi') if "whole" in f]
print(tissue_boundary_roi)
if len(tissue_boundary_roi)!=1:
    raise RuntimeError("STOP: there should be one and only whole tissue roi file !!!!!!!!!!!!!!!!!!!!!!!")
tmp = read_roi_file(tissue_boundary_roi[0])
val = list(tmp.values()) [0]
whole_tissue = Polygon(list(zip(val["x"], val["y"])))


###############################################################################
# make XX_region_area_length file

# compute the area of each region
region_roi_files = [f for f in glob.glob('*.roi') if "epi_region" in f or "dermis_region" in f]
region_roi_files.sort()
print(region_roi_files)
areas=[]
for region_roi_file in region_roi_files:
    region = read_roi_file(region_roi_file)    
    val = list(region.values()) [0]
    p1 = Polygon(list(zip(val["x"], val["y"])))
    if whole_tissue is None:
        areas.append(p1.area)
    else:
        areas.append(p1.intersection(whole_tissue).area)



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


# save the area and length of each region to a file

# remove .roi from names
region_roi_files=[f.replace(".roi","") for f in region_roi_files]
length_roi_files=[f.replace(".roi","") for f in length_roi_files]

if args.XX == 'CJ_Les':
    # add the two line segments length to get tissue length 
    data = pd.DataFrame({"name": region_roi_files+["CJ_Les_length_1_normal"], 'size': areas+[lengths[0]+lengths[1]]})
else: 
    data = pd.DataFrame({"name": region_roi_files+length_roi_files, 'size': areas+lengths})
print(data)

data.to_csv("../"+args.XX+"_region_area_length.csv", header=True, index=False, sep=',')


###############################################################################
# Create a Grid Within a Polygon

# Define the polygon ROI (example: a triangle)
polygon_coords = [(0, 0), (10, 0), (5, 10)]
polygon = Polygon(polygon_coords)

# Define grid spacing
grid_spacing = 1.0

# Get bounds of the polygon
min_x, min_y, max_x, max_y = polygon.bounds

# Generate grid points
x_coords = np.arange(min_x, max_x, grid_spacing)
y_coords = np.arange(min_y, max_y, grid_spacing)
grid_points = [(x, y) for x in x_coords for y in y_coords]

# Filter points that are inside the polygon
points_in_polygon = [Point(p) for p in grid_points if polygon.contains(Point(p))]

# Display results
print("Number of points inside polygon:", len(points_in_polygon))
for point in points_in_polygon:
    print(point)
