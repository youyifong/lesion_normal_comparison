# compute area of intersection of polygons

from shapely.geometry import Polygon
import glob
from read_roi import read_roi_file 
import pandas as pd
import numpy as np

s="DD"


tmp = read_roi_file(s+"_Les_ROI_regions/"+s+"_Les_whole_tissue.roi")    
val = list(tmp.values()) [0]
tissue = Polygon(list(zip(val["x"], val["y"])))

region_roi_files = glob.glob(s+'_Les_ROI_regions/'+s+'_Les_*_region_*.roi')
region_roi_files.sort()
print(region_roi_files)

areas=[]
for region_roi_file in region_roi_files:
    region = read_roi_file(region_roi_file)    
    val = list(region.values()) [0]
    p1 = Polygon(list(zip(val["x"], val["y"])))
    areas.append(p1.intersection(tissue).area)

data = pd.DataFrame({"name": region_roi_files, 'size': areas})

data.to_csv(s+"_region_size.csv", header=True, index=None, sep=',')



# compute the length of the smooth curve to measure the size of the lesion

tmp = read_roi_file(s+"_Les_ROI_regions/"+s+"_Les_lesion_skin_boundary_smooth.roi")    
val = list(tmp.values()) [0]
points = np.array(list(zip(val["x"], val["y"])))

# Calculate the length of each segment
segment_lengths = np.linalg.norm(points[:-1] - points[1:], axis=1)

# Sum the lengths of all segments to get the length of the polyline
polyline_length = np.sum(segment_lengths)


data.to_csv(s+"_region_size.csv", header=True, index=None, sep=',')


