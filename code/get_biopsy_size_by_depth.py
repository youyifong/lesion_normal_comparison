"""
Compute the size of biopsy, defined as the number of points in a 5x5 px grid, by depth in 200 px (equivalently, 100 um) increments for all biopsies.
Save the results in Chipcytometry_images/biopsy_size_by_depth.csv.
"""



from shapely.geometry import Polygon, Point
import glob,os
from read_roi import read_roi_file
import pandas as pd
import numpy as np
from tsp.geom import pnt2line


def do_biopsy(roi_dir):
    """
    read whole tissue boundary ROI file and surface ROI file
    """

    biopsy_roi = [f for f in glob.glob(roi_dir+'*.roi') if "whole" in f]
    print(biopsy_roi)
    if len(biopsy_roi)!=1:
        raise RuntimeError("There should be one and only one whole tissue roi file.")

    val = list(read_roi_file(biopsy_roi[0]).values()) [0]
    biopsy = Polygon(list(zip(val["x"], val["y"])))


    skin_roi = [f for f in glob.glob(roi_dir+'*.roi') if "skin_boundary" in f]
    print(skin_roi)
    if len(skin_roi)!=1:
        raise RuntimeError("There should be one and only one skin boundary roi file.")

    val = list(read_roi_file(skin_roi[0]).values()) [0]
    surface = np.vstack((val["x"], val["y"])).T
    n = len(val["x"])



    """
    Create a Grid Within the Biopsy
    """

    # Define grid spacing
    grid_spacing = 5.0

    # Get bounds of the biopsy
    min_x, min_y, max_x, max_y = biopsy.bounds

    # Generate grid points
    x_coords = np.arange(min_x, max_x, grid_spacing)
    y_coords = np.arange(min_y, max_y, grid_spacing)
    grid_points = [(x, y) for x in x_coords for y in y_coords]

    # Filter points that are inside the biopsy
    points_in_polygon = [list(p) for p in grid_points if biopsy.contains(Point(p))]

    # Display results
    print("Number of points inside biopsy:", len(points_in_polygon))


    """
    Count the number of grid points by depth in 200 px (equivalently, 100 um) increments. 
    """

    # compute depth as measured by the shortest distance to the surface
    depth = []
    for point in points_in_polygon:
        depth.append(np.min(pnt2line(P0=surface[0:(n - 1), :], P1=surface[1:n, :], p=point)))


    # bin the depth into 200 px increments
    depth = np.array(depth)
    depth_bins = np.arange(0, np.ceil(np.max(depth)), 200)

    return np.histogram(depth, bins=depth_bins)


# list the directories under ../Chipcytometry_images
chipcytometry_images_dir = '../Chipcytometry_images/'
directories = [d for d in os.listdir(chipcytometry_images_dir) if os.path.isdir(os.path.join(chipcytometry_images_dir, d))]

res=[]
for folder in directories:
    tmp = chipcytometry_images_dir + folder + '/ZKWstitched_alignedtoCD45_cropped/ROI_regions/'
    if os.path.exists(tmp):
        print(folder)
        dist, depth_hist = do_biopsy(tmp)

        # add results to the dataframe. convert dist from px to um by dividing by 2
        for i in range(len(dist)):
            res.append({'biopsy': folder, 'dist': dist[i]/2, 'depth_hist': depth_hist[i]})

# save the results to a csv file
df = pd.DataFrame(res)
df.to_csv(chipcytometry_images_dir+'/tissue_size_by_depth_from_skin_for_all_biopsies.csv', index=False)