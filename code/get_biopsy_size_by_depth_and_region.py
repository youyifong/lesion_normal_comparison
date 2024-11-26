"""
Compute the size of biopsy, defined as the number of points in a 5x5 px grid, by depth in 200 px (equivalently, 100 um) increments for all biopsies.
Save the results in Chipcytometry_images_from8bit/biopsy_size_by_depth_and_region.csv.
"""



from shapely.geometry import Polygon, Point
import glob,os
from read_roi import read_roi_file
import pandas as pd
import numpy as np
from tsp.geom import pnt2line


def do_biopsy(roi_dir):
    """
    read whole tissue boundary ROI file and skin ROI file
    """

    skin_roi = glob.glob(roi_dir+'*skin_boundary*.roi')
    print(skin_roi)
    if len(skin_roi)!=1:
        raise RuntimeError("There should be one and only one skin boundary roi file.")
    val = list(read_roi_file(skin_roi[0]).values()) [0]
    skin = np.vstack((val["x"], val["y"])).T
    n_skin = len(val["x"])

    DEJ_roi = glob.glob(roi_dir+'*DEJ*.roi')
    print(DEJ_roi)
    if len(DEJ_roi)!=1:
        raise RuntimeError("There should be one and only one skin boundary roi file.")
    val = list(read_roi_file(DEJ_roi[0]).values()) [0]
    DEJ = np.vstack((val["x"], val["y"])).T
    n_DEJ = len(val["x"])

    biopsy_roi = glob.glob(roi_dir+'*whole*.roi')
    print(biopsy_roi)
    if len(biopsy_roi)!=1:
        raise RuntimeError("There should be one and only one whole tissue roi file.")
    val = list(read_roi_file(biopsy_roi[0]).values()) [0]
    biopsy = Polygon(list(zip(val["x"], val["y"])))

    normal_dermis_regions_roi = glob.glob(roi_dir+'*_dermis_region_*normal*.roi')
    print("normal_dermis_regions_rois")
    print(normal_dermis_regions_roi)
    normal_dermis_regions = []
    for x in normal_dermis_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        normal_dermis_regions.append (Polygon(list(zip(val["x"], val["y"]))))

    normal_epi_regions_roi = glob.glob(roi_dir+'*_epi_region_*normal*.roi')
    print("normal_epi_regions_rois")
    print(normal_epi_regions_roi)
    normal_epi_regions = []
    for x in normal_epi_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        normal_epi_regions.append (Polygon(list(zip(val["x"], val["y"]))))

    lesion_dermis_regions_roi = glob.glob(roi_dir+'*_dermis_region_*lesion*.roi')
    print("lesion_dermis_regions_rois")
    print(lesion_dermis_regions_roi)
    lesion_dermis_regions = []
    for x in lesion_dermis_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        lesion_dermis_regions.append (Polygon(list(zip(val["x"], val["y"]))))

    lesion_epi_regions_roi = glob.glob(roi_dir+'*_epi_region_*lesion*.roi')
    print("lesion_epi_regions_rois")
    print(lesion_epi_regions_roi)
    lesion_epi_regions = []
    for x in lesion_epi_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        lesion_epi_regions.append (Polygon(list(zip(val["x"], val["y"]))))


    """
    Create a grid of points and classify each point as normal or lesion based on the region it is in.
    """

    # Define grid spacing
    grid_spacing = 5.0
    depth_spacing = 200.0

    # Get bounds of the biopsy
    min_x, min_y, max_x, max_y = biopsy.bounds

    # Generate grid points
    x_coords = np.arange(min_x, max_x, grid_spacing)
    y_coords = np.arange(min_y, max_y, grid_spacing)
    grid_points = [(x, y) for x in x_coords for y in y_coords]

    points_in_normal_epi=[]
    points_in_lesion_epi=[]
    points_in_normal_dermis=[]
    points_in_lesion_dermis=[]

    for p in grid_points:
        classified = False

        for polygon in normal_epi_regions:
            if polygon.contains(Point(p)):
                points_in_normal_epi.append(p)
                classified = True
                break
        if classified:
            continue

        for polygon in normal_dermis_regions:
            if polygon.contains(Point(p)):
                points_in_normal_dermis.append(p)
                classified = True
                break
        if classified:
            continue

        for polygon in lesion_dermis_regions:
            if polygon.contains(Point(p)):
                points_in_lesion_dermis.append(p)
                classified = True
                break
        if classified:
            continue

        for polygon in lesion_epi_regions:
            if polygon.contains(Point(p)):
                points_in_lesion_epi.append(p)
                break

    points_in_normal = points_in_normal_epi + points_in_normal_dermis
    points_in_lesion = points_in_lesion_epi + points_in_lesion_dermis


    """
    Count the number of grid points by distance to skin in depth_spacing px increments. 
    """

    # normal regions: compute depth as measured by the shortest distance to the skin for each point in the normal regions
    depth1 = []
    for point in points_in_normal:
        depth1.append(np.min(pnt2line(P0=skin[0:(n_skin - 1), :], P1=skin[1:n_skin, :], p=point)))

    # lesion regions: compute depth as measured by the shortest distance to the skin for each point in the lesion regions
    depth2 = []
    for point in points_in_lesion:
        depth2.append(np.min(pnt2line(P0=skin[0:(n_skin - 1), :], P1=skin[1:n_skin, :], p=point)))

    # bin all depths into depth_spacing px increments
    depths = np.array(depth1+depth2)
    depth_bins = np.arange(0,
                           np.ceil(np.max(depths) / depth_spacing) * depth_spacing + depth_spacing,
                           depth_spacing)

    # compute histograms
    skin_hist1 = np.histogram(np.array(depth1), bins=depth_bins)
    skin_hist2 = np.histogram(np.array(depth2), bins=depth_bins)
    skin_hist0 = np.histogram(depths, bins=depth_bins)
    
    
    """
    Count the number of grid points by distance to DEJ in depth_spacing px increments. 
    """
    
    # normal regions: compute depth as measured by the shortest distance to the DEJ for each point in the normal regions
    depth1 = []
    # epidermis cells count as negative depth
    for point in points_in_normal_epi:
        depth1.append( - np.min(pnt2line(P0=DEJ[0:(n_DEJ - 1), :], P1=DEJ[1:n_DEJ, :], p=point)))
    for point in points_in_normal_dermis:
        depth1.append(np.min(pnt2line(P0=DEJ[0:(n_DEJ - 1), :], P1=DEJ[1:n_DEJ, :], p=point)))

    # lesion regions: compute depth as measured by the shortest distance to the DEJ for each point in the lesion regions
    depth2 = []
    # epidermis cells count as negative depth
    for point in points_in_lesion_epi:
        depth2.append( - np.min(pnt2line(P0=DEJ[0:(n_DEJ - 1), :], P1=DEJ[1:n_DEJ, :], p=point)))
    for point in points_in_lesion_dermis:
        depth2.append(np.min(pnt2line(P0=DEJ[0:(n_DEJ - 1), :], P1=DEJ[1:n_DEJ, :], p=point)))

    # bin all depths into depth_spacing px increments
    depths = np.array(depth1+depth2)
    depth_bins = np.arange(np.floor(np.min(depths)) // depth_spacing * depth_spacing,
                           np.ceil(np.max(depths) / depth_spacing) * depth_spacing + depth_spacing,
                           depth_spacing)

    # compute histograms
    DEJ_hist1 = np.histogram(np.array(depth1), bins=depth_bins)
    DEJ_hist2 = np.histogram(np.array(depth2), bins=depth_bins)
    DEJ_hist0 = np.histogram(depths, bins=depth_bins)
        

    return [skin_hist0, skin_hist1, skin_hist2, DEJ_hist0, DEJ_hist1, DEJ_hist2]



# list the directories under ../Chipcytometry_images_from8bit
chipcytometry_images_dir = '../Chipcytometry_images_from8bit/'
directories = [d for d in os.listdir(chipcytometry_images_dir) if os.path.isdir(os.path.join(chipcytometry_images_dir, d))]

res=[]
for folder in directories:
    roi_dir = chipcytometry_images_dir + folder + '/ZKWstitched_alignedtoCD45_cropped/ROI_regions/'
    if os.path.exists(roi_dir):
        print(folder)
        out = do_biopsy(roi_dir)

        # add results to the dataframe, convert depth from px to um by dividing by 2

        count, depth = out[0]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "biopsy", 'count': count[i], 'depth_skin': depth[i] / 2})

        count, depth = out[1]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "normal", 'count': count[i], 'depth_skin': depth[i] / 2})

        count, depth = out[2]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "lesion", 'count': count[i], 'depth_skin': depth[i] / 2})

        count, depth = out[3]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "biopsy", 'count': count[i], 'depth_epi': depth[i] / 2})

        count, depth = out[4]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "normal", 'count': count[i], 'depth_epi': depth[i] / 2})

        count, depth = out[5]
        for i in range(len(count)):
            res.append({'biopsy': folder, 'region': "lesion", 'count': count[i], 'depth_epi': depth[i] / 2})


# save the results to a csv file
df = pd.DataFrame(res)
df.to_csv(chipcytometry_images_dir+'/tissue_size_by_depth.csv', index=False)