"""
Compute the size of biopsy in the unit of 100um x 100um by depth in 200 px (equivalently, 100 um) increments for all biopsies.
Save the results in Chipcytometry_images_from32bit/biopsy_size_by_depth_and_region.csv. The depth column is the starting point of the bin.
It takes a long time to run at the default grid size of 5 px, so you can change it to 10 px or 20 px to speed up the process.
"""

from shapely.geometry import Polygon, Point
import glob, os, re
from read_roi import read_roi_file
import pandas as pd
import numpy as np
from tsp.geom import pnt2line

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--grid', type=float, default=5.0, help='size of patches in px')
args = parser.parse_known_args()[0]
grid_spacing = args.grid

depth_spacing = 200.0 # px


def extract_region_id(filepath):
    """Extract the region identifier from a ROI filename, e.g. '1' from '*_region_1_normal.roi'."""
    basename = os.path.basename(filepath)
    m = re.search(r'_region_([^_]+)_(?:normal|lesion)', basename)
    if m:
        return m.group(1)
    return os.path.splitext(basename)[0]


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
        raise RuntimeError("There should be one and only one DEJ boundary roi file.")
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
    print(normal_dermis_regions_roi)
    if len(normal_dermis_regions_roi)==0:
        raise RuntimeError("There are no normal dermis regions.")
    normal_dermis_regions = {}
    for x in normal_dermis_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        normal_dermis_regions[extract_region_id(x)] = Polygon(list(zip(val["x"], val["y"])))

    normal_epi_regions_roi = glob.glob(roi_dir+'*_epi_region_*normal*.roi')
    print(normal_epi_regions_roi)
    normal_epi_regions = {}
    for x in normal_epi_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        normal_epi_regions[extract_region_id(x)] = Polygon(list(zip(val["x"], val["y"])))

    lesion_dermis_regions_roi = glob.glob(roi_dir+'*_dermis_region_*lesion*.roi')
    print(lesion_dermis_regions_roi)
    lesion_dermis_regions = {}
    for x in lesion_dermis_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        lesion_dermis_regions[extract_region_id(x)] = Polygon(list(zip(val["x"], val["y"])))

    lesion_epi_regions_roi = glob.glob(roi_dir+'*_epi_region_*lesion*.roi')
    print(lesion_epi_regions_roi)
    lesion_epi_regions = {}
    for x in lesion_epi_regions_roi:
        val = list(read_roi_file(x).values()) [0]
        lesion_epi_regions[extract_region_id(x)] = Polygon(list(zip(val["x"], val["y"])))


    """
    Create a grid of points and classify each point as normal or lesion based on the region it is in.
    Per-region point lists track epi and dermis separately (needed for signed DEJ depth).
    """

    all_normal_region_ids = sorted(set(normal_dermis_regions) | set(normal_epi_regions))
    all_lesion_region_ids = sorted(set(lesion_dermis_regions) | set(lesion_epi_regions))

    normal_points_by_region = {rid: {'epi': [], 'dermis': []} for rid in all_normal_region_ids}
    lesion_points_by_region = {rid: {'epi': [], 'dermis': []} for rid in all_lesion_region_ids}

    # Get bounds of the biopsy
    min_x, min_y, max_x, max_y = biopsy.bounds

    # Generate grid points
    x_coords = np.arange(min_x, max_x, grid_spacing)
    y_coords = np.arange(min_y, max_y, grid_spacing)
    grid_points = [(x, y) for x in x_coords for y in y_coords]

    for p in grid_points:
        classified = False

        for rid, polygon in normal_epi_regions.items():
            if polygon.contains(Point(p)):
                normal_points_by_region[rid]['epi'].append(p)
                classified = True
                break
        if classified:
            continue

        for rid, polygon in normal_dermis_regions.items():
            if polygon.contains(Point(p)):
                normal_points_by_region[rid]['dermis'].append(p)
                classified = True
                break
        if classified:
            continue

        for rid, polygon in lesion_dermis_regions.items():
            if polygon.contains(Point(p)):
                lesion_points_by_region[rid]['dermis'].append(p)
                classified = True
                break
        if classified:
            continue

        for rid, polygon in lesion_epi_regions.items():
            if polygon.contains(Point(p)):
                lesion_points_by_region[rid]['epi'].append(p)
                break


    """
    Count the number of grid points by distance to skin in depth_spacing px increments.
    """

    def skin_depth(point):
        return np.min(pnt2line(P0=skin[0:(n_skin - 1), :], P1=skin[1:n_skin, :], p=point))

    normal_skin_depths = {}
    for rid, pts in normal_points_by_region.items():
        normal_skin_depths[rid] = [skin_depth(p) for p in pts['epi'] + pts['dermis']]

    lesion_skin_depths = {}
    for rid, pts in lesion_points_by_region.items():
        lesion_skin_depths[rid] = [skin_depth(p) for p in pts['epi'] + pts['dermis']]

    all_skin_depths = ([d for depths in normal_skin_depths.values() for d in depths] +
                       [d for depths in lesion_skin_depths.values() for d in depths])

    skin_bins = np.arange(0,
                          np.ceil(np.max(all_skin_depths) / depth_spacing) * depth_spacing + depth_spacing,
                          depth_spacing)

    skin_hist_biopsy = np.histogram(np.array(all_skin_depths), bins=skin_bins)
    skin_hists_normal = {rid: np.histogram(np.array(depths), bins=skin_bins)
                         for rid, depths in normal_skin_depths.items()}
    skin_hists_lesion = {rid: np.histogram(np.array(depths), bins=skin_bins)
                         for rid, depths in lesion_skin_depths.items()}


    """
    Count the number of grid points by distance to DEJ in depth_spacing px increments.
    Epidermis points count as negative depth.
    """

    def DEJ_depth(point, is_epi):
        d = np.min(pnt2line(P0=DEJ[0:(n_DEJ - 1), :], P1=DEJ[1:n_DEJ, :], p=point))
        return -d if is_epi else d

    normal_DEJ_depths = {}
    for rid, pts in normal_points_by_region.items():
        normal_DEJ_depths[rid] = ([DEJ_depth(p, True) for p in pts['epi']] +
                                   [DEJ_depth(p, False) for p in pts['dermis']])

    lesion_DEJ_depths = {}
    for rid, pts in lesion_points_by_region.items():
        lesion_DEJ_depths[rid] = ([DEJ_depth(p, True) for p in pts['epi']] +
                                   [DEJ_depth(p, False) for p in pts['dermis']])

    all_DEJ_depths = ([d for depths in normal_DEJ_depths.values() for d in depths] +
                      [d for depths in lesion_DEJ_depths.values() for d in depths])

    DEJ_bins = np.arange(np.floor(np.min(all_DEJ_depths)) // depth_spacing * depth_spacing,
                         np.ceil(np.max(all_DEJ_depths) / depth_spacing) * depth_spacing + depth_spacing,
                         depth_spacing)

    DEJ_hist_biopsy = np.histogram(np.array(all_DEJ_depths), bins=DEJ_bins)
    DEJ_hists_normal = {rid: np.histogram(np.array(depths), bins=DEJ_bins)
                        for rid, depths in normal_DEJ_depths.items()}
    DEJ_hists_lesion = {rid: np.histogram(np.array(depths), bins=DEJ_bins)
                        for rid, depths in lesion_DEJ_depths.items()}

    return {
        'skin_biopsy': skin_hist_biopsy,
        'skin_normal': skin_hists_normal,
        'skin_lesion': skin_hists_lesion,
        'DEJ_biopsy': DEJ_hist_biopsy,
        'DEJ_normal': DEJ_hists_normal,
        'DEJ_lesion': DEJ_hists_lesion,
    }



chipcytometry_images_dir = '../Chipcytometry_images_from32bit/'
children_dirs = [d for d in os.listdir(chipcytometry_images_dir) if os.path.isdir(os.path.join(chipcytometry_images_dir, d))]
directories = []
for child in children_dirs:
    child_path = os.path.join(chipcytometry_images_dir, child)
    directories.extend(
        [os.path.join(child, g) for g in os.listdir(child_path) if os.path.isdir(os.path.join(child_path, g)) and g.endswith('_Les')]
    )

res=[]

for folder in directories:
    roi_dir = chipcytometry_images_dir + folder + '/TSP_stitched/ROI_regions/'
    if os.path.exists(roi_dir):
        print(roi_dir)

        try:
            out = do_biopsy(roi_dir)
        except RuntimeError as e:
            print(f"Error processing {folder}: {e}")
            continue

        # convert count from number of 5x5 px patches to number of 100x100 um by dividing by 40**2 because each px is 1/2 um
        # convert depth from px to um by dividing by 2

        count, depth = out['skin_biopsy']
        for i in range(len(count)):
            res.append({'biopsy': folder,
                        'region': "biopsy",
                        'region_id': None,
                        'dist_type': 'dist2skin',
                        'area': count[i] / (depth_spacing / grid_spacing)**2,
                        'depth': depth[i] / 2})

        for rid, (count, depth) in out['skin_normal'].items():
            for i in range(len(count)):
                res.append({'biopsy': folder,
                            'region': "intact",
                            'region_id': rid,
                            'dist_type': 'dist2skin',
                            'area': count[i] / (depth_spacing / grid_spacing)**2,
                            'depth': depth[i] / 2})

        for rid, (count, depth) in out['skin_lesion'].items():
            for i in range(len(count)):
                res.append({'biopsy': folder,
                            'region': "ulcerative",
                            'region_id': rid,
                            'dist_type': 'dist2skin',
                            'area': count[i] / (depth_spacing / grid_spacing)**2,
                            'depth': depth[i] / 2})

        count, depth = out['DEJ_biopsy']
        for i in range(len(count)):
            res.append({'biopsy': folder,
                        'region': "biopsy",
                        'region_id': None,
                        'dist_type': 'dist2DEJ',
                        'area': count[i] / (depth_spacing / grid_spacing)**2,
                        'depth': depth[i] / 2})

        for rid, (count, depth) in out['DEJ_normal'].items():
            for i in range(len(count)):
                res.append({'biopsy': folder,
                            'region': "intact",
                            'region_id': rid,
                            'dist_type': 'dist2DEJ',
                            'area': count[i] / (depth_spacing / grid_spacing)**2,
                            'depth': depth[i] / 2})

        for rid, (count, depth) in out['DEJ_lesion'].items():
            for i in range(len(count)):
                res.append({'biopsy': folder,
                            'region': "ulcerative",
                            'region_id': rid,
                            'dist_type': 'dist2DEJ',
                            'area': count[i] / (depth_spacing / grid_spacing)**2,
                            'depth': depth[i] / 2})


df = pd.DataFrame(res)
df.to_csv(chipcytometry_images_dir+'/tissue_size_by_depth_and_region_2.csv', index=False)