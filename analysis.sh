#### DD

# cell segmentation
nohup python -m tsp runcellpose --f 'DD_Les_CD*.png' --l [3,0] --model cytotrain7 --s &

# measure distance to skin and DEJ
python -m tsp dist2boundary --f 'DD_Les_CD*_masks.csv' --boundaryroi [DD_Les_ROI_regions/DD_Les_skin_boundary.roi,DD_Les_ROI_regions/DD_Les_DEJviaDAPI_v2.roi] 

# determine membership in regions
python -m tsp regionmembership --f 'DD_Les_CD*_masks_d2b.csv' --regionroi 'DD_Les_ROI_regions/DD_Les_*_region_*.roi'




#### CF

# cell segmentation
nohup python -m tsp runcellpose --f 'CF_Les_CD*.png' --l [3,0] --model cytotrain7 --s &

# measure distance to skin and DEJ
python -m tsp dist2boundary --f 'CF_Les_CD*_masks.csv' --boundaryroi [CF_Les_ROI_regions/CF_Les_skin_boundary.roi,CF_Les_ROI_regions/CF_Les_DEJ_v2.roi] 

# determine membership in regions
python -m tsp regionmembership --f 'CF_Les_CD*_masks_d2b.csv' --regionroi 'CF_Les_ROI_regions/CF_Les_*_region_*.roi'



