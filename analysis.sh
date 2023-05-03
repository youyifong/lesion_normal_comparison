#### DD

# cell segmentation
nohup python -m tsp runcellpose --f 'DD_Les_CD*.png' --l [3,0] --model cytotrain7 --s &

    # for CD123 DD, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
    python -m tsp runcellpose --f 'DD_Les_CD123.png' --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s

# measure distance to skin and DEJ
python -m tsp dist2boundary --f 'DD_Les_CD*_masks.csv' --boundaryroi [DD_Les_ROI_regions/DD_Les_skin_boundary.roi,DD_Les_ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 

# determine membership in regions
python -m tsp regionmembership --f 'DD_Les_CD*_masks_d2b.csv' --regionroi 'DD_Les_ROI_regions/DD_Les_*_region_*.roi'




#### CF

# cell segmentation
nohup python -m tsp runcellpose --f 'CF_Les_CD*.png' --l [3,0] --model cytotrain7 --s &

    # CD123 is different in the lesion
    python -m tsp runcellpose --f 'CF_Les_CD123.png' --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
    
    # CD279. one cell is given two masks
    python -m tsp roifiles2mask --f CF_Les_CD279.roi --height 5296 --width 4040
    python -m tsp runcellpose --f 'CF_Les_CD279.png' --l [3,0] --model cytotrain7 --cellprob -5

# measure distance to skin and DEJ
python -m tsp dist2boundary --f 'CF_Les_CD*_masks.csv' --boundaryroi [CF_Les_ROI_regions/CF_Les_skin_boundary.roi,CF_Les_ROI_regions/CF_Les_DEJ_v2.roi] 

# determine membership in regions
python -m tsp regionmembership --f 'CF_Les_CD*_masks_d2b.csv' --regionroi 'CF_Les_ROI_regions/CF_Les_*_region_*.roi'



