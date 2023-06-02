# Workbook for Lesion Normal Comparison

## Methods
A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis. For each region, we expect three roi files under the ROI_regions folder:
- epi region boundary. Type: polygon. File name example: AS_Les_**epi_region**_1_normal.roi
- dermis region boundary. Type: polygon. File name example: AS_Les_**dermis_region**_1_normal.roi
- region length. Type: line or polyline. File name example: AS_Les_**length**_1_normal.roi. Note that **region** is not part of the file name, which allows us to pick up all region roi (both epi and dermis) files by matching with 'region'.

In addition to the three roi files per region, each ROI_regions folder also contains:
- skin boundary. Type: line or polyline. File name example: AS_Les_**skin_boundary**.roi
- DEJ. Type: line or polyline. File name example: AS_Les_**DEJ**_v2.roi
- whole tissue boundary roi. Type: polygon. File name example: AS_Les_**whole_tissue**.roi. Note that this file is not needed if the regions have tight boundaries.

The default command for cell segmentation is:
>python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
- use [3,0] for RGB images (CF and DD) and [0,0] for grayscale images (all others)
- normalize100 is set because 1) it avoids errors when the number of cells in the image is small, 2) it leads to more specific results, e.g. some fibers will not be predicted masks. The downside is that it can be less sensitive to dim signals
- need to remove *masks* files before running
- \*_Les_CD*.png is interpreted by glob, which does not support regular expressions and only supports Unix shell-style wildcards


## Commands

Since special care are needed for some markers and file names do not always follow the same naming convention, we list the commands for each tissue separately.


Steps:
- Define an environment variable XX. In shell commands, ${XX} can be used to represent a generic biopsy, e.g., CF, DD, AS. Note that it needs to be enclosed within double quotes not single quotes for it to work.
- Enter ROI_regions folder
- Make a png file with region boundaries
- Create region_area_length.csv
- Enter the parent folder
- Cell segmentation
- Compute distance to skin and to DEJ for each mask
- Compute region membership

Special handling for segmentation:
- DD CD123: increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
- CF CD123: increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
- CF CD279, one cell is wrongly segmented as two masks under the default setting and 
- JM CD68: a high background makes it necessary to set a higher min_avgintensity


```bash
# load the applications:
ml Python/3.9.6-GCCcore-11.2.0
ml cuDNN/8.2.2.26-CUDA-11.4.1
source envs/tspenv/bin/activate

# determines wehther to run segmentation
do_segmentation = 1

export XX=DD
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/DD_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f "DD_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/${XX}_Les_skin_boundary.roi,ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=CF
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/CF_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f "CF_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
    python -m tsp runcellpose --f 'CF_Les_CD279.png' --l [3,0] --model cytotrain7 --cellprob -5
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/CF_Les_skin_boundary.roi,ROI_regions/CF_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=JM
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/JM_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f 'JM_Les_CD68.png' --l [0,0] --model cytotrain7 --normalize100 --min_size 100 --min_avgintensity 30 --flow 0.5
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/JM_Les_skin_boundary.roi,ROI_regions/JM_Les_DEJ_v3.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=CJ
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/${XX}_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/AS_Les_skin_boundary_v2.roi,ROI_regions/AS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=AS
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/${XX}_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/CJ_Les_skin_boundary_v3.roi,ROI_regions/CJ_Les_DEJ_v3.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=LG
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/${XX}_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LG_Les_skin_boundary.roi,ROI_regions/LG_Les_DEJ.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"

export XX=LS
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP/${XX}_Les_stitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "../${XX}_Les_regions" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/ROI_work.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LS_Les_skin_boundary.roi,ROI_regions/LS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_Les_*_region_*.roi"
```