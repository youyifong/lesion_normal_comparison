# Lesion vs Normal Comparison

## Setup
Load modules on Volta or the Hutch cluster.
```{bash}
ml cuDNN/8.2.2.26-CUDA-11.4.1
ml Python/3.9.6-GCCcore-11.2.0
ml IPython/7.26.0-GCCcore-11.2.0
```

Load python virtual environment. 
```{bash}
source ~/envs/tspenv/bin/activate
cd /fh/fast/fong_y/images_from_K/lesion_normal_comparison
export CUBLAS_WORKSPACE_CONFIG=:4096:8 # for GPU computation reproducibility
```


## Data

Images are described in Chipcytometry_images_from32bit.md

A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis. For each region, we expect three roi files under the ROI_regions folder:
- epi region boundary. Type: polygon. File name example: AS_Les_**epi_region**_1_normal.roi
- dermis region boundary. Type: polygon. File name example: AS_Les_**dermis_region**_1_normal.roi
- region length. Type: line or polyline. File name example: AS_Les_**length**_1_normal.roi. Note that **region** is not part of the file name, which allows us to pick up all region roi (both epi and dermis) files by matching with 'region'.



## Methods 

Generate tissue size data.
```{bash}
python code/get_biopsy_size_by_depth_and_region_32bit.py --grid 5
```


The default command for cell segmentation is:
>python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
- use [3,0] for RGB images (CF and DD) and [0,0] for grayscale images (all others)
- normalize100 is set because 1) it avoids errors when the number of cells in the image is small, 2) it leads to more specific results, e.g. some fibers will not be predicted masks. The downside is that it can be less sensitive to dim signals
- need to remove *masks* files before running
- \*_Les_CD*.png is interpreted by glob, which does not support regular expressions and only supports Unix shell-style wildcards


Since special care are needed for some markers and file names do not always follow the same naming convention, we list the commands for each tissue separately.

Steps:
- Define an environment variable XX. In shell commands, ${XX} can be used to represent a generic biopsy, e.g., CF, DD, AS. Note that it needs to be enclosed within double quotes not single quotes for it to work.
- Enter ROI_regions folder
- Make a png file with region boundaries
- Create region_area_length.csv by running get_region_stats.py.
- Enter the parent folder
- Cell segmentation
- Compute distance to skin and to DEJ for each mask by running python -m tsp dist2boundary
- Compute region membership by running python -m tsp regionmembership

Special handling for segmentation:
- DD CD123: increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
- CF CD123: increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
- CF CD279, one cell is wrongly segmented as two masks under the default setting and 
- JM CD68: a high background makes it necessary to set a higher min_avgintensity

R scripts can be run from Rstudio for the Hutch cluster (need VPN).

```bash
# determines whether to run segmentation
do_segmentation=1


export XX=AS/AS1_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/TSP_stitched

cd ROI_regions
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 

cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f "DD_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/${XX}_skin_boundary.roi,ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region_*.roi"


export XX=CF_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [3,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f "CF_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
    python -m tsp runcellpose --f 'CF_Les_CD279.png' --l [3,0] --model cytotrain7 --cellprob -5
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/CF_Les_skin_boundary.roi,ROI_regions/CF_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=JM_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
    python -m tsp runcellpose --f 'JM_Les_CD68.png' --l [0,0] --model cytotrain7 --normalize100 --min_size 100 --min_avgintensity 30 --flow 0.5
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/JM_Les_skin_boundary.roi,ROI_regions/JM_Les_DEJ_v3.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=CJ_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/AS_Les_skin_boundary_v2.roi,ROI_regions/AS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=AS_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/CJ_Les_skin_boundary_v3.roi,ROI_regions/CJ_Les_DEJ_v3.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=LG_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LG_Les_skin_boundary.roi,ROI_regions/LG_Les_DEJ.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=LS_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${XX}/ZKWstitched_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LS_Les_skin_boundary.roi,ROI_regions/LS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=CF_1Heal
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/CF_post_heal/${XX}_alignedtoCD45_cropped
cd ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LS_Les_skin_boundary.roi,ROI_regions/LS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"


export XX=CF_2w
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/CF_post_heal/${XX}_alignedtoCD279_cropped
cd ${XX}_ROI
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --savemasksonly --saveas "../${XX}_region_boundary" 
python /fh/fast/fong_y/images_from_K/lesion_normal_comparison/code/get_region_stats.py
cd ..
if [ $do_segmentation -eq 1 ]; then
    python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
fi
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/LS_Les_skin_boundary.roi,ROI_regions/LS_Les_DEJ_v2.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region*.roi"

```



## Report

[Report](report/report.pdf).

