# Workbook for Lesion Normal Comparison


Environment requirement:
```
ml Python/3.9.6-GCCcore-11.2.0
ml cuDNN/8.2.2.26-CUDA-11.4.1
source envs/tspenv/bin/activate
```

Define an environment variable XX. In shell commands, ${XX} can be used to represent a generic biopsy, e.g., CF, DD, AS. Note that it needs to be enclosed within double quotes not single quotes for it to work.

```bash
export XX=JM
```


## Compute region area and length


A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis. For each region, we expect three roi files:
- epi region boundary. Type: polygon. File name example: AS_Les_**epi_region**_1_normal.roi
- dermis region boundary. Type: polygon. File name example: AS_Les_**dermis_region**_1_normal.roi
- region length. Type: line or polyline. File name example: AS_Les_**length**_1_normal.roi. Note that **region** is not part of the file name, which allows us to pick up all region roi (both epi and dermis) files by matching with 'region'.

In addition to the three roi files per region, each biopsy ROI folder also contains:
- skin boundary. Type: line or polyline. File name example: AS_Les_**skin_boundary**.roi
- DEJ. Type: line or polyline. File name example: AS_Les_**DEJ**_v2.roi
- whole tissue boundary roi. Type: polygon. File name example: AS_Les_**whole_tissue**.roi. Note that this file is not needed if the regions have tight boundaries.

All the ROI files are under a folder named ROI_regions.


```bash
cd /fh/fast/fong_y/images_from_K/ROI_lesion_vs_normal_TSP
cd ${XX}_Les_stitched_alignedtoCD45_cropped/ROI_regions
```

Make a png file with region boundaries. 
```bash
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_Les_CD3.png" --savemasksonly --saveas "${XX}_regions" 
```

Create region_area_length.csv.
```bash
ipython # run ROI_work.py
```


## Cell segmentation 


The default command is:
```bash
python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s --normalize100  --min_size 100 --flow 0.5 
```
- normalize100 is set because 1) it avoids errors when the number of cells in the image is small, 2) it leads to more specific results, e.g. some fibers will not be predicted masks. The downside is that it can be less sensitive to dim signals.


Special care are needed for the following:

For CD123 CF/DD, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
```bash
python -m tsp runcellpose --f "CF_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
python -m tsp runcellpose --f "DD_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
```

For CD279 CF, one cell is wrongly segmented as two masks under the default setting and 
```bash
python -m tsp runcellpose --f 'CF_Les_CD279.png' --l [3,0] --model cytotrain7 --cellprob -5
```

For CD68 JM, a high background makes it necessary to set a higher min_avgintensity
```bash
python -m tsp runcellpose --f 'JM_Les_CD68.png' --l [3,0] --model cytotrain7 --normalize100 --min_size 100 --min_avgintensity 30 --flow 0.5
```


## Compute spatial statistics

Measure distance to skin and DEJ
```bash
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [${XX}_Les_ROI_regions/${XX}_Les_skin_boundary.roi,${XX}_Les_ROI_regions/*DEJ*.roi] 
```

Determine membership in regions
```
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "${XX}_Les_ROI_regions/${XX}_Les_*_region_*.roi"
```
