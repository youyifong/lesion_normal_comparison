# Lesion Mormal Comparison


## Methods

### Compute region area and length

Environment requirement:
```
 ml IPython/7.26.0-GCCcore-11.2.0
```


A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis. For each region, we expect three roi files:
- epi region boundary. Type: polygon. File name example: AS_Les_**epi_region**_1_normal.roi
- dermis region boundary. Type: polygon. File name example: AS_Les_**dermis_region**_1_normal.roi
- region length. Type: line or polyline. File name example: AS_Les_**length**_1_normal.roi. Note that **region** is not part of the file name, which allows us to pick up the region boundary files by matching with 'region'.

In addition to the three roi files per region, each biopsy ROI folder also contains:
- whole tissue boundary roi. Type: polygon. File name example: AS_Les_**whole_tissue**.roi
- skin boundary. Type: line or polyline. File name example: AS_Les_**skin_boundary**.roi
- DEJ. Type: line or polyline. File name example: AS_Les_**DEJ**_v2.roi

With these ROI files ready, enter the directory where the ROI files are and run analysis.py in IPython. The script writes region_area_length.csv in the parent directory.


### Segmentation and spatial statistics

Environment requirement:
```
ml Python/3.9.6-GCCcore-11.2.0
ml cuDNN/8.2.2.26-CUDA-11.4.1
source envs/tspenv/bin/activate
```

#### XX

XX is the label, e.g. CF, DD, AS etc. In the shell command, ${XX} can be used to denote the biopsy if used inside double quotes and not single quotes.

```bash
export XX=AS
```

```bash
nohup python -m tsp runcellpose --f "*_Les_CD*.png" --l [0,0] --model cytotrain7 --s &
```

For CD123 from CF and DD, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
```bash
python -m tsp runcellpose --f "DD_Les_CD123.png" --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
```

Measure distance to skin and DEJ
```bash
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [${XX}_Les_ROI_regions/${XX}_Les_skin_boundary.roi,${XX}_Les_ROI_regions/*DEJ*.roi] 
```

Determine membership in regions
```
python -m tsp regionmembership --f "${XX}_Les_CD*_masks_d2b.csv" --regionroi "${XX}_Les_ROI_regions/${XX}_Les_*_region_*.roi"
```


#### DD

```bash
nohup python -m tsp runcellpose --f 'DD_Les_CD*.png' --l [3,0] --model cytotrain7 --s &
```

For CD123, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area

```bash
python -m tsp runcellpose --f 'DD_Les_CD123.png' --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
```

Measure distance to skin and DEJ
```bash
python -m tsp dist2boundary --f 'DD_Les_CD*_masks.csv' --boundaryroi [DD_Les_ROI_regions/DD_Les_skin_boundary.roi,DD_Les_ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 
```

Determine membership in regions
```
python -m tsp regionmembership --f 'DD_Les_CD*_masks_d2b.csv' --regionroi 'DD_Les_ROI_regions/DD_Les_*_region_*.roi'
python -m tsp regionmembership --f "${XX}_Les_CD3_masks_d2b.csv" --regionroi "${XX}_Les_ROI_regions/*_region_*.roi"
```


#### CF

```bash
nohup python -m tsp runcellpose --f 'CF_Les_CD*.png' --l [3,0] --model cytotrain7 --s &
```

For CD123, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area
```bash
python -m tsp runcellpose --f 'CF_Les_CD123.png' --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
```
    
For CD279, one cell is given two masks
```bash
python -m tsp roifiles2mask --f CF_Les_CD279.roi --height 5296 --width 4040
python -m tsp runcellpose --f 'CF_Les_CD279.png' --l [3,0] --model cytotrain7 --cellprob -5
```

Measure distance to skin and DEJ
```bash
python -m tsp dist2boundary --f 'CF_Les_CD*_masks.csv' --boundaryroi [CF_Les_ROI_regions/CF_Les_skin_boundary.roi,CF_Les_ROI_regions/CF_Les_DEJ_v2.roi] 
```

determine membership in regions
```bash
python -m tsp regionmembership --f 'CF_Les_CD*_masks_d2b.csv' --regionroi 'CF_Les_ROI_regions/CF_Les_*_region_*.roi'
```


## Results

[Report](report.pdf).