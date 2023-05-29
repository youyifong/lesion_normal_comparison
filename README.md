# Lesion Mormal Comparison


## Methods

### Compute region area and length
A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis by the DEJ line. For each region, we expect three ROI files:
- region roi for epi. Type polygon, file name example AS_Les_epi_region_1_normal
- region roi for dermis. Type polygon, file name example AS_Les_dermis_region_1_normal
- length roi for the region. Type line or polyline, file name example AS_Les_length_region_1_normal
To compute the area of each region and the length of the region, run analysis.py in the ROI folder. 
```
 ml IPython/7.26.0-GCCcore-11.2.0
```


### Segmentation and spatial statistics

```
ml Python/3.9.6-GCCcore-11.2.0
ml cuDNN/8.2.2.26-CUDA-11.4.1
```

#### XX

XX is the label, e.g. CF, DD, AS etc.

```bash
export XX=AS
```

```bash
nohup python -m tsp runcellpose --f '$XX_Les_CD*.png' --l [3,0] --model cytotrain7 --s &
```

For CD123, increase flow to allow more masks in high density regions and set a min total intensity threshold to remove masks in the debris field in the lesion area

```bash
python -m tsp runcellpose --f 'DD_Les_CD123.png' --l [3,0] --model cytotrain7 --flow 0.5 --min_totalintensity 5000 --s
```

Measure distance to skin and DEJ
```bash
python -m tsp dist2boundary --f '$XX_Les_CD*_masks.csv' --boundaryroi [_Les_ROI_regions/DD_Les_skin_boundary.roi,DD_Les_ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 
```

Determine membership in regions
```
python -m tsp regionmembership --f 'DD_Les_CD*_masks_d2b.csv' --regionroi 'DD_Les_ROI_regions/DD_Les_*_region_*.roi'
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