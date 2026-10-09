# Analysis of 32-bit images

```{bash}
ml cuDNN/8.2.2.26-CUDA-11.4.1
ml Python/3.9.6-GCCcore-11.2.0
source ~/envs/tspenv/bin/activate
export CUBLAS_WORKSPACE_CONFIG=:4096:8 # for GPU computation reproducibility
cd /fh/fast/fong_y/images_from_K/lesion_normal_comparison
```

Analysis of images that are processed through a 32-bit import from the instrument. Older analysis files using 8-bit import are saved under the folder 8bit.



### Directory structure

The only clone of lesion_normal_comparison is on the cluster. 
- py. For image processing and cell segmentation. 
- R. For spatial analysis. Use Rstudio on laptop to write on the files from X drive (too slow to run from X drive). Use Rstudio online to write or run scripts. 

../Chipcytometry_images_from32bit contains the ZellKraftWerk tiff images of biopsy samples
- Each subfolder contains images from one subject, e.g., CF
  - Each subfolder contains one slide, representing one slice of a biopsy at one of the five timepoints of HSV infection and healing: Lesion, Healed, 2 weeks post heal (wph), 4wph, 8wph, and 8wph control. In addition, 8wph control is a contralateral biopsy of 8wph, so it is biopsied at the same time as 8wph.

  - For each slide folder, e.g. AS/AS1_Les, the content is organized as follows:

    - TSP_stitched
      - stitched full-size images of the sample 
      - aligned to CD45 if CD45 is available, if not, I just align CD3 to CD8.  They are not cropped, just stitched the active positions, or positions with images

      - ROI_regions folder. A biopsy is divided into contiguous regions. A region is either lesion (ulcerative) or normal (intact). Each region is further divided into epidermis and dermis. For each region, we expect three roi files under the ROI_regions folder:
         - skin boundary. Type: line or polyline. File name example: AS_Les_**skin_boundary**.roi
         - DEJ. Type: line or polyline. File name example: AS_Les_**DEJ**_v2.roi
         - whole tissue boundary roi. Type: polygon. File name example: AS_Les_**whole_tissue**.roi. Note that this file is not needed if the regions have tight boundaries.
         - epi region boundary. Type: polygon. File name example: AS_Les_**epi_region**_1_normal.roi
         - dermis region boundary. Type: polygon. File name example: AS_Les_**dermis_region**_1_normal.roi


### Background

Marker groups

- T cell: CD3 CD4 CD8 CD69 CD279 Foxp3
- Neutrophil: CD15
- Dendritic cells: CD11c CD123.  CD123 also stains blood vessels, so need human-in-loop to remove the segmented outline of blood vessels.
- Macrophage: CD68.  CD68 staining often result in bright epidermis and thus segmented.  Can also stain dendritic cells/dendrites in skin.  Require human-in-loop to remove these segmented signals.
- NK cells: CD56. This marker is expressed on both NK and nerve fibers, need to check if the staining on the latter affects the NK cell count.  Remove nerve fiber segmentation via human-in-loop.


### Tissue stats

```{bash}
# makes Chipcytometry_images_from32bit/tissue_size_by_depth_and_region.csv
nohup python py/get_biopsy_size_by_depth_and_region_32bit.py --grid 5 &


# count each intact/ulcer region within a biopsy separately
# makes Chipcytometry_images_from32bit/tissue_size_by_depth_and_region_2.csv
nohup python py/get_biopsy_size_by_depth_and_region_32bit_2.py --grid 5 &

```


### Segmentation

Performed by K. From K:
- There are problems with some CD68 segmentation because the antibody also stained the epidermis of some tissues. I tried to reduce the epidermis brightness with some success. Thus, some samples will have *CD68.png (original) and *CD68_2.png (my attempt to adjust epidermis).
- CD4T (CD3+CD4+) and CD8T (CD3+CD8+) cell populations are obtained by the following:
  - Segment CD3 with python -m tsp runcellpose
  - Shrink CD3 masks by 1 pixel via python -m tsp dilatemasks
    - Shrunken masks files are saved in the format *_m_c7_d-1.png
  - Use shrunken CD3 masks files in python -m tsp cellphenotyping to identify coexpressions in respective CD4 or CD8 images by fluorescence intensity.
    - Intensity_median analysis for (--m)
    - Set threshold to 9.8 for (--c)
    - 'True' value for (--p)
  - Cellphenotyping results include a csv file as one of the outputs.
    - The csv file is saved in the format *_masks.csv, and can be used in 'Extraction of spatial features'.


### Extraction of spatial features

Slide-specific xxx_region_area_length files are no longer needed because we now use tissue_size_by_depth_and_region.csv.

K performs the following steps to extract spatial features from the stitched images:
- Define an environment variable XX. In shell commands, ${XX} can be used to represent a generic biopsy, e.g., CF, DD, AS. Note that it needs to be enclosed within double quotes not single quotes for it to work.
- Enter ROI_regions folder
- Make a png file with region boundaries
- Enter the parent folder
- Cell segmentation
- Compute distance to skin and to DEJ for each mask by running python -m tsp dist2boundary
- Compute region membership by running python -m tsp regionmembership
for example, 

```{bash}
export YY=AS
export XX=AS1_Les
cd /fh/fast/fong_y/images_from_K/Chipcytometry_images_from32bit/${YY}/${XX}/TSP_stitched/ROI_regions
python -m tsp roifiles2mask --f "*region*.roi" --imagefile "../${XX}_CD3.png" --saveas "../${XX}_region_boundary.png" --saveoutlineonly

cd ..
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/${XX}_skin_boundary.roi,ROI_regions/${XX}_Les_DEJ.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region_*.roi"
```

When a slide has multiple discontiguous lesion or intact regions, the region membership is computed by passing all the region ROIs to the regionmembership command.


Human_In_Loop editing (HIL)

- Required when antibodies/markers like CD56, CD68, CD123 are also present in other structures (CD123 on blood vessels) or cells (CD56 on NK and neuronal cells)
- Also when tissue regions such as the epidermis is bright enough to get segmented
- This leads to over-counting the cell population of interest
- These unintended segmented masks should be removed before continuing with spatial analyses
  - For example, CD123
      - Segment CD123 via runcellpose  python -m tsp runcellpose --f 'XX_CD123.png' --model c7 --min_size 100 --flow 0.5 --l [0,0] --s --saveroi (required to generate the necessary .txt file)
      - Convert XX_CD123_cp_outlines.txt to ImageJ ROI's using 'imagej_roi_converter.py'
          - Drag and drop 'image_roi.converter.py' into ImageJ to open a conversion window
          - Open CD123 image
          - Click 'Run' in the window to find and open XX_CD123_cp_outlines.txt to generate ImageJ ROI's
      - Look through the image and ROI's to delete unintended ROI's (HIL)
      - When done, save all ROI's as a .zip file
      - Convert .zip back to mask files using 'roifiles2masks' command to obtain .csv file for the HIL masks
      - Run 'dist2boundary' with the .csv file, then 'regionmembership' 

  
### Spatial analysis

Run Rstudio online on Gizmo and source _run.R.

The density results saved in output/csv/ are in the unit of count/mm2, not count/100umx100um as in the density depth plots. 

Report pdf is rendered by a line of R code in _run.R.
