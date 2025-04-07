# Analysis of 32-bit images

Analysis of images that are processed through a 32-bit import from the instrument. Older analysis files using 8-bit import are saved under the folder 8bit.

### Setup
Load modules on Volta or the Hutch cluster.
```{bash}
ml cuDNN/8.2.2.26-CUDA-11.4.1
ml Python/3.9.6-GCCcore-11.2.0
ml IPython/7.26.0-GCCcore-11.2.0
source ~/envs/tspenv/bin/activate
cd /fh/fast/fong_y/images_from_K/lesion_normal_comparison
export CUBLAS_WORKSPACE_CONFIG=:4096:8 # for GPU computation reproducibility
```


### Directory structure

- ../Chipcytometry_images_from32bit contains the ZellKraftWerk tiff images of biopsy samples
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

  
### Tissue stats

```{bash}
python code/get_biopsy_size_by_depth_and_region_32bit.py --grid 5
```


### Segmentation

Performed by K. From K:
- There are problems with some CD68 segmentation because the antibody also stained the epidermis of some tissues. I tried to reduce the epidermis brightness with some success. Thus, some samples will have *CD68.png (original) and *CD68_2.png (my attempt to adjust epidermis).



### Extraction of spatial features

This is done by K. Note that xxx_region_area_length is no longer needed because we now use tissue_size_by_depth_and_region.csv.

Steps:
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
python -m tsp dist2boundary --f "*_Les_CD*_masks.csv" --boundaryroi [ROI_regions/${XX}_skin_boundary.roi,ROI_regions/DD_Les_DEJviaDAPI-WGA_alternate.roi] 
python -m tsp regionmembership --f "${XX}_CD*_masks_d2b.csv" --regionroi "ROI_regions/${XX}_*_region_*.roi"
```

