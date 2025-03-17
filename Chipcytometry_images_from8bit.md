# Chipcytometry_images_from8bit

This folder contains images that are processed through a 8-bit import from the instrument.

The Chipcytometry_images_from8bit folder should be on the same level as the clone of lesion_normal_comparison repo.

- tissue_size_by_depth_from_skin_for_all_biopsies.csv. This file is created by 
```{bash}
python code/get_biopsy_size_by_depth_and_region_8bit.py
```

- Each folder under Chipcytometry_images_from8bit contains ZellKraftWerk tiff images of biopsy samples at one of the five timepoints of HSV infection and healing: Lesion, Healed, 2 weeks post heal (wph), 4wph, 8wph, and 8wph control. In addition, 8wph control is a contralateral biopsy of 8wph, so it is biopsied at the same time as 8wph.

- Under each folder, e.g. CF_Les, the content is organized as follows:

  - Individual_positions_preprocessed
    - AS15065_Les_BG_4x4.png
    - CD1a, CD3, CD4, etc. folders
      - Images from individual positions (position1_CD1a, position2_CD1a, etc.) of that marker. These positions can be arranged according to their respective reference images to create a stitched image, or full-size image, of each marker for each sample.

  - ZKWstitched_alignedtoCD45_cropped
    - stitched full-size images of the sample
    - ROI_regions folder
       - skin boundary. Type: line or polyline. File name example: AS_Les_**skin_boundary**.roi
       - DEJ. Type: line or polyline. File name example: AS_Les_**DEJ**_v2.roi
       - whole tissue boundary roi. Type: polygon. File name example: AS_Les_**whole_tissue**.roi. Note that this file is not needed if the regions have tight boundaries.



## Methods

Stitching command:
>python -m tsp alignimages --ref_image ../Chipcytometry_images_from8bit/CF_Les/2_individual_positions_PE/CD45/M926910_Position1_CD45  --image2  ../Chipcytometry_images_from8bit/CF_Les/1_individual_positions_AF488/CD8/M926910_Position1_CD8a-AF488 
- Stitching using the stitching class https://github.com/OpenStitching/stitching

