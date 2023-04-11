# Distance to skin boundary

## Methods
For each sample, we will compare the distribution of each markers in the normal and lesion regions.

- Markers for DD:       CD3, CD4, CD8a, CD11c, CD15, CD45, CD68, CD69, CD123, CD279 

- Markers for CF: CD1a, CD3, CD4, CD8a, CD11c, CD15, CD45,       CD69, CD123, CD279 
(CD68 did not pass QC)

Todo

- compute areas in normal/lesion, epidermis/dermis 
- normalize cell count by areas
- save cell counts as tables
- refine cell segmentation results


## Results

<object data="R/tables/DD_counts.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="R/tables/DD_counts.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="R/tables/DD_counts.pdf">Download PDF</a>.</p>
    </embed>
</object>

### DEJ DD
![plot](R/figures/DD_dist2DEJ_compare_markers.png?raw=true)

### Skin DD  
![plot](R/figures/DD_dist2skin_compare_markers.png?raw=true)


### DEJ CF 
![plot](R/figures/CF_dist2DEJ_compare_markers.png?raw=true)

### Skin CF 
![plot](R/figures/CF_dist2skin_compare_markers.png?raw=true)


