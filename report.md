# Distance to skin boundary

## Methods
For each sample, we will compare the distribution of each markers in the normal and lesion regions.

- Samples: DD, CF
- Markers: CD3, CD4, CD8a, CD15, CD45, CD68, CD69, CD123, CD279, and others (to be added)


Considerations on region boundary and distance metric:
1. Should we draw the boundary between normal and lesion based on WGA only or all markers? 
2. Shall each ROI be convex?

For (1), the two choices lead are answering different questions. Drawing boundary based on WGA only answers a stronger question than drawing the boundary based on all markers. 

For (2), a convex shape is a shape that a line connecting any two points in the shape completely lies within the shape. Practically speaking, in this context it is equivalent to drawing a straight line to separate normal from lesion. Pros for straight boundary: the span between any cell and the closest point on the skin lies completely within the region. So the distance measure is well defined. In contrast, if a line from a cell to its closest point on the skin surface cuts through both normal and lesion regions, the interpretation of distance to skin surface is muddled. Pros for curved boundary: Regions (which is a type of neighborhoods) do not necessarily have geometric shapes. 

These two questions are also intertwined. If we use WGA only to draw boundaries, we can only use straight line boundaries because the difference between normal and lesion is only visible in epidermis in WGA staining. A straight line can be drawn in epidermis and extend to dermis, but not a curve.

## Results

![plot](R/figures/CF_dist2skin_compare_regions_CD1a.png?raw=true)
   
![plot](R/figures/DD_dist2skin_compare_regions_CD3.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD3.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD4.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD4.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD8a.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD8a.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD11c.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD15.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD15.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD45.png?raw=true)
  
![plot](R/figures/CF_dist2skin_compare_regions_CD45.png?raw=true)
  
![plot](R/figures/DD_dist2skin_compare_regions_CD68.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD69.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD123.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD123.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_regions_CD279.png?raw=true)
    
![plot](R/figures/CF_dist2skin_compare_regions_CD279.png?raw=true)
    
![plot](R/figures/DD_dist2skin_compare_CD4_CD15.png?raw=true)
  