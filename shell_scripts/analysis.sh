nohup python -m tsp runcellpose --f 'DD_Les_CD*.png' --l [3,0] --model cytotrain7 &

python -m tsp dist2boundary --f 'DD_Les_CD*_masks.csv' --boundaryroi DD_Les_skin_boundary.roi

python -m tsp regionmembership --f 'DD_Les_CD*_masks_d2b.csv' --regionroi [DD_Les_normal_side.roi,DD_Les_lesion_side.roi]


python -m tsp dist2boundary --f 'CF_Les_CD*_masks.csv' --boundaryroi CF_Les_skin_boundary.roi

python -m tsp regionmembership --f 'CF_Les_CD*_masks_d2b.csv' --regionroi [CF_Les_normal_side.roi,CF_Les_lesion_side.roi]
