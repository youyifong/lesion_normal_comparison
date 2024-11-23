# Disease Severity Correlates

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

## Methods 

- Define an environment variable XX to denote the biopsy folder, e.g.,
```{bash}
export XX=AS_Les
cd ${XX}
```

- Perform cell segmentation

- Compute distance to skin and to DEJ for each mask 
```{bash}
python -m tsp dist2boundary
```

- Run R scripts perform statistical analyses and generate plots.
