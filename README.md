# Motion correction

The scripts within this repository utilize both two different approachs for aligning movies collected with a K2 Summit direct electron detector: 

* Whole movie frame alignment:
  * runMotionCorrection.py
     * Which runs dosefgpu_driftcorr published by [Li et al. 2013]  (http://www.ncbi.nlm.nih.gov/pubmed/23644547) 

* Per-particle alignment: 
  * runLMBFGS_relion.py
     * Which runs alignparts_lmbfgs.exe published by [Rubinstein & Brubaker 2015] (http://arxiv.org/abs/1409.6789)

__Table of contents:__

1. Software dependencies 
2. Whole frame alignment using runMotionCorrection.py
3. Per-particle alignment using runLMBFGS_relion.py

## Software dependencies

#### runMotionCorrection.py (dosefgpu_motioncorr)

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-5.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

```PATH: /usr/local/cuda-5.0/bin```
```LD_LIBRARY_PATH: /usr/local/cuda-5.0/lib64```

#### runLMBFGS_relion.py (alignparts_lmbfgs.exe)

Work in progress

## Whole frame alignment using runMotionCorrection.py

Work in progress

## Per-particle alignment using runLMBFGS_relion.py

Work in progress
