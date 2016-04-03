# Motion correction

The scripts within this repository utilize both whole frame alignments (runMotionCorrection.py --> dosef_motioncorr) [Li et al. 2013] (http://www.ncbi.nlm.nih.gov/pubmed/23644547) as well as per-particle movie alignment (runLMBFGS_relion.py --> lm-bfgs) (Source). 

While both of these programs utilize underlying GPU or fortan code, we have written python wrappers to improve the workflow in using these programs.

# Software dependencies

## dosef_motioncorr

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-5.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

   PATH: /usr/local/cuda-5.0/bin
   LD_LIBRARY_PATH: /usr/local/cuda-5.0/lib64

## LM-BFGS

##runMotionCorr.py

This is a python wrapper script that will run this motion correction software. 

