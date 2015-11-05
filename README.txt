#Motion correction

Python scripts for using 'dosefgpu_driftcorr' alignment procedure for direct electron
detector movies.

##Software dependencies

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-5.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

   PATH: /usr/local/cuda-5.0/bin
   LD_LIBRARY_PATH: /usr/local/cuda-5.0/lib64

##runMotionCorr.py

This is a python wrapper script that will run this motion correction software. 
