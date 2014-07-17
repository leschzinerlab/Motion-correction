Motion Correction for Dose-Fractionation Stack
Version 1.0, April 11, 2013

Features:
GPU enabled. The computations are accelerated by GPU. Processing a stack with 20
~30frames only takes ~20s, and usually less than 1 minute for more frames, hence
 is suit for on-the-fly processing.


System requirements:
1) NVIDIA GPU with large device memory. The memory size limits the maximum image
   size which can be Fourier transformed. Current K2 camera(Gatan) has nearly 
   7.5k*7.5k pixles output at super-resolution mode, it requires ~1.5Gb device 
   memory.For a 4k*4k image, it requires ~1Gb device memory. we had tested this
   program on several GPUs: GTX 480, GTX 580, GTX 690, and Tesla C2070.
2) Avilable system memory must be at least ~1.5x larger than the image stack.
3) Linux operating system, such as Fedora and Centos. gcc/g++ 4.1 or higher
   version must be installed. 
4) The newest NVIDIA CUDA driver and CUDA toolkit v5. 


Usage:
Running the program without following any parameters will give the help information.
Except the input stack filename, all other parameters are optional, which control 
the processing parameters and outputs. The outputs include corrected/uncorrected
frame sum, aligned image stack, and high-contrast diagnosis images(output into 
dosef_quick folder). The program has a build-in output filename convention, or
you can use your convention by specifying them in the parameters. Because the
algorithm is designed to handle very weak signal in the subframe, it is strongly
recommended to dose-fractionate an exposure to more subframes(larger than 10, 
typically 30~40 for 30e/A^2 total dose) for better correction.


The Motion Correction program was written by Xueming Li at Yifan Cheng Laboratory
at UCSF. For more details, source code and future updates, please visit the lab 
website at http://cryoem.ucsf.edu , or contact the author directly.

The Motion Correction program are licensed under the terms of the GNU Public 
License version 3 (GPLv3).

