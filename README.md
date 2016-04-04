# Motion correction

The scripts within this repository utilize both two different approachs for aligning movies collected with a K2 Summit direct electron detector: 

* Whole movie frame alignment:
  * runMotionCorrection.py
     * Which runs dosefgpu_driftcorr published by [Li et al. 2013]  (http://www.ncbi.nlm.nih.gov/pubmed/23644547) 

* Per-particle alignment: 
  * runLMBFGS_relion.py
     * Which runs alignparts_lmbfgs.exe published by [Rubinstein & Brubaker 2015] (http://arxiv.org/abs/1409.6789)

__Table of contents:__

1. [Software dependencies] (https://github.com/leschzinerlab/Motion-correction#software-dependencies) 
2. [Whole frame alignment using runMotionCorrection.py] (https://github.com/leschzinerlab/Motion-correction#whole-frame-alignment-using-runmotioncorrectionpy) 
3. [Per-particle alignment using runLMBFGS_relion.py] (https://github.com/leschzinerlab/Motion-correction#per-particle-alignment-using-runlmbfgs_relionpy)

## Software dependencies

#### runMotionCorrection.py (dosefgpu_motioncorr)

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-5.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

<pre>PATH: /usr/local/cuda-5.0/bin
LD_LIBRARY_PATH: /usr/local/cuda-5.0/lib64</pre>

#### runLMBFGS_relion.py (alignparts_lmbfgs.exe)

Within this repository we are including a pre-compiled alignparts_lmbfgs.exe program that was compiled on CentOS 7. To learn more about downloading the fortan source code and compiling, check out: 

[Rubinstein Lab Google Sites Scripts Repository] (https://sites.google.com/site/rubinsteingroup/direct-detector-align_lmbfgs)

## Whole frame alignment using runMotionCorrection.py

### Inputs

As the program is currently written, it will normalize and align UN-normalized movies that were collected using Leginon. This means that it expects the following inputs: 

* Movies with the extension '.frames.mrcs'
* Gain reference '.mrc' file that has the same dimensions as the movies

### Running the program

<pre>$  Motion-correction/runMotionCorr.py 
Usage: runMotionCorr.py --dir=<folder with mrc frames> --gain_ref=<gain reference in mrc format with full path;input the *_norm* file from the leginon reference directory> --save_bin <save binned mic> --save_norm <save normalized frames>

This program takes movies with .frames.mrcs extensions and creates aligned movies with .mrc extension, along with the option to create normalized movies with the .mrcs extension.

Options:
  -h, --help       show this help message and exit
  --dir=FILE       Directory containing direct detector movies (.frames.mrcs
                   extension)
  --gain_ref=FILE  Gain reference file from Leginon with the full path (.mrc)
  --save_bin       Save binned image for quick inspection
  --save_norm      Save normalized movie frames as .mrcs
  --bin=INT        Binning factor to use during movie alignment, 1 or 2.
                   (Default=1)
  -d               debug</pre>
  
Running notes: 

* To run this program, make sure that only the movies you want aligned are within a specified directory and have the '.frames.mrcs' extension.
* Specify the absolute paths to the directory with movies and to the gain reference .mrc file

Example command: 

<pre>$ Motion-correction/runMotionCorr.py --dir=/data/frames/leginon/15sep30a/rawdata/ --gain_ref=/data/frames/leginon/15sep30a_ref/rawdata/15sep30a_31115207_07_7676x7420_norm_1.mrc --bin=2</pre>

### Outputs

The program will then align each movie using its GPU cores to produce files with the extension '.mrc'. 

For example, if you had the input micrograph: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en.frames.mrcs</pre>

You would get: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en.mrc</pre>

And, if you asked for normalized movie frames as an output with the <pre>--save_norm</pre> option,  you would also get: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en.mrcs</pre>
Work in progress

## Per-particle alignment using runLMBFGS_relion.py

Work in progress
