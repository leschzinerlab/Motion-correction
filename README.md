# Motion correction

The scripts within this repository utilize both two different approachs for aligning movies collected with a K2 Summit direct electron detector: 

* Whole movie frame alignment:
  * runMotionCorrection.py
     * Which runs dosefgpu_driftcorr published by [Li et al. 2013]  (http://www.ncbi.nlm.nih.gov/pubmed/23644547) 

* Whole movie frame alignment with sub-frame alignment:
  * runMotionCorr2.py
     * Which runs MotionCor2 published by [Zheng et al. 2015] (http://biorxiv.org/content/early/2016/07/04/061960)

* Per-particle alignment: 
  * runLMBFGS_relion.py
     * Which runs alignparts_lmbfgs.exe published by [Rubinstein & Brubaker 2015] (http://arxiv.org/abs/1409.6789)

__Table of contents:__

1. [Software dependencies] (https://github.com/leschzinerlab/Motion-correction#software-dependencies) 
2. [Whole frame alignment using runMotionCorrection.py] (https://github.com/leschzinerlab/Motion-correction#whole-frame-alignment-using-runmotioncorrectionpy)
3. [Whole frame alignment with sub-frame alignment using runMotionCorr2.py] (https://github.com/leschzinerlab/Motion-correction#whole-frame-and-sub-frame-alignment-using-motioncorr2)
4. [Per-particle alignment using runLMBFGS_relion.py] (https://github.com/leschzinerlab/Motion-correction#per-particle-alignment-using-runlmbfgs_relionpy)

## Software dependencies

#### runMotionCorrection.py (dosefgpu_motioncorr)

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-5.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

<pre>PATH: /usr/local/cuda-5.0/bin
LD_LIBRARY_PATH: /usr/local/cuda-5.0/lib64</pre>

#### runMotionCorr2.py (MotionCorr2)

You must have a GPU processor installed on your workstation with at least 4 GB of GPU memory (e.g. NVIDIA Tesla K10).

CUDA-7.0 libraries must be installed on the GPU and be placed into your path before you can execute this program:

<pre>PATH: /usr/local/cuda-7.0/bin
LD_LIBRARY_PATH: /usr/local/cuda-7.0/lib64</pre>

#### runLMBFGS_relion.py (alignparts_lmbfgs.exe)

Within this repository we are including a pre-compiled alignparts_lmbfgs.exe program that was compiled on CentOS 7. To learn more about downloading the fortan source code and compiling, check out: 

[Rubinstein Lab Google Sites Scripts Repository] (https://sites.google.com/site/rubinsteingroup/direct-detector-align_lmbfgs)

Minimally, it will require: 
* Compiled versions of lm-bfgs from the Rubinstein lab
* EMAN2
* Gnuplot
 
## Whole frame alignment using runMotionCorrection.py

### Inputs

As the program is currently written, it will normalize and align UN-normalized movies that were collected using Leginon. This means that it expects the following inputs: 

* Movies with the extension '.mrcs'
* Gain reference '.mrc' file that has the same dimensions as the movies

### Running the program

<pre>$  Motion-correction/runMotionCorr.py 
Usage: runMotionCorr.py --dir=<folder with mrc frames> --gain_ref=<gain reference in mrc format with full path;input the *_norm* file from the leginon reference directory> --save_bin <save binned mic> --save_norm <save normalized frames>

This program takes movies with .mrcs extensions and creates aligned movies with .mrc extension, along with the option to create normalized movies with the .mrcs extension.

Options:
  -h, --help       show this help message and exit
  --dir=FILE       Directory containing direct detector movies (.mrcs extension)
  --gain_ref=FILE  Gain reference file from Leginon with the full path (.mrc)
  --save_bin       Save binned image for quick inspection
  --save_norm      Save normalized movie frames as _norm.mrcs
  --bin=INT        Binning factor to use during movie alignment, 1 or 2. (Default=1)
  -d               debug</pre>
  
Running notes: 

* To run this program, make sure that only the movies you want aligned are within a specified directory and have the '.mrcs' extension. NOTE: This program will run on ALL .mrcs files in the directory provided
* Specify the absolute paths to the directory with movies and to the gain reference .mrc file

Example command: 

<pre>$ Motion-correction/runMotionCorr.py --dir=/data/frames/leginon/15sep30a/rawdata/ --gain_ref=/data/frames/leginon/15sep30a_ref/rawdata/15sep30a_31115207_07_7676x7420_norm_1.mrc --bin=2</pre>

### Outputs

The program will then align each movie using its GPU cores to produce files with the extension '.mrc'. 

For example, if you had the input micrograph: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en.frames.mrcs</pre>

You would get: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en.frames.mrc</pre>

And, if you asked for normalized movie frames as an output with the <pre>--save_norm</pre> option,  you would also get: 

<pre>/data/frames/leginon/15sep30a/rawdata/15sep30a_b1_1e_00009gr_00003sq_00007hl_00001en_norm.mrcs</pre>

## Whole frame and sub-frame alignment using MotionCorr2

This program can align movies in addition to gain correcting raw movies before movie alignment. As described by the publication, MotionCor2 will also perform sub-movie frame alignment, as specified by the user. Defaults currently suggested by Zheng et al. are incorporated as defaults in this python wrapper.

The wrapper script runMotionCorr2.py will provide users with the ability to run the program over all movies within a given directory, or within a given text file. The program will also skip any movies that have output files already created.

### Inputs

This program will run on two types of movies: 
* Gain-corrected movies - These movies must have .mrcs extension
* Raw movies - These movies must have .frames.mrc extension along with a gain reference image from Leginon.

### Running the program

<pre>$ Motion-correction/runMotionCorr2.py 
Usage: runMotionCorr2.py --dir=[folder with mrc frames] --gain_ref=[gain reference in mrc format with full path;input the _norm file from the leginon reference directory] --save_bin [save binned mic] 

This program takes movies with .mrcs or .frames.mrc extensions and will create aligned movies with -a.mrc extension.
 If dose weighting, output aligned movie will be named -a_weighted.mrc


Options:
  -h, --help        show this help message and exit
  --microlist=FILE  Provide list of movies instead of directory location
  --dir=FILE        Or provide directory containing direct detector movies
                    (.mrcs extension unless using gain reference, then
                    .frames.mrc)
  --gain_ref=FILE   Gain reference file from Leginon with the full path (.mrc)
  --throw=INT       Number of initial frames to discard from alignment.
                    (Default=2)
  --binning=INT     Scaling factor for output aligned movie (Default=2)
  --patchsize=INT   Number of patches for local alignment (Default=5, which is
                    a 5 x 5 tiling)
  --dose=FLOAT      Optional: Input dose rate for dose weighting (electrons
                    per Angstrom-squared
  --kev=INT         Optional: If dose weighting, provide accelerating voltage
                    in keV
  --apix=FLOAT      Optional: If dose weighting, provide pixel size in
                    Angstroms/pixel
  -d                debug</pre>
  
  Users can provide either a directory or an input file. 
  
The program will decide if the input files are .mrcs or .frames.mrc based on whether the user input a gain reference. If it is a gain reference then it looks for .frames.mrc files, otherwise it will look for .mrcs.

**Dose weighting:** Users can dose weight their data by providing the dose rate, accelerating voltage, and pixel size. Otherwise the program will not weight the frames according to dose.

## Per-particle alignment using runLMBFGS_relion.py

In order to use alignparts_lmbfgs.exe, there is a certain workflow that you will need to follow: 

* Align direct detector movies
* Pick particles
* Extract particles using Relion
* Run alignparts_lmbfgs.exe to track individual particle trajectories

### Inputs

The script described below assumes that you have the following set up (based upon Relion directory and file formatting):

* Micrographs/ - directory containing aligned movie frames (.mrc), normalized UNALIGNED movies with .mrcs extension, CTF log files, and particle coordinates 
* Particles/Micrographs - directory containing extracted particles by Relion, which includes particle stacks per micrograph and associated STAR files.
* [input].star - STAR file generated by Relion during particle extraction

### Running the program

<pre>Usage: runLMBFGS_relion.py --star=<relion star> --radius=<radius>

Options:
  -h, --help            show this help message and exit
  --star=FILE           Relion starfile with particles. IMPORTANT: Particles
                        must be the same pixel size as the movies for this
                        program to work.
  --radius=INTEGER      Radius of particles in pixels
  --trialrun            Flag to run lm-bfgs on a single movie to check
                        particle trajectories and to show vector field plot
                        (Default=False)
  --overwrite           Flag to over write any existing LM-BFGS runs
                        (Default=False)
  --invert=INTEGER      Indicate whether to invert contrast (=1) or not (=0)
                        for output particle stack (Default=1)
  --exepath=PATH        Optional: Path to executable files. (Default=Motion-
                        correction/lm-bfgs_v3.0/)
  --movieNAME=Movie extension
                        Optional: Additional name for movies. (Default=_movie)
  --movieEXT=Movie extension
                        Optional: Movie extension. (Default=mrcs)
  --firstFrame=INTEGER  Optional: First frame of movies to use for alignment.
                        (Default=1)
  --lastFrame=INTEGER   Optional: Last frame of movies to use for alignment.
                        (Default=Last Frame)
  --smoothening=STRING  Optional: Amount of smoothening forced onto
                        trajectories of particles. (Default=1.0d4)
  --localsigma=INTEGER  Optional: Local sigma factor to increase correlation
                        between trajectories. (Default=500)
  --exaggerate=INTEGER  Optional: Factor by which particle trajectories should
                        be exaggerated in vector file. (Default=5)
  --apix=FLOAT          Optional: Provide pixel size instead of calculating
                        from star file.
  --exposureweight      Optional: Flag to exposureweight particles based upon
                        dose. (Default=False)
  --dose=FLOAT          IF EXPOSURE WEIGHTING: Dose per frame in electrons per
                        Angstroms-squared.
  --moviedimx=INT       Optional: Input movie dimensions - X axis. (By default
                        this is read from input file)
  --moviedimy=INT       Optional: Input movie dimensions - Y axis. (By default
                        this is read from input file)
  --maxframes=INT       Optional: Input maximum number of movie frames (By
                        default this is read from input file)
  --boxsize=INT         Optional: Input box size for particles(By default this
                        is read from input file)
  -d                    debug</pre>
  
 ___Required inputs:___  
  * --star - STAR file from Relion particle extraction
  * --radius - radius of particle in pixels
  * --trialrun - IMPORTANT: As described on [Rubinstein Lab Google site] (https://sites.google.com/site/rubinsteingroup/direct-detector-align_lmbfgs), they recommend check the smoothened trajectories to check if the program is work satisfactorily. To check this on a single micorgraph, include this flag. At the end of the particle trajectory alignment, you will be shown a plot of particle trajectories.
  
 ___Optional inputs:___
  * --overwrite - Ignores previous LM-BFGS run output files and removes them, allowing users to bypass the normal checks and error messages that would prevent the user from overwriting existing data.
  * --execpath - Input absolute path to compiled program 'alignparts_lmbfgs.exe'. By default, it will look in the folder provided by this Github repo: Motion-correction/lm-bfgs_v3.0/. 
  * --movieName - Additional name for movies. This is NOT the movie file extension, just additional text that differentiates the movie (if at all) from the aligned micrograph.
  * --movieEXT - Movie file extension (e.g. .mrcs)
  * --firstFrame - First frame to use in the movie alignment
  * --lastFrame - last frame to use in alignment. By default this will be the last frame in the movie.
  * --smoothening - parameter for smoothening function. NEEDS TO BE CHANGED IF YOU WANT TO ALTER OUTPUT PARTICLE TRAJECTORIES. 
  * --localsigma - parameter for enforcing local correlations for trajectories. NEEDS TO BE CHANGED IF YOU WANT TO ALTER OUTPUT PARTICLE TRAJECTORIES.
  * --exaggerate - can be changed to alter output trajectories. 
  * --apix - override pixel size calculated from STAR file
  * --exposureweight - flag to weight each frame of the movie according to the dose-dependent decay of resolution information 
  * --dose - if exposure weighting, include the dose on each frame of the movie in electrons per Angstroms-squared per frame.
  * --moviedimx - Movie dimensions on the X-axis. Allows user to override automatic reading of dimensions with e2iminfo.py within routine.
  * --moviedimy - Movie dimensions on the Y-axis. Allows user to override automatic reading of dimensions with e2iminfo.py within routine. 
  * --maxframes - Maximum number of frames in movies. Allows user to override automatic reading of dimensions with e2iminfo.py within routine.
  * --boxsize - Box size for particles. Allows user to override automatic reading of dimensions with e2iminfo.py within routine.

Example command: 

<pre>Motion-correction/runLMBFGS_relion.py --star=particles.star --radius=50</pre>

### Outputs

The program will provide the following outputs: 

* [input]_lmbfgs.star - output star file containing paths to the particles aligned per trajectory
* Particles/Micrographs/[file name]_lmbfgs.mrcs - particle stack for aligned particles
* Particles/Micrographs/[file name]_lmbfgs.vec - vector file to plot per particle trajectories 
