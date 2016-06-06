#!/usr/bin/env python 

from __future__ import division
import shutil
import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache
import time


#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog --star=<relion star> --radius=<radius>")
        parser.add_option("--star",dest="starfile",type="string",metavar="FILE",
                help="Relion starfile with particles. IMPORTANT: Particles must be the same pixel size as the movies for this program to work.")
        parser.add_option("--radius",dest="radius",type="int",metavar="INTEGER",
                help="Radius of particles in pixels")
	parser.add_option("--trialrun",action='store_true',dest="trialrun",default=False,
                help="Flag to run lm-bfgs on a single movie to check particle trajectories and to show vector field plot (Default=False)")
	parser.add_option("--invert",dest="invertcontrast",default=1,type="int",metavar="INTEGER",
                help="Indicate whether to invert contrast (=1) or not (=0) for output particle stack (Default=1)")
	parser.add_option("--nprocs",dest="nprocs",type="int",metavar="INTEGER",default=1,
                help="Number of CPUs for parallelization. (Default=1)")
	parser.add_option("--exepath",dest="exepath",type="string",metavar="PATH",default='lm-bfgs_v3.0/',
                help="Optional: Path to executable files. (Default=Motion-correction/lm-bfgs_v3.0/)")
	parser.add_option("--movieNAME",dest="movieEXT1",type="string",metavar="Movie extension",default='.frames',
                help="Optional: Additional name for movies. (Default=.frames)")
	parser.add_option("--movieEXT",dest="movieEXT2",type="string",metavar="Movie extension",default='mrcs',
                help="Optional: Movie extension. (Default=mrcs)")
	parser.add_option("--firstFrame",dest="firstframe",type="int",metavar="INTEGER",default=1,
                help="Optional: First frame of movies to use for alignment. (Default=1)")
	parser.add_option("--lastFrame",dest="lastframe",type="int",metavar="INTEGER",default=0,
                help="Optional: Last frame of movies to use for alignment. (Default=Last Frame)")
        parser.add_option("--smoothening",dest="smooth",type="string",metavar="STRING",default='1.0d4',
                help="Optional: Amount of smoothening forced onto trajectories of particles. (Default=1.0d4)")
	parser.add_option("--exaggerate",dest="exaggerate",type="int",metavar="INTEGER",default='5',
                help="Optional: Factor by which particle trajectories should be exaggerated in vector file. (Default=5)")
	parser.add_option("--apix",dest="pixelsize",type="float",metavar="FLOAT",
                help="Optional: Provide pixel size instead of calculating from star file.")
	parser.add_option("--exposureweight",action='store_true',dest="exposureweight",default=False,
                help="Optional: Flag to exposureweight particles based upon dose. (Default=False)")
	parser.add_option("--dose",dest="dose",type='float',metavar='FLOAT',default=1,
                help="IF EXPOSURE WEIGHTING: Dose per frame in electrons per Angstroms-squared.")
	parser.add_option("--moviedimx",dest="moviedimx",type="int",metavar="INT",default=-1,
                help="Optional: Input movie dimensions - X axis. (By default this is read from input file)")
	parser.add_option("--moviedimy",dest="moviedimy",type="int",metavar="INT",default=-1,
                help="Optional: Input movie dimensions - Y axis. (By default this is read from input file)")
	parser.add_option("--maxframes",dest="maxframes",type="int",metavar="INT",default=-1,
                help="Optional: Input maximum number of movie frames (By default this is read from input file)")
	parser.add_option("--boxsize",dest="boxsize",type="int",metavar="INT",default=-1,
                help="Optional: Input box size for particles(By default this is read from input file)")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#=============================
def checkConflicts(params):

        if not os.path.exists(params['starfile']):
                print "\nError: Starfile %s doesn't exists, exiting.\n" %(params['starfile'])
                sys.exit()

	if os.path.exists('%s_lmbfgs.star' %(params['starfile'][:-5])):
		print "\nError: Output starfile %s_lmbfgs.star already exists. Exiting\n" %(params['starfile'][:-5])
		sys.exit()

	if not os.path.exists('Micrographs/'):
		print "\nError: Could not find folder containing movies and micrographs 'Micrographs/'. Exiting\n"
		sys.exit()

	if not os.path.exists('Particles/Micrographs'):
                print "\nError: Could not find folder containing particles 'Particles/Micrographs'. Exiting\n"
                sys.exit()

	if os.path.exists('align_lmbfgs.bash'):
		print "\nError: script align_lmbfgs.bash already exists. Exiting."
		sys.exit()

	if os.path.exists('coord.txt'):
		print "\nError: coordinate text file coord.txt already exists. Exiting."
		sys.exit()

	if os.path.exists('movie.txt'):
		print "\nErorr: movie text file movie.txt already exists. Exiting."
		sys.exit()

#==============================
def getPath(pathcheck,executable):

	if os.path.exists('%s/%s'%(pathcheck,executable)):
		return '%s/%s'%(pathcheck,executable)

#==============================
def checkExists(starfile,debug):

	#get relion indices for microraph and particles
	microcol=int(getRelionColumnIndex(starfile,'_rlnMicrographName'))-1
	particlecol=int(getRelionColumnIndex(starfile,'_rlnImageName'))-1
	flag=0

	basename=len(starfile.split('.')[0])

	for line in open(starfile,'r'):

		if len(line) < 50:
			continue

		micro=line.split()[microcol]
		part=line.split()[particlecol]

		if debug is True:
			print line
			print micro	
			print part

		if os.path.exists('%s.coord' %(micro[:-4])):
			print '\nError: LM-BFGS file %s.coord already exists.' %(micro[:-4])
			flag=1

		if os.path.exists('%s_lmbfgs.mrcs' %(part.split('@')[-1][:-(5+basename+1)])):
			print '\nError: LM-BFGS output particles already exist: %s_lmbfgs.mrcs' %(part.split('@')[-1][:-5])
			flag=1

		if os.path.exists('%s_lmbfgs.vec' %(part.split('@')[-1][:-(5+basename+1)])):
                        print '\nError: LM-BFGS output particle vector already exists: %s_lmbfgs.vec' %(part.split('@')[-1][:-5])
                        flag=1
		
		return flag
#===============================
def getRelionColumnIndex(star,rlnvariable):

    counter=50
    i=1

    while i<=50:

        line=linecache.getline(star,i)

        if len(line)>0:
            if len(line.split())>1:
                if line.split()[0] == rlnvariable:
                    return line.split()[1][1:]

        i=i+1

#================================
def splitSTAR(starfile,nprocs,debug):

	#Check if other files exist
	i=1
	while i<=nprocs:
		if os.path.exists('%s_set%i.star' %(starfile[:-5],nprocs)):
			os.remove('%s_set%i.star' %(starfile[:-5],nprocs))	
		i=i+1

	#Get number of micrographs
	counter=0
	currentMicro='blank.mrc'
	for line in open(starfile,'r'):
		if len(line)>50:
		
			if line.split()[0] != currentMicro:
				currentMicro=line.split()[0]
				counter=counter+1

	if debug is True:
		print 'Number of micrographs = %i' %(counter)

	if counter < nprocs:
		nprocs=counter

	if debug is True:
		print 'Nprocs=%i' %(nprocs)
	
	#Check if even division:
	if counter % nprocs == 0: 
		group=float(counter)/float(nprocs)

	if counter % nprocs != 0:
		group=counter/nprocs + nprocs

	if debug is True:
		print 'Splitting star file %s into %i groups with %i micrographs in each file' %(starfile,nprocs,int(group))

	groupcount=1
	grouptot=nprocs

	while groupcount <= grouptot:
		if debug is True:
			print 'Working on creating group %i' %(groupcount)

		#Open file & write header, keeping track of number o flines in header
		o1=open('%s_set%i.star' %(starfile[:-5],groupcount),'w')

		headerlength=0
		totlines=0
		for line in open(starfile,'r'):

			if len(line)<50:
				o1.write(line)
				headerlength=headerlength+1

			totlines=1+totlines

		counterpart=1
		currentMicro='blank.mrc'
		countermicro=0
		begingroupnum=float(group)*(groupcount)-(float(group))+1
		endgroupnum=float(group)*(groupcount)

		if debug is True:
			print 'Totlines=%i' %(totlines)
			print 'Begin group num=%f' %(begingroupnum)
			print 'End group num=%f' %(endgroupnum)

		while counterpart < totlines:

			line=linecache.getline(starfile,counterpart)

			if len(line) < 50:
				counterpart=counterpart+1
				continue

			micro=line.split()[0]

			if micro != currentMicro:
				if debug is True:
					print 'Micro query %s does not match currentMicro %s' %(micro, currentMicro) 
					print 'Incrementing countermicro from %i to %i' %(countermicro,countermicro+1)
				currentMicro=micro	
				countermicro=countermicro+1

			if countermicro >= begingroupnum:
				if countermicro < endgroupnum:
					o1.write(line)
			
			counterpart=counterpart+1

		o1.close()

		groupcount=groupcount+1

	return nprocs	

#=========================
def combineSTARfiles(basename,nprocs,outfile):

	nproc=1

	o1=open(outfile,'w')
	#write header
	for line in open('%s%i.star' %(basename,1)):
		if len(line)<50:
			o1.write(line) 

	while nproc<=nprocs:

		for line in open('%s%i.star' %(basename,nproc)):
			if len(line)>50:
				o1.write(line)
		
		nproc=nproc+1

#==============================
def getEMANPath():
        ### get the eman2 directory        
        emanpath = subprocess.Popen("env | grep EMAN2DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()

        if emanpath:
                emanpath = emanpath.replace("EMAN2DIR=","")
        if os.path.exists(emanpath):
                return emanpath
        print "EMAN2 was not found, make sure it is in your path"
        sys.exit()

#==============================
def getGNUPLOTPath():
        ### get the eman2 directory        
        emanpath = subprocess.Popen("which gnuplot", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if emanpath.split('/')[-1] == 'gnuplot':
                return emanpath
        print "GNUPLOT was not found, please install in order to visualize trajectories."
        sys.exit()


#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
	getEMANPath()
	getGNUPLOTPath()
		
	checkexec='%s/lm-bfgs_v3.0/' %('/'.join(sys.argv[0].split('/')[:-1]))

	if params['exepath'] != 'lm-bfgs_v3.0/':	
	        lmbfgs=getPath(params['exepath'],'alignparts_lmbfgs.exe')
		lmbstar=getPath(params['exepath'],'alignparts_starfilehandler.exe')

	if params['exepath'] == 'lm-bfgs_v3.0/':
		lmbfgs=getPath(checkexec,'alignparts_lmbfgs.exe')
                lmbstar=getPath(checkexec,'alignparts_starfilehandler.exe')

	if not lmbfgs or not lmbstar:
		print 'Error: Could not find executables alignparts_lmbfgs.exe and alignparts_starfilehandler.exe in %s' %(params['exepath'])
		sys.exit()

	if params['debug'] is True:
		print lmbfgs
		print lmbstar

	#Check if any outputs exist in the Micrographs and Particles/Micrographs folders:
	checkFlag=checkExists(params['starfile'],params['debug'])
	if checkFlag ==1:
		print '\nError: one or more LM-BFGS files already exist. Please remove these files, and resubmit.\n'
		sys.exit()

	#Get dimensions of movies
	if params['debug'] is True:
		print 'Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2'])	
		print glob.glob('Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2']))[0] 
	
	if params['moviedimx'] <0:
		iminfo=subprocess.Popen("e2iminfo.py %s" %(glob.glob('Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2']))[0]), shell=True, stdout=subprocess.PIPE).stdout.read().strip().split('\t')[-1].split('\n')[0].split('x')
		if params['debug'] is True:
			print iminfo	
	if params['moviedimy'] <0:
                iminfo=subprocess.Popen("e2iminfo.py %s" %(glob.glob('Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2']))[0]), shell=True, stdout=subprocess.PIPE).stdout.read().strip().split('\t')[-1].split('\n')[0].split('x')
                if params['debug'] is True:
                        print iminfo

	if params['maxframes'] <0:
                iminfo=subprocess.Popen("e2iminfo.py %s" %(glob.glob('Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2']))[0]), shell=True, stdout=subprocess.PIPE).stdout.read().strip().split('\t')[-1].split('\n')[0].split('x')
                if params['debug'] is True:
                        print iminfo

	if params['moviedimx'] < 0:
		movieDIMX=iminfo[0]
	if params['moviedimx'] > 0:
		movieDIMX=params['moviedimx']

	if params['moviedimy'] < 0:
                movieDIMY=iminfo[1]
        if params['moviedimy'] > 0:
                movieDIMY=params['moviedimy']

	if params['maxframes'] <0:	
		if len(iminfo) == 2:
			maxframe=subprocess.Popen("e2iminfo.py %s" %(glob.glob('Micrographs/*%s.%s'%(params['movieEXT1'],params['movieEXT2']))[0]), shell=True, stdout=subprocess.PIPE).stdout.read().strip().split('\n')[1].split()[0]

		if len(iminfo) ==3:
			maxframe=iminfo[2]
	if params['maxframes'] >0:
		maxframe=params['maxframes']

	if params['debug'] is True:
		print movieDIMX
		print movieDIMY
		print maxframe
	
	#Get boxsize of particles
	try:
    		testMic=linecache.getline(params['starfile'],20).split()[3].split('@')[-1]
	except IndexError:
    		testMic=linecache.getline(params['starfile'],40).split()[3].split('@')[-1]
	if params['boxsize']<0:
		boxsize=subprocess.Popen("e2iminfo.py %s" %(testMic), shell=True, stdout=subprocess.PIPE).stdout.read().strip().split('\t')[-1].split('\n')[0].split('x')[0] 

	if params['boxsize']>0:
		boxsize=params['boxsize']
	if params['debug'] is True:
		print boxsize
	
	#First, last frames 

	if params['lastframe'] == 0:
		#Set last frame equal to max number 
		params['lastframe']=maxframe

	#Zero frame
	zeroframe=int(params['lastframe'])/2
	#======>Need to check if integer after dividing by 2

	#Get parameters from relion file
	voltage=int(getRelionColumnIndex(params['starfile'],'_rlnVoltage'))-1
	pixel=int(getRelionColumnIndex(params['starfile'],'_rlnDetectorPixelSize'))-1
	mag=int(getRelionColumnIndex(params['starfile'],'_rlnMagnification'))-1
	
	#Get random particle & grab voltage
	kev=int(float(linecache.getline(params['starfile'],100).split()[voltage]))
	if params['debug'] is True:
		print kev
	
	#Get random particle and get pixle size and mag
	detectorsize=linecache.getline(params['starfile'],100).split()[pixel]
	magnification=linecache.getline(params['starfile'],100).split()[mag]
	pixelSize=float(detectorsize)/float(magnification)*10000

	if params['debug'] is True:
		print detectorsize
		print magnification
		print pixelSize

	if params['pixelsize']:
		pixelSize=params['pixelsize']

	print 'Using pixel size of %f Angstroms/pixel' %(pixelSize)

	#Movie list
	movielist='movie.txt'
	coordlist='coord.txt'
	
	#Running for a trial run
	if params['nprocs'] == 1 or params['trialrun'] is True:

		cmd='#!/bin/bash\n'
		cmd+='# Things you need to adjust:\n'
		cmd+='#********************************************************************\n'
		cmd+='# locations of executables provided by user\n'
		cmd+='alignparts_starfilehandler=%s\n' %(lmbstar)
		cmd+='alignparts_lmbfgs=%s\n' %(lmbfgs)
		cmd+='#********************************************************************\n'
		cmd+='# Input and output files, paths, extensions\n'
		cmd+='instar=%s             # input star file name\n' %(params['starfile'])
		cmd+='outstar=%s_lmbfgs.star     # output star file name\n' %(params['starfile'][:-5])
		cmd+='moviepath=Micrographs/                                # directory where movies are located (input)\n' 
		cmd+='particlepath=Particles/Micrographs/                   # directory where particles are located (output)\n'
		cmd+='movieflag=%s                                 # what to add to the micrograph name to get the movie name (use \"\" to indicate no flag)\n' %(params['movieEXT1'])
		cmd+='movieext=%s                                    # what to change the micrograph extension to to get the movie name\n' %(params['movieEXT2'])
		cmd+='#********************************************************************\n'
		cmd+='# particle and movie information \n'
		cmd+='boxsize=%i                    # boxsize for particles (in pixels)\n' %(int(boxsize))
		cmd+='particleradius=%i              # radius for particles within box (in pixels)\n' %(params['radius'])
		cmd+='pixelsize=%f                 # pixel size (in Angstroms)\n' %(pixelSize)
		cmd+='framex=%i                    # x-dimension of movies (in pixels)\n' %(int(movieDIMX))
		cmd+='framey=%i                    # y-dimension of movies (in pixels)\n' %(int(movieDIMY))
		cmd+='#********************************************************************\n'
		cmd+='# frame information provided\n'
		cmd+='framefirstali=%i                # first frame to be used in alignment\n' %(params['firstframe'])
		cmd+='framelastali=%i                # last frame to be used in alignment\n' %(int(params['lastframe']))
		cmd+='framefirstave=%i                # first frame to be used in average of frames\n' %(params['firstframe'])
		cmd+='framelastave=%i                # last frame to be used in average of frames\n'%(int(params['lastframe']))
		cmd+='#********************************************************************\n'
		cmd+='# trajectory smoothing information provided by user\n'
		cmd+='smooth=%s                   # specifies the amount of smoothing forced on trajectories of particles\n' %(params['smooth'])
		cmd+='zeroframe=%i                   # which frame is considers as unshifted (MUST MATCH THE ZEROFRAME FOR WHOLE FRAME ALIGNMENT)\n' %(zeroframe)
		cmd+='exaggerate=%i                   # factor by which particle trajectories should be exaggerated in vector file\n' %(params['exaggerate'])
		cmd+='invertoutput=%i                 # 1 inverts output from movie densities, 0 does not\n'%(params['invertcontrast'])
		cmd+='localavg=1                     # 1 performs local averaging of trajectories, 0 turns off this feature\n'
		cmd+='localavgsigma=500              # the standard deviation used to weight local averaging\n'
		cmd+='#********************************************************************\n'
		cmd+='# Exposure weighting (optional)\n'
		if params['exposureweight'] is True:
			weight=1
		if params['exposureweight'] is False:
			weight=0
		cmd+='expweight=%i                    # 1 turns on exposure weighting, 0 turns off exposure weighting\n' %(weight)
		cmd+='akv=%i                        # microscope accelerating voltage in kV\n' %(kev)
		cmd+='expperframe=%f                # Exposure per frame in electrons per Angstrom squared\n' %(params['dose'])
		cmd+='# Things you do not need to adjust:\n'
		cmd+='#=================================================================\n'
		cmd+='# Files and info generated by alignparts_starfilehandler and used by alignparts_lmbfgs\n'	
		cmd+='lmbfgsflag=_lmbfgs             # file name modifier for locally aligned particle stacks\n' 
		cmd+='lmbfgsext=mrcs                 # extension for locally aligned particle stacks\n'
		cmd+='movielist=%s           # name of list of movie files\n' %(movielist)
		cmd+='coordlist=%s           # name of list of coordinate files\n'%(coordlist)
		cmd+='vecext=vec                     # extension for output vector files\n'
		cmd+='# Alignparts_lmbfgs variables that do not usually need adjustment\n'
		cmd+='maxparts=100000                # maximum number of particles that may be in a micrograph\n'
		cmd+='nsigma=5                       # standard deviations above the mean to remove outlier pixels\n'
		cmd+='rmax1=500                      # Low resolution cutoff (in Angstroms) used for alignment\n'
		cmd+='rmax2=40                       # High resolution cutoff (in Angstroms) used for alignment\n'
		cmd+='bfactor=2000                   # B-factor (in A**2) used for alignment\n'
		cmd+='factr=1d7                      # Accuracy of minimizer (eg.1.0d7=1.0x10^7 in dble precision)\n'
		cmd+='#=====================================================================\n'
		cmd+='time $alignparts_starfilehandler << eot1\n'
		cmd+='$instar\n'	
		cmd+='$outstar\n'
		cmd+='$movielist\n'
		cmd+='$coordlist\n'
		cmd+='$boxsize,$framex,$framey\n'
		cmd+='$moviepath\n'
		cmd+='$particlepath\n'
		cmd+='$movieflag\n'
		cmd+='$movieext\n'
		cmd+='$lmbfgsflag\n'
		cmd+='$lmbfgsext\n'	
		cmd+='eot1\n'
		cmd+='\n'
		if params['trialrun'] is True:
			if os.path.exists('tmp12222.txt'):
				os.remove('tmp12222.txt')
			cmd+='mv $movielist tmp12222.txt\n'
			cmd+='head -1 tmp12222.txt > $movielist\n'
			if os.path.exists('tmp12223.txt'):
       	        	        os.remove('tmp12223.txt')
                	cmd+='mv $coordlist tmp12223.txt\n'
                	cmd+='head -1 tmp12223.txt > $coordlist\n'
			cmd+='rm tmp12222.txt\n'
			cmd+='rm tmp12223.txt\n'
		cmd+='time $alignparts_lmbfgs << eot2\n'
		cmd+='$movielist\n'
		cmd+='$coordlist\n'
		cmd+='$boxsize,0,$particleradius,$pixelsize,$nsigma,$rmax1,$rmax2\n'
		cmd+='$expweight,$akv,$expperframe\n'
		cmd+='$bfactor,$smooth,$exaggerate,$zeroframe,$invertoutput\n'
		cmd+='$localavg,$maxparts,$localavgsigma\n'
		cmd+='$framefirstali,$framelastali,$framefirstave,$framelastave\n'
		cmd+='$factr\n'
		cmd+='$moviepath\n'
		cmd+='$particlepath\n'
		cmd+='$lmbfgsflag\n'
		cmd+='$vecext\n'
		cmd+='$lmbfgsext\n'
		cmd+='eot2\n'

		run = open('align_lmbfgs.bash','w')
        	run.write(cmd)
        	run.close()

        	cmd2 = 'chmod +x align_lmbfgs.bash'
        	subprocess.Popen(cmd2,shell=True).wait()	

		print '\nRunning LM-BFGS ...\n'

		cmd3 = './align_lmbfgs.bash'
		subprocess.Popen(cmd3,shell=True).wait()

		if params['debug'] is False:
			os.remove('align_lmbfgs.bash')

		print '\nLM-BFGS finished.\n'

		if params['trialrun'] is True:
			print '\nPlotting vector field...\n'
			print '\nPress [Enter] to advance\n'

			#Get only vector file created
			particle=linecache.getline('coord.txt',1).split('/')[-1]	

			vec='set size square\n'
			vec+='set xrange [1:%i]\n' %(int(movieDIMX))
       		 	vec+='set yrange [1:%i]\n' %(int(movieDIMY))
			vec+='set palette defined (0 "black", 1 "red", 2 "orange", 3 "green", 4 "blue")\n'
			vec+='file1 = "Particles/Micrographs/%s_lmbfgs.vec"\n' %(particle[:-(7)])
			vec+='pointsize = 0.5\n'
			vec+='plot file1 u 2:3:1 w l lc palette\n'
			vec+='pause -1\n'

			if os.path.exists('gnuplotvectorfield.script'):
				os.remove('gnuplotvectorfield.script')

			run = open('gnuplotvectorfield.script','w')
	        	run.write(vec)
        		run.close()

        		cmd3 = 'gnuplot gnuplotvectorfield.script'
        		subprocess.Popen(cmd3,shell=True).wait()
	
			print '\nWriting post-script file of vector trajectories\n'
			print '\nPress [Enter] to finish\n'

			if os.path.exists("Particles/Micrographs/%s_vectorTrajectories.ps" %(particle[:-(7)])):
				os.remove("Particles/Micrographs/%s_vectorTrajectories.ps" %(particle[:-(7)]))

			vec='set term postscript\n'
        		vec+='set output "Particles/Micrographs/%s_vectorTrajectories.ps"\n' %(particle[:-(7)])
			vec+='set size square\n'
                	vec+='set xrange [1:%i]\n' %(int(movieDIMX))
                	vec+='set yrange [1:%i]\n' %(int(movieDIMY))
                	vec+='set palette defined (0 "black", 1 "red", 2 "orange", 3 "green", 4 "blue")\n'
                	vec+='file1 = "Particles/Micrographs/%s_lmbfgs.vec"\n' %(particle[:-(7)])
	                vec+='pointsize = 0.5\n'
	   	        vec+='plot file1 u 2:3:1 w l lc palette\n'
       	         	vec+='pause -1\n'
			vec+='set term x11\n'

                	if os.path.exists('gnuplotvectorfield.script'):
                        	os.remove('gnuplotvectorfield.script')

                	run = open('gnuplotvectorfield.script','w')
                	run.write(vec)
                	run.close()

                	cmd3 = 'gnuplot gnuplotvectorfield.script'
               		subprocess.Popen(cmd3,shell=True).wait()

			if params['debug'] is False:
				cmd='rm Micrographs/*.coord Particles/Micrographs/*.vec align_lmbfgs.bash movie.txt coord.txt %s_set*star' %(params['starfile'][:-5])
                        	subprocess.Popen(cmd2,shell=True).wait()
                        
	if params['trialrun'] is False and params['nprocs'] >1:

		cmd='#!/bin/bash\n'
		cmd+='# Things you need to adjust:\n'
		cmd+='#********************************************************************\n'
		cmd+='# locations of executables provided by user\n'
		cmd+='alignparts_starfilehandler=%s\n' %(lmbstar)
		cmd+='alignparts_lmbfgs=%s\n' %(lmbfgs)
		cmd+='#********************************************************************\n'
		cmd+='# Input and output files, paths, extensions\n'
		cmd+='instar=$1             # input star file name\n' 
		cmd+='outstar=${1%.*}_lmbfgs_set${2}.star     # output star file name\n' 
		cmd+='moviepath=Micrographs/                                # directory where movies are located (input)\n' 
		cmd+='particlepath=Particles/Micrographs/                   # directory where particles are located (output)\n'
		cmd+='movieflag=%s                                 # what to add to the micrograph name to get the movie name (use \"\" to indicate no flag)\n' %(params['movieEXT1'])
		cmd+='movieext=%s                                    # what to change the micrograph extension to to get the movie name\n' %(params['movieEXT2'])
		cmd+='#********************************************************************\n'
		cmd+='# particle and movie information \n'
		cmd+='boxsize=%i                    # boxsize for particles (in pixels)\n' %(int(boxsize))
		cmd+='particleradius=%i              # radius for particles within box (in pixels)\n' %(params['radius'])
		cmd+='pixelsize=%f                 # pixel size (in Angstroms)\n' %(pixelSize)
		cmd+='framex=%i                    # x-dimension of movies (in pixels)\n' %(int(movieDIMX))
		cmd+='framey=%i                    # y-dimension of movies (in pixels)\n' %(int(movieDIMY))
		cmd+='#********************************************************************\n'
		cmd+='# frame information provided\n'
		cmd+='framefirstali=%i                # first frame to be used in alignment\n' %(params['firstframe'])
		cmd+='framelastali=%i                # last frame to be used in alignment\n' %(int(params['lastframe']))
		cmd+='framefirstave=%i                # first frame to be used in average of frames\n' %(params['firstframe'])
		cmd+='framelastave=%i                # last frame to be used in average of frames\n'%(int(params['lastframe']))
		cmd+='#********************************************************************\n'
		cmd+='# trajectory smoothing information provided by user\n'
		cmd+='smooth=%s                   # specifies the amount of smoothing forced on trajectories of particles\n' %(params['smooth'])
		cmd+='zeroframe=%i                   # which frame is considers as unshifted (MUST MATCH THE ZEROFRAME FOR WHOLE FRAME ALIGNMENT)\n' %(zeroframe)
		cmd+='exaggerate=%i                   # factor by which particle trajectories should be exaggerated in vector file\n' %(params['exaggerate'])
		cmd+='invertoutput=%i                 # 1 inverts output from movie densities, 0 does not\n'%(params['invertcontrast'])
		cmd+='localavg=1                     # 1 performs local averaging of trajectories, 0 turns off this feature\n'
		cmd+='localavgsigma=500              # the standard deviation used to weight local averaging\n'
		cmd+='#********************************************************************\n'
		cmd+='# Exposure weighting (optional)\n'
		if params['exposureweight'] is True:
			weight=1
		if params['exposureweight'] is False:
			weight=0
		cmd+='expweight=%i                    # 1 turns on exposure weighting, 0 turns off exposure weighting\n' %(weight)
		cmd+='akv=%i                        # microscope accelerating voltage in kV\n' %(kev)
		cmd+='expperframe=%f                # Exposure per frame in electrons per Angstrom squared\n' %(params['dose'])
		cmd+='# Things you do not need to adjust:\n'
		cmd+='#=================================================================\n'
		cmd+='# Files and info generated by alignparts_starfilehandler and used by alignparts_lmbfgs\n'	
		cmd+='lmbfgsflag=_%s_lmbfgs             # file name modifier for locally aligned particle stacks\n' %(params['starfile'][:-5])
		cmd+='lmbfgsext=mrcs                 # extension for locally aligned particle stacks\n'
		cmd+='movielist=%s           # name of list of movie files\n' %(movielist)
		cmd+='coordlist=%s           # name of list of coordinate files\n'%(coordlist)
		cmd+='vecext=vec                     # extension for output vector files\n'
		cmd+='# Alignparts_lmbfgs variables that do not usually need adjustment\n'
		cmd+='maxparts=100000                # maximum number of particles that may be in a micrograph\n'
		cmd+='nsigma=5                       # standard deviations above the mean to remove outlier pixels\n'
		cmd+='rmax1=500                      # Low resolution cutoff (in Angstroms) used for alignment\n'
		cmd+='rmax2=40                       # High resolution cutoff (in Angstroms) used for alignment\n'
		cmd+='bfactor=2000                   # B-factor (in A**2) used for alignment\n'
		cmd+='factr=1d7                      # Accuracy of minimizer (eg.1.0d7=1.0x10^7 in dble precision)\n'
		cmd+='#=====================================================================\n'
		cmd+='time $alignparts_starfilehandler << eot1\n'
		cmd+='$instar\n'	
		cmd+='$outstar\n'
		cmd+='$movielist\n'
		cmd+='$coordlist\n'
		cmd+='$boxsize,$framex,$framey\n'
		cmd+='$moviepath\n'
		cmd+='$particlepath\n'
		cmd+='$movieflag\n'
		cmd+='$movieext\n'
		cmd+='$lmbfgsflag\n'
		cmd+='$lmbfgsext\n'	
		cmd+='eot1\n'
		cmd+='\n'
		
		cmd+='time $alignparts_lmbfgs << eot2\n'
		cmd+='$movielist\n'
		cmd+='$coordlist\n'
		cmd+='$boxsize,0,$particleradius,$pixelsize,$nsigma,$rmax1,$rmax2\n'
		cmd+='$expweight,$akv,$expperframe\n'
		cmd+='$bfactor,$smooth,$exaggerate,$zeroframe,$invertoutput\n'
		cmd+='$localavg,$maxparts,$localavgsigma\n'
		cmd+='$framefirstali,$framelastali,$framefirstave,$framelastave\n'
		cmd+='$factr\n'
		cmd+='$moviepath\n'
		cmd+='$particlepath\n'
		cmd+='$lmbfgsflag\n'
		cmd+='$vecext\n'
		cmd+='$lmbfgsext\n'
		cmd+='eot2\n'

		run = open('align_lmbfgs.bash','w')
        	run.write(cmd)
        	run.close()

        	cmd2 = 'chmod +x align_lmbfgs.bash'
        	subprocess.Popen(cmd2,shell=True).wait()	

		print '\nRunning LM-BFGS ...\n'

		actualnprocs=splitSTAR(params['starfile'],params['nprocs'],params['debug'])

		if params['debug'] is True:
			print 'Actual nprocs=%i' %(actualnprocs)

		nproc=1

		while nproc <= actualnprocs:

			cmd3 = './align_lmbfgs.bash %s_set%i.star %i'%(params['starfile'][:-5],nproc,nproc)
			if params['debug'] is True:
				print cmd3
			subprocess.Popen(cmd3,shell=True)
			time.sleep(5)
			nproc=nproc+1
		combineSTARfiles('%s_set' %(params['starfile'][:-5]),actualnprocs,'%s_lmbfgs.star' %(params['starfile'][:-5]))

		#check to see if completed

		time.sleep(10)

		testnum=0
		testnumcounter=1

		while testnumcounter<=500:

			time.sleep(20)
	
			testnum2=len(glob.glob('Particles/Micrographs/*.vec'))

			if testnum2 != testnum:
				testnum=testnum2
				testnumcounter=testnumcounter+1
				continue

			if testnum2 == testnum:
				time.sleep(20)
				testnum3=len(glob.glob('Particles/Micrographs/*.vec'))
				if testnum3 == testnum:
					time.sleep(60)
					testnum4=len(glob.glob('Particles/Micrographs/*.vec'))

					if testnum4 == testnum:
						print 'LM-BFGS finished.\n'

						#Clean up
						cmd='rm Micrographs/*.coord Particles/Micrographs/*.vec align_lmbfgs.bash movie.txt coord.txt %s_set*star' %(params['starfile'][:-5])
						subprocess.Popen(cmd2,shell=True).wait()			
						sys.exit()	
