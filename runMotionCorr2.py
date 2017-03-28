#!/usr/bin/env python

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
        parser.set_usage("%prog --dir=<folder with mrc frames> --gain_ref=<gain reference in mrc format with full path;input the *_norm* file from the leginon reference directory> --save_bin <save binned mic> \n\nThis program takes movies with .mrcs or .frames.mrc extensions and will create aligned movies with -a.mrc extension.\n If dose weighting, output aligned movie will have the extension -a_DW.mrc\n")
        parser.add_option("--microlist",dest="inputlist",type="string",metavar="FILE",default='empty',
                    help="Provide list of movies instead of directory location")
	parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Or provide directory containing direct detector movies (.mrcs extension unless using gain reference, then .frames.mrc)")
        parser.add_option("--gain_ref",dest="gain_ref",type="string",metavar="FILE",default='empty',
                    help="Gain reference file from Leginon with the full path (.mrc)")
	parser.add_option("--throw",dest="throw",type="int",metavar="INT",default=2,
                    help="Number of initial frames to discard from alignment. (Default=2)")
	parser.add_option("--binning",dest="binning",type="int",metavar="INT",default=2,
                    help="Scaling factor for output aligned movie (Default=2) ")
	parser.add_option("--patchsize",dest="patch",type="int",metavar="INT",default=5,
                    help="Number of patches for local alignment (Default=5, which is a 5 x 5 tiling)")
	parser.add_option("--bfactor",dest="bfactor",type="int",metavar="INT",default=100,
                    help="Bfactor to use during movie alignment. Must be positive value (Default=100)")
	parser.add_option("--dose",dest="doserate",type="float",metavar="FLOAT",default=0,
                    help="Optional: Input dose rate for dose weighting (electrons per Angstrom-squared")
	parser.add_option("--kev",dest="kev",type="int",metavar="INT",
                    help="Optional: If dose weighting, provide accelerating voltage in keV")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
                    help="Optional: If dose weighting, provide pixel size in Angstroms/pixel")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                    help="debug")

        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))
        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#=============================
def checkConflicts(params):

	if params['bfactor'] < 0: 
		print 'Error: Bfactor must be positive. Exiting'
		sys.exit()

        if params['inputlist'] == 'empty':
		if not os.path.exists(params['dir']):
                	print "\nError: input directory %s doesn't exists, exiting.\n" %(params['dir'])
                	sys.exit()
	cudapath = subprocess.Popen("env | grep cuda-7.5", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if not cudapath:
                print '\nError: cuda-7.5 not found in PATH. Please load cuda/cuda-7.5 and try again.'
		sys.exit()

#==============================
def alignmovies(params,motionCor2Path):

    if params['gain_ref'] == 'empty':
          print '\nInput movies are assumed to be already gain corrected.\n'
	  if params['debug'] is True:
	       print 'Input will be .mrcs'
          if params['inputlist'] is 'empty':
	       mrcList = sorted(glob.glob('%s/*.mrcs'%(params['dir'])))
	  if params['inputlist'] is not 'empty':
	       mrcList = open(params['inputlist'],'r')
    if params['gain_ref'] != 'empty':
	  if params['debug'] is True:
               print 'Input will be .frames.mrc'
          if params['inputlist'] is 'empty':
               mrcList = sorted(glob.glob('%s/*.frames.mrc'%(params['dir'])))
          if params['inputlist'] is not 'empty':
               mrcList = open(params['inputlist'],'r')

    for mrcs in mrcList:
        if params['debug'] is True:
            print 'Working on movie %s' %(mrcs)
	mrcs=mrcs.split()[0]
	extension=mrcs.split('.')[-1]
	inmovie=mrcs
	if extension == 'mrcs':
		outmicro='%s-a.mrc' %(mrcs[:-5])

	if extension != 'mrcs':
		outmicro='%s-a.mrc' %(mrcs[:-11])

	if params['debug'] is True:
		print 'Input movie: %s' %(inmovie)
		print 'Output aligned movie filename: %s' %(outmicro)

        if os.path.exists(outmicro):
            print "\n"
            print 'Micrograph %s already exists, skipping %s' %(outmicro,inmovie)
            print "\n"
            continue
        
	if params['doserate'] > 0:
		doseinfo='-PixSize %f -kV %i' %(params['apix'],params['kev'])
	if params['doserate'] == 0:
		doseinfo=''

	if params['gain_ref'] != 'empty':
		gainref='-Gain %s' %(params['gain_ref'])
	if params['gain_ref'] == 'empty':
		gainref=''

	if params['debug'] is True:
		print params['doserate']

	cmd = '%s -InMrc %s -OutMrc %s -Throw %i -Bft %i -Iter 10 -Patch %i %i -FtBin %i -FmDose %f %s %s' %(motionCor2Path,inmovie,outmicro,params['throw'],params['bfactor'],params['patch'],params['patch'],params['binning'],params['doserate'],doseinfo,gainref)
 
	if params['debug'] is True:
		print cmd
        subprocess.Popen(cmd,shell=True).wait()

#=============================
def getimodPath():

        imodpath = subprocess.Popen("env | grep IMOD_DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if imodpath:
                return imodpath
        print "IMOD was not found, make sure IMOD is in your path"
        sys.exit()
#=============================
def getmotionCor2Path(params):

        #first try same directory as this script
        currentPath = sys.argv[0]
	currentPath=currentPath[:-18]
	if os.path.exists('%s/MotionCor2/MotionCor2-03-16-2016' %(currentPath)):
                return '%s/MotionCor2/MotionCor2-03-16-2016' %(currentPath)

        print "\nError: path to MotionCor2-03-16-2016 does not exist, exiting.\n"
        sys.exit()
#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
        imodpath=getimodPath()
        motionCor2Path=getmotionCor2Path(params)
        alignmovies(params,motionCor2Path)

	#CLeanup 
	if os.path.exists('%s_rot-180.mrc' %(params['gain_ref'])):
		os.remove('%s_rot-180.mrc' %(params['gain_ref']))
		os.remove('%s_rot-180_flipy.mrc' %(params['gain_ref']))
