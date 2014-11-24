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
        parser.set_usage("%prog --dir=<folder with micrographs>")
        parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                help="Directory containing direct detector movies with .mrcs extension")
        parser.add_option("--save", action="store_true",dest="save",default=False,
                help="Saved aligned, binned image for easy viewing")
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

        if not os.path.exists(params['dir']):
                print "\nError: directory %s doesn't exists, exiting.\n" %(params['dir'])
                sys.exit()

	#Check if CUDA is loaded
	cudapath = subprocess.Popen("env | grep cuda", shell=True, stdout=subprocess.PIPE).stdout.read().strip()	
	if not cudapath:
		print "\nError: CUDA libraries are not loaded. Load CUDA and try again.\n"
		sys.exit()
#==============================
def alignDDmovies(params,motionCorrPath):

	mrcsList = sorted(glob.glob('%s/*.mrcs'%(params['dir'])))
	currentPath = os.getcwd()

	for mrcs in mrcsList:
		
		if params['debug'] is True:
			print mrcs

		if os.path.exists('%s.mrc'%(mrcs[:-5])):
			continue

		print 'Motion correcting movie %s --> %s.mrc' %(mrcs,mrcs[:-5])

		cmd = '%s %s -fcs %s.mrc' %(motionCorrPath,mrcs,mrcs[:-5]) 
		subprocess.Popen(cmd,shell=True).wait()
	
		mrcsOnly = mrcs.split('/')

		if params['save'] is True:
			print 'Saving binned image'
			shutil.move('%s/dosef_quick/%s_CorrSum.mrc'%(currentPath,mrcsOnly[-1][:-5]),'%s_CorrSum.mrc'%(mrcs[:-5]))
	
		#Clean up	
		os.remove('%s_Log.txt' %(mrcs[:-5]))
		shutil.rmtree('dosef_quick')	
		
#=============================
def getMotionCorrPath(params):

	#first try same directory as this script
	currentPath = sys.argv[0]
	currentPath = currentPath[:-30]

	if os.path.exists('%s/motioncorr' %(currentPath)):
		return '%s/motioncorr' %(currentPath)

	if os.path.exists('motioncorr'):
		return 'motioncorr'

	if os.pathexists('/usr/local/bin/motioncorr'):
		return '/usr/local/bin/motioncorr'
		
#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
	motionCorrPath=getMotionCorrPath(params)
	alignDDmovies(params,motionCorrPath)
