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

#This script runs IMOD and motioncorrection on a large number of movies in .mrc format
#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog --dir=<folder with mrc frames> --gain_ref=<gain reference in mrc format with full path;input the *_norm* file from the leginon reference directory> --save_bin <save binned mic> --save_norm <save normalized frames>\n\nThis program takes movies with .frames.mrcs extensions and creates aligned movies with .mrc extension, along with the option to create normalized movies with the .mrcs extension.")
        parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Directory containing direct detector movies (.frames.mrcs extension)")
        parser.add_option("--gain_ref",dest="gain_ref",type="string",metavar="FILE",
                    help="Gain reference file from Leginon with the full path (.mrc)")
        parser.add_option("--save_bin", action="store_true",dest="save_bin",default=False,
                    help="Save binned image for quick inspection")
        parser.add_option("--save_norm", action="store_true",dest="save_norm",default=False,
                    help="Save normalized movie frames as .mrcs")
        parser.add_option("--bin",dest="bin",type="int",metavar="INT",default=1,
                    help="Binning factor to use during movie alignment, 1 or 2. (Default=1)")
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

        if not os.path.exists(params['dir']):
                print "\nError: input directory %s doesn't exists, exiting.\n" %(params['dir'])
                sys.exit()
        if not os.path.exists(params['gain_ref']):
                print "\nError: input gain reference file %s dosen't exist, exiting.\n" %(params['gain_ref'])
                sys.exit()

#==============================
def normalize(params,unblurPath):

    if not os.path.exists('%s_rot-180.mrc'%(params['gain_ref'][:-4])):
    	cmd = 'newstack -rot -180 %s %s_rot-180.mrc' %(params['gain_ref'], params['gain_ref'][:-4])
    	print cmd
    	subprocess.Popen(cmd,shell=True).wait()

    if not os.path.exists('%s_rot-180_flipy.mrc'%(params['gain_ref'][:-4])):
    	cmd = 'clip flipy %s_rot-180.mrc %s_rot-180_flipy.mrc' %(params['gain_ref'][:-4], params['gain_ref'][:-4])
    	print cmd
    	subprocess.Popen(cmd,shell=True).wait()

    mrcList = sorted(glob.glob('%s/*.frames.mrcs'%(params['dir'])))

    if len(mrcList) == 0: 
	print 'Error: No .frames.mrcs movies found in directory %s' %(params['dir'])
	sys.exit()

    for mrcs in mrcList:
        if params['debug'] is True:
            print 'Working on movie %s' %(mrcs)
        

        filename = '%s' %(mrcs)
        split_file = filename.split('_')
        split_line = filename.split('/')
        actual_file_name = '%s' %(split_line[-1])
        name_a = '%s' %(split_file[-1])

        if name_a=='normalized.mrc':
            continue 
        if name_a=='SumCorr.mrc':
            continue
        if name_a=='CorrSum.mrc':
            continue
        if name_a=='1.mrc':
            continue
        if name_a=='0.mrc':
            continue
        if name_a=='flipy.mrc':
            continue
        if name_a=='rot-180.mrc':
            continue
        
	inmovie=mrcs
	normmovie='%s.mrcs' %(mrcs[:-12])
	norm='%s.mrc' %(mrcs[:-12])

	if params['debug'] is True:
		print 'Input movie: %s' %(inmovie)
		print 'Normalized output movie filename: %s' %(normmovie)
		print 'Output aligned movie filename: %s' %(norm)

        if os.path.exists(norm):
            print "\n"
            print 'Micrograph %s already exists, skipping %s' %(norm,inmovie)
            print "\n"
            continue
        
	cmd = '%s %s -fcs %s -hgr 0 -fgr %s_rot-180_flipy.mrc -bin %i' %(unblurPath,inmovie,norm,params['gain_ref'][:-4],params['bin'])
	if params['debug'] is True:
		print cmd
        subprocess.Popen(cmd,shell=True).wait()

        cmd = 'mv %s.mrc %s/' %(actual_file_name[:-12],params['dir'])
        if params['debug'] is True:
		print cmd
        subprocess.Popen(cmd,shell=True).wait()

        if params['save_bin'] is True:
            cmd = 'mv dosef_quick/%s_normalized_CorrSum.mrc %s/%s_dosefQuick.mrc' %(actual_file_name[:-12],params['dir'],actual_file_name[:-12])
            if params['debug'] is True:
		print cmd
            subprocess.Popen(cmd,shell=True).wait()

        if params['save_norm'] is True:  
            cmd = 'clip mult -n 16 %s %s_rot-180_flipy.mrc %s' %(mrcs,params['gain_ref'][:-4],normmovie)
	    if params['debug'] is True:
                print cmd
            subprocess.Popen(cmd,shell=True).wait()

    cmd = 'rm -rf dosef_quick %s_Log.txt' %(actual_file_name[:-12])
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
def getunblurPath(params):

        #first try same directory as this script
        currentPath = sys.argv[0]
        currentPath = currentPath[:-23]

        if os.path.exists('%s/motioncorr_v2.1/bin/dosefgpu_driftcorr' %(currentPath)):
                return '%s/motioncorr_v2.1/bin/dosefgpu_driftcorr' %(currentPath)

        if os.path.exists('motioncorr_v2.1/bin/dosefgpu_driftcorr'):
                return 'motioncorr_v2.1/bin/dosefgpu_driftcorr'

        if os.path.exists('/home/indrajit/Desktop/Scripts/motioncorr_v2.1/bin/dosefgpu_driftcorr'):
                return '/home/indrajit/Desktop/Scripts/motioncorr_v2.1/bin/dosefgpu_driftcorr'


        print "\nError: path to motioncorr not found, exiting.\n"
        sys.exit()
#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
        imodpath=getimodPath()
        unblurPath=getunblurPath(params)
        normalize(params,unblurPath)

	#CLeanup 
	if os.path.exists('%s_rot-180.mrc' %(params['gain_ref'])):
		os.remove('%s_rot-180.mrc' %(params['gain_ref']))
		os.remove('%s_rot-180_flipy.mrc' %(params['gain_ref']))
