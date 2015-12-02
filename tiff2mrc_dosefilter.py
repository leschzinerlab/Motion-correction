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

#This script runs unblur on a large number of movies in .mrc format
#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog --dir=<folder with tif frames> --gain_normalize <use only if gain normalization required> --gain_ref=<gain reference needed only for gain normalization>")
        parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Directory containing direct detector movies with .tiff extension")
        parser.add_option("--gain_normalize", action="store_true",dest="gain_normalize",default=False,
                    help="Do gain normalization")
	parser.add_option("--gain_ref",dest="gain_ref",type="string",metavar="FILE",
                    help=" .dm4 gain reference file name with the full path")
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
	if params['gain_normalize'] is True:
		if not os.path.exists(params['gain_ref']):
			print "\nError: input gain reference file %s dosen't exist, exiting.\n" %(params['gain_ref'])
			sys.exit()

#==============================
def convtif2mrc(params):


    if params['gain_normalize'] is True:
    	cmd = 'dm2mrc %s %s/gatan_gain_Ref.mrc' %(params['gain_ref'], params['dir'])
        print cmd
        subprocess.Popen(cmd,shell=True).wait()


    tifList = sorted(glob.glob('%s/*.tif'%(params['dir'])))

    for mrcs in tifList:
        if params['debug'] is True:
            print mrcs

        if os.path.exists('%s.mrc' %(mrcs[:-4])):
	    print "\n"
            print ' Normalized micrograph %s.mrc already exists, skipping %s' %(mrcs[:-4],mrcs)
	    print "\n"
            continue

	if os.path.exists('%s_normalized.mrc' %(mrcs[:-4])):
            print "\n"
            print ' Normalized micrograph %s_normalized.mrc already exists, skipping %s' %(mrcs[:-4], mrcs)
            print "\n"
            continue

	print 'Working on %s' %(mrcs)
	cmd = 'tif2mrc %s %s.mrc' %(mrcs, mrcs[:-4])
	print cmd
        subprocess.Popen(cmd,shell=True).wait()
	

	if params['gain_normalize'] is True:
		cmd = 'clip mult -n 16 %s.mrc %s/gatan_gain_Ref.mrc %s_normalized.mrc' %(mrcs[:-4],params['dir'],mrcs[:-4])
		print cmd
        	subprocess.Popen(cmd,shell=True).wait()
		cmd = 'rm -rf %s.mrc' %(mrcs[:-4])
		print cmd
                subprocess.Popen(cmd,shell=True).wait()
		
#=============================
def getimodPath():

	imodpath = subprocess.Popen("env | grep IMOD_DIR", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if imodpath:
                return imodpath
        print "IMOD was not found, make sure IMOD is in your path"
        sys.exit()

#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
	imodpath=getimodPath()
	convtif2mrc(params)
