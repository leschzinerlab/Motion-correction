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
        parser.set_usage("%prog --dir=<folder with mrc frames> --gain_ref=<gain reference in mrc format with full path> ")
        parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Directory containing direct detector movies with .mrc extension")
        parser.add_option("--gain_ref",dest="gain_ref",type="string",metavar="FILE",
                    help=" .mrc gain reference file name with the full path")
        parser.add_option("--outext",dest="outext",type="string",metavar="STRING",default='-a.mrcs',
                    help="Output extension for normalized movied (Default=-a.mrcs)")
        parser.add_option("--bin",dest="bin",type="int",metavar="INTEGER",default=1,
                    help="Binning factor for output movie (Default=1)")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                    help="debug")

        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))
        if len(sys.argv) <= 3:
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
def normalize(params):

    if not os.path.exists('%s_rot-180.mrc'%(params['gain_ref'][:-4])):
    	cmd = 'newstack -rot -180 %s %s_rot-180.mrc' %(params['gain_ref'], params['gain_ref'][:-4])
    	print cmd
    	subprocess.Popen(cmd,shell=True).wait()

    if not os.path.exists('%s_rot-180_flipy.mrc'%(params['gain_ref'][:-4])):
    	cmd = 'clip flipy %s_rot-180.mrc %s_rot-180_flipy.mrc' %(params['gain_ref'][:-4], params['gain_ref'][:-4])
    	print cmd
    	subprocess.Popen(cmd,shell=True).wait()

    if not os.path.exists('%s_rot-180_flipy_bin.mrc'%(params['gain_ref'][:-4])):
        cmd = 'proc2d %s_rot-180_flipy.mrc %s_rot-180_flipy_bin.mrc meanshrink=%i' %(params['gain_ref'][:-4], params['gain_ref'][:-4],params['bin'])
        print cmd
        subprocess.Popen(cmd,shell=True).wait()

    mrcList = sorted(glob.glob('%s/*.mrc'%(params['dir'])))

    for mrcs in mrcList:
        if params['debug'] is True:
            print mrcs

        inputfile = '%s' %(mrcs)
        outputfile = '%s%s' %(mrcs[:-11],params['outext'])

	if os.path.exists(outputfile):
		continue

        if inputfile[-9:]=='flipy.mrc':
            continue
        if inputfile[-11:]=='rot-180.mrc':
            continue
       	if inputfile[-22:] == '_rot-180_flipy_bin.mrc':
	    continue
	if inputfile == params['gain_ref']:
	    continue

        print 'Working on %s' %(mrcs)
        
	if os.path.exists('tmp.mrc'):
		os.remove('tmp.mrc')

	cmd='e2proc2d.py %s tmp.mrcs --meanshrink %i' %(inputfile,params['bin'])
	print cmd
        subprocess.Popen(cmd,shell=True).wait()

        cmd = 'clip mult -n 16 tmp.mrcs %s_rot-180_flipy_bin.mrc %s%s' %(params['gain_ref'][:-4],mrcs[:-11],params['outext'])
        if params['debug'] is True:
		print cmd
        subprocess.Popen(cmd,shell=True).wait()

	if os.path.exists('tmp.mrcs'):
		os.remove('tmp.mrcs')

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

        if os.path.exists('%s/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe' %(currentPath)):
                return '%s/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe' %(currentPath)

        if os.path.exists('unblur_1.0.2/bin/unblur_openmp_7_17_15.exe'):
                return 'unblur_1.0.2/bin/unblur_openmp_7_17_15.exe'

        if os.path.exists('/home/indrajit/Desktop/Scripts/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe'):
                return '/home/indrajit/Desktop/Scripts/unblur_1.0.2/bin/unblur_openmp_7_17_15.exe'


        print "\nError: path to unblur not found, exiting.\n"
        sys.exit()
#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
        imodpath=getimodPath()
        normalize(params)
