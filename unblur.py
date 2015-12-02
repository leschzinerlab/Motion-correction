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
        parser.set_usage("%prog --dir=<folder with micrographs> --frame=<number of frames per movie> --apix=<pixel sixe> --dose_filt <apply dose filter> --exp=<exposure per frame> --volt=<accelerating voltage> --preexp=<preexposure amount>")
        parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Directory containing direct detector movies with .mrcs extension")
        parser.add_option("--frame",dest="frame",type="int",metavar="INTEGER",default=50,
                    help="number of frames per movie (Default = 50)")
        parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",
                    help="pixel size at the object level.")
        parser.add_option("--dose_filt", action="store_true",dest="dose_filt",default=False,
                    help="Apply dose dilter")
        parser.add_option("--exp",dest="exp",type="float",metavar="FLOAT",
                    help="exposure per frame.")
        parser.add_option("--volt",dest="volt",type="int",metavar="INTEGER",default=300,
                    help="accelerating voltage (Default = 300)")
        parser.add_option("--preexp",dest="preexp",type="float",metavar="FLOAT",default=0,
                    help="preexposure amount in e/A^2 (Default = 0).")
        parser.add_option("--save", action="store_true",dest="save",default=False,
                    help="Save aligned frames")
        parser.add_option("--expert", action="store_true",dest="expert",default=False,
                    help="Set expert options")
        parser.add_option("--minshift",dest="minshift",type="int",metavar="INTEGER",default=2,
                    help="expert option; minimum shift for initial search in A (Default = 0).")
        parser.add_option("--outrad",dest="outrad",type="int",metavar="INTEGER",default=200,
                    help="expert option; outer raduis shift limit in A (Default = 200).")
        parser.add_option("--Bfact",dest="Bfact",type="int",metavar="INTEGER",default=1500,
                    help="expert option; B factor to apply in A^2 (Default = 1500).")
        parser.add_option("--hw_vertical",dest="hw_vertical",type="int",metavar="INTEGER",default=1,
                    help="expert option; half-width of central vertical line of fourier mask (Default = 1).")
        parser.add_option("--hw_horizontal",dest="hw_horizontal",type="int",metavar="INTEGER",default=1,
                    help="expert option; half-width of central horizontal line of fourier mask (Default = 1).")
        parser.add_option("--termination",dest="termination",type="float",metavar="FLOAT",default=0.1,
                    help="expert option; termination shift threshold (Default = 0.1).")
        parser.add_option("--iter",dest="iter",type="int",metavar="INTEGER",default=10,
                    help="expert option; maximum number of iterations (Default = 10).")
        parser.add_option("--no_restore", action="store_true",dest="no_restore",default=False,
                    help="expert option; do not restore noise power (Default is restore)")
        parser.add_option("--verbose", action="store_true",dest="verbose",default=False,
                    help="expert option; verbose output (Default is no)")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                    help="debug")

        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 3:
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

#==============================
def alignDDmovies(params,unblurPath):

    mrcsList = sorted(glob.glob('%s/*.mrc'%(params['dir'])))
    currentPath = os.getcwd()

    for mrcs in mrcsList:
        if params['debug'] is True:
            print mrcs

        if os.path.exists('%s_aligned.mrc' %(mrcs[:-4])):
	    print "\n"
            print 'Micrograph %s.mrc already exists, skipping %s' %(mrcs[:-5],mrcs)
	    print "\n"
            continue

	    print 'Motion correcting movie %s --> %s_aligned.mrc' %(mrcs,mrcs[:-4])

        #Marker for dose filter
        dose='no'
        if params['dose_filt'] is True:
            dose='yes'

	#Marker for restore noise power. Note that default option is restore
	restore='yes'
	if params['no_restore'] is True:
	    restore='no'
	
	#Marker for verbose
	verb='no'
        if params['verbose'] is True:
            verb='yes'

       #Saving aligned frames
        if params['save'] is True:
            #Not modifying expert options
            if params['expert'] is False:
                shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
                shell=shell.split('/')[-1][:-1]

                do_dosefilt_save='#!/bin/%s -x\n' %(shell)
                do_dosefilt_save+='%s << eof\n' %(unblurPath)
                do_dosefilt_save+='%s\n' %(mrcs)
                do_dosefilt_save+='%s\n' %(params['frame'])
                do_dosefilt_save+='%s_aligned.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='%s_shift.txt\n' %(mrcs[:-4])
                do_dosefilt_save+='%s\n' %(params['apix'])
                do_dosefilt_save+='%s\n' %(dose)
                do_dosefilt_save+='%s\n' %(params['exp'])
                do_dosefilt_save+='%s\n' %(params['volt'])
                do_dosefilt_save+='%s\n' %(params['preexp'])
                do_dosefilt_save+='yes\n'
                do_dosefilt_save+='%s_aligned_frames.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='no\n'
                do_dosefilt_save+='eof\n'

	    #Modifying expert options	
	    if params['expert'] is True:
                shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
                shell=shell.split('/')[-1][:-1]

                do_dosefilt_save='#!/bin/%s -x\n' %(shell)
                do_dosefilt_save+='%s << eof\n' %(unblurPath)
                do_dosefilt_save+='%s\n' %(mrcs)
                do_dosefilt_save+='%s\n' %(params['frame'])
                do_dosefilt_save+='%s_aligned.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='%s_shift.txt\n' %(mrcs[:-4])
                do_dosefilt_save+='%s\n' %(params['apix'])
                do_dosefilt_save+='%s\n' %(dose)
                do_dosefilt_save+='%s\n' %(params['exp'])
                do_dosefilt_save+='%s\n' %(params['volt'])
                do_dosefilt_save+='%s\n' %(params['preexp'])
                do_dosefilt_save+='yes\n'
                do_dosefilt_save+='%s_aligned_frames.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='yes\n'
		do_dosefilt_save+='%s_frc.txt\n' %(mrcs[:-4])
		do_dosefilt_save+='%s\n' %(params['minshift'])
		do_dosefilt_save+='%s\n' %(params['outrad'])
		do_dosefilt_save+='%s\n' %(params['Bfact'])
		do_dosefilt_save+='%s\n' %(params['hw_vertical'])
		do_dosefilt_save+='%s\n' %(params['hw_horizontal'])
		do_dosefilt_save+='%s\n' %(params['termination'])
		do_dosefilt_save+='%s\n' %(params['iter'])
		do_dosefilt_save+='%s\n' %(restore)
		do_dosefilt_save+='%s\n' %(verb)
                do_dosefilt_save+='eof\n'

        #Not saving aligned frames
        if params['save'] is False:
            #Not modifying expert options
            if params['expert'] is False:
                shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
                shell=shell.split('/')[-1][:-1]

                do_dosefilt_save='#!/bin/%s -x\n' %(shell)
                do_dosefilt_save+='%s << eof\n' %(unblurPath)
                do_dosefilt_save+='%s\n' %(mrcs)
                do_dosefilt_save+='%s\n' %(params['frame'])
                do_dosefilt_save+='%s_aligned.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='%s_shift.txt\n' %(mrcs[:-4])
                do_dosefilt_save+='%s\n' %(params['apix'])
                do_dosefilt_save+='%s\n' %(dose)
                do_dosefilt_save+='%s\n' %(params['exp'])
                do_dosefilt_save+='%s\n' %(params['volt'])
                do_dosefilt_save+='%s\n' %(params['preexp'])
                do_dosefilt_save+='no\n'
                do_dosefilt_save+='no\n'
                do_dosefilt_save+='eof\n'

	    #Modifying expert options
            if params['expert'] is True:
                shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
                shell=shell.split('/')[-1][:-1]

                do_dosefilt_save='#!/bin/%s -x\n' %(shell)
                do_dosefilt_save+='%s << eof\n' %(unblurPath)
                do_dosefilt_save+='%s\n' %(mrcs)
                do_dosefilt_save+='%s\n' %(params['frame'])
                do_dosefilt_save+='%s_aligned.mrc\n' %(mrcs[:-4])
                do_dosefilt_save+='%s_shift.txt\n' %(mrcs[:-4])
                do_dosefilt_save+='%s\n' %(params['apix'])
                do_dosefilt_save+='%s\n' %(dose)
                do_dosefilt_save+='%s\n' %(params['exp'])
                do_dosefilt_save+='%s\n' %(params['volt'])
                do_dosefilt_save+='%s\n' %(params['preexp'])
                do_dosefilt_save+='no\n'
		do_dosefilt_save+='yes\n'
                do_dosefilt_save+='%s_frc.txt\n' %(mrcs[:-4])
                do_dosefilt_save+='%s\n' %(params['minshift'])
                do_dosefilt_save+='%s\n' %(params['outrad'])
                do_dosefilt_save+='%s\n' %(params['Bfact'])
                do_dosefilt_save+='%s\n' %(params['hw_vertical'])
                do_dosefilt_save+='%s\n' %(params['hw_horizontal'])
                do_dosefilt_save+='%s\n' %(params['termination'])
                do_dosefilt_save+='%s\n' %(params['iter'])
                do_dosefilt_save+='%s\n' %(restore)
                do_dosefilt_save+='%s\n' %(verb)
                do_dosefilt_save+='eof\n'
               


        if os.path.exists('do_dosefilt_save.com'):
		    os.remove('do_dosefilt_save.com')

        do_dosefilt_saveFile = open('do_dosefilt_save.com','w')
        do_dosefilt_saveFile.write(do_dosefilt_save)
        do_dosefilt_saveFile.close()

        cmd = 'chmod +x do_dosefilt_save.com'
	subprocess.Popen(cmd,shell=True).wait()

	cmd = './do_dosefilt_save.com'
	subprocess.Popen(cmd,shell=True).wait()

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
	unblurPath=getunblurPath(params)
	alignDDmovies(params,unblurPath)
