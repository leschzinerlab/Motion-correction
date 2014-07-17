#!/usr/bin/env python 

#This script will convert an imagic or spider stack into the appropriate file format for XMIPP/Relion,
#which is single particle images in a new folder
#
#TEM|pro - mcianfrocco

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
        parser.set_usage("%prog -i <stack.img> -o <output folder name> --num=[number of particles in stack]")
        parser.add_option("-i",dest="stack",type="string",metavar="FILE",
                help="Particle stack in .img or .spi format")
        parser.add_option("-o",dest="folder",type="string",metavar="FILE",
                help="Output folder name for single particles")
        parser.add_option("--num",dest="numParts",type="int", metavar="INT",
                help="Number of particles in stack")
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
        if not params['stack']:
                print "\nWarning: no stack specified\n"
        elif not os.path.exists(params['stack']):
                print "\nError: stack file '%s' does not exist\n" % params['stack']
                sys.exit()
        if params['stack'][-4:] != '.img':
                if params['stack'][-4:] != '.spi':
                        if params['stack'][-4:] != '.hed':
                                print 'Stack extension %s is not recognized as .spi, .hed or .img file' %(params['stack'][-4:])
                                sys.exit()

        if os.path.exists(params['folder']):
                print "\nError: output folder already exists, exiting.\n"
                sys.exit()

#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)
