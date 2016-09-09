#!/usr/bin/env python

import glob 
import os

folder='/data/microscopes/frames/16aug11c/rawdata/'
extension='.frames.mrc'
outfilename='micrographs_'

miclist=sorted(glob.glob('%s/*%s' %(folder,extension)))

numgroups=4
nummics=924
micpergroup=232

totmiccounter=1
currentgroup=1
groupmiccounter=1

for mic in miclist:

	if not os.path.exists('%s%i.txt' %(outfilename,currentgroup)):
		o1=open('%s%i.txt' %(outfilename,currentgroup),'w')

	o1.write('%s\n' %(mic))
		

	groupmiccounter=groupmiccounter+1
	totmiccounter=totmiccounter+1

	if groupmiccounter > micpergroup:
		currentgroup=currentgroup+1
		groupmiccounter=1
