#!/usr/bin/env python

import glob 
import shutil

mrclist = glob.glob('*.mrc')

for mrc in mrclist:

	print 'Moving %s to %s.mrcs' %(mrc,mrc[:-4])
	shutil.move(mrc,'%s.mrcs'%(mrc[:-4]))	
