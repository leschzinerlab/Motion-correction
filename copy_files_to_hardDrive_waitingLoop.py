#!/usr/bin/env python

import glob
import os
import time
import datetime
import shutil

extDrivePath = '/media/SCOPEDRIVE2'

filesComp = sorted(glob.glob("*.mrc"))

i = 1

while i < 100:

	for file in filesComp: 

		if os.path.exists('%s/%s' %(extDrivePath,file)):
			print 'File %s already copied' %(file)
			continue
		if not os.path.exists('%s/%s' %(extDrivePath,file)): 
			print 'Copying %s to %s/%s' %(file,extDrivePath,file)
			shutil.copy('%s'%(file),'%s/%s' %(extDrivePath,file))

	ts = time.time()
	timeStamp = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

	
	print "Finished current folder at %s. Waiting 10 min." %(timeStamp)
	time.sleep(600)
	i = i + 1
