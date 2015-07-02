# card 1: inlist (list of files containing movies to be aligned)
# card 2: boxsie (ftsize, smaller than movie), psize (e.g. 1.45 A) ,nsigma (e.g. 5) ,rmax1 (500A),rmax2 (100A) 
# card 3: bfactor (e.g. 2000 A**2)
# card 4: framefirst (e.g. 1), framelast(0=nframes), zeroframe (e.g. 15)
# card 5: factr (1d1 [v. high precision=not necessary], 1d7 [med=good], 1d14 [low=bad])
# card 6: inpath (path for input movies)
# card 7: outpath (path for output shift files)
# card 8: vecext (extension for output shift files)
time /home/jlr/programs/k2_lmbfgs/alignframes_lmbfgs.exe << eot
movielist.txt
3600,1.45,5,500,100
2000
1,0,15
1d7
./
./
shf
eot
