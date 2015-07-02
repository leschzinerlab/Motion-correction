# card 1: ftsize (e.g. 4096, bigger than movie size), nsigma (e.g. 5; removes outliers)
# card 2: input movie list (e.g. movielist.txt)
# card 3: input shift lists (e.g. shiftlist.txt; file order must correspond with input movie list)
# card 4: input path (e.g. './'; location of input movie and shift files
# card 5: output path (e.g. './'; location output files will be placed
# card 6: rootname modifying flag (e.g. ' "" ' inorder to _not_ modify rootname; otherwise '_aligned')
# card 7: aliflag(1/0; turn on or off output of aligned movies),aliext (extension for aligned movies. e.g. 'mrcs')
# card 8: avgflag(1/0; turn on or off output of aligned movies),avgext (extension for output average, e.g. 'mrc')
# card 9: framefirst,framelast (set frame last to 0 to use all frames; e.g. "1,30", or "5,25", or "5,0")
time /home/jlr/programs/shiftframes/shiftframes_list.exe <<eot
4096,5
movielist.txt
shiftlist.txt
./
./
"_aligned"
0,mrcs
1,mrc
1,0
eot
