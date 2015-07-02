#  SB: original script 2013-11
# JLR: create a file called shiftlist.txt containing all of the 
#      output files from alignframes_lmbfgs
#Script to display all shifts from full frame alignment using Gnuplot
#!/bin/csh -x
rm -f Gnuplot.com
clear

foreach line (`cat shiftlist.txt`)
echo $line

set shifttxt = ${line}

#echo "set yrange [0:1.5]" > Gnuplot.com
echo "set xrange [1:]" > Gnuplot.com

echo 'shift_file ="'$shifttxt'"' >> Gnuplot.com

echo 't1 ="x-shift"' >> Gnuplot.com
echo 't2 = "y-shift"' >> Gnuplot.com


echo set title ' " '$line' " ' >> Gnuplot.com
echo set xlabel ' "Frame Number" ' >> Gnuplot.com
echo set ylabel ' "Shift" ' >> Gnuplot.com
echo plot shift_file using 1:2 with lines lc rgb \'black\' title t1, \\ >> Gnuplot.com
echo      shift_file using 1:3 with lines lc rgb \'red\' title t2  >> Gnuplot.com
echo pause -1 '"'Hit any key to continue'"' >> Gnuplot.com

gnuplot Gnuplot.com

rm -f Gnuplot.com

end

