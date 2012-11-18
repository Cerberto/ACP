#! /bin/bash

ls *.data | sed "s/.data//" > list

for i in `cat list` ; do
	sed -e "s/INPUTFILE/$i/" -e "s/OUTPUTFILE/$i/" \
	plot.gnu | gnuplot
done

rm list
