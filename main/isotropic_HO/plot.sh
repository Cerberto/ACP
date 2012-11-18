#!bin/bash

echo "call \"gnuscript.txt\" \"'.'\" \"'data.dat'\"" | gnuplot

PATH = $0
DATAFILE = $1

plot PATH."/".DATAFILE

