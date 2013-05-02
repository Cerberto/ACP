#!/usr/bin/perl -w

###############################################################
#
# Script per binning di dati
#
# eseguire con:
# ./gaussian_binning.pl input.dat > output.dat
#
# input.dat ha una sola colonna, contenente i dati,
#
# output.dat ha 3 colonne che sono:
#	_ estremo inferiore del bin;
#	_ estremo superiore del bin;
#	_ numero di dati nel bin.
#
###############################################################


$step = 100;
$x_min = -4;
$x_max = 4;

$dx = ($x_max - $x_min)/$step;

@HISTO = (0) x $step;

$i = 0;

while(<>){
	chomp;
	$i++;
	$j = int (($_ - $x_min)/$dx);
	$HISTO[$j]++;
}

for($i = 0; $i < $step; $i++){
	$x1 = $x_min + $i*$dx;
	$x2 = $x_min + ($i + 1)*$dx;
	printf "%e\t%e\t%d\n", $x1, $x2, $HISTO[$i];
}
