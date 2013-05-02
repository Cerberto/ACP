#!/bin/bash

Nmax=100

for ((N=10; N<=Nmax; L+=10))
do
	./GP_eq 1.0e-06 $N 0.8
done

for ((N=10; N<=Nmax; L+=10))
do
	./GP_eq 1.0e-06 $N -0.8
done


exit
