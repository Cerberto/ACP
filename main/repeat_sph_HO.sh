#!/bin/bash

Lmax=5

for ((L=0; L<=Lmax; L+=1))
do
	./sph_HO $L
done

	echo ""

exit
