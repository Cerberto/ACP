#!/usr/bin/perl -w

for($alpha = 0; $alpha <= 0.5; $alpha+=0.05){
	system("./GP_eq","1.0e-4","1000", "$alpha");
}

printf "\nFine!\n\n"
