
reset

#set terminal epslatex size 12.5cm,7.5cm color colortext
set terminal postscript eps enhanced color font 'Helvetica,10'

#set logscale x 2
#set log x
set grid x
set grid y
#set format '$%g$'
#set ylabel "$y(E,x_{max})$"
#set xlabel "$E$"

#unset colorbox

#
# define line styles using explicit rgbcolor names
#
set style line 1 lt 1 lc rgb "red" lw 3
set style line 4 lt 2 lc rgb "orange" lw 1
set style line 2 lt 1 lc rgb "green" lw 2
set style line 3 lt 2 lc rgb "green" lw 1

unset key


######################################
# PLOT EIGENVALUES
######################################
#set xrange[9.9:23.1]
set output 'ev_vs_iteration.eps'
plot 'eigenvalues.dat' with linespoints ls 1


######################################
# PLOT DENSITY
######################################
set output 'density_effpot.eps'
set xrange [.01:5]

alpha = 0.5
scale = 10
npar = 10
f(x) = npar*exp(-x*x)/(pi**1.5)

plot 'ultimate.dat' u 1:(0.5*$1*$1 + 4*pi*alpha*$3) with lines ls 2, \
	'ultimate.dat' u 1:(scale*($3)) with lines ls 1, \
	0.5*x*x with lines ls 3, \
	'ultimate.dat' u 1:(scale*f($1)) w l ls 4
