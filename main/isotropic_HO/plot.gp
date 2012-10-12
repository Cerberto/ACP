
reset

#set terminal epslatex size 12.5cm,7.5cm color colortext
set terminal postscript eps enhanced color font 'Helvetica,10'

#set logscale x 2
#set log x
set grid x
set grid y
set format '$%g$'
set ylabel "$y(E,x_{max})$"
set xlabel "$E$"

#unset colorbox

#
# define line styles using explicit rgbcolor names
#
set style line 1 lt 2 lc rgb "red" lw 3
#set style line 2 lt 2 lc rgb "green" lw 2
#set style line 3 lt 2 lc rgb "orange" lw 2


unset key
#set output 'plot_yRmaxE.tex'
set output 'plot_yRmaxE.eps'

#show style line

plot 'yRmaxE_0.dat' with linespoints ls 1
