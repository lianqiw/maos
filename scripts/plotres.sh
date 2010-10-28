#!/usr/bin/env gnuplot
#set commands has to be before plot
set term x11
set grid xtics ytics 
set xlabel "Time Steps"
set ylabel "Wavefront Error"

plot \
    "Res_1.txt" using 1:3 with lines lt 2 lw 1 title "OL TT" ,\
    "Res_1.txt" using 1:4 with lines lt 3 lw 1 title "OL PTTR" ,\
    "Res_1.txt" using 1:6 with lines lt 2 lw 2 title "CL TT" ,\
    "Res_1.txt" using 1:7 with lines lt 3 lw 2 title "CL PTTR",\
    "Res_1.txt" using 1:8 with lines lt 4 lw 1 title "CL LGS" ,\
    "Res_1.txt" using 1:10 with lines lt 5 lw 2 title "CL NGS"
pause -1 "Hit return to continue"


set terminal postscript eps enhanced color size 6,4 solid linewidth 2 'Arial' 24
set output "Res.eps"

replot