#!/usr/bin/env gnuplot
#set commands has to be before plot
set term x11
set grid xtics ytics 
set xlabel "Time Steps"
set ylabel "Wavefront Error"

plot \
    "Resupt.txt" using 1:2 with lines lt 2 lw 1 title "1" ,\
    "Resupt.txt" using 1:3 with lines lt 2 lw 2 title "2" ,\
    "Resupt.txt" using 1:4 with lines lt 3 lw 1 title "3" ,\
    "Resupt.txt" using 1:5 with lines lt 3 lw 2 title "4",\
    "Resupt.txt" using 1:6 with lines lt 3 lw 1 title "5" ,\
    "Resupt.txt" using 1:7 with lines lt 3 lw 2 title "6",\
    "Resupt.txt" using 1:8 with lines lt 3 lw 1 title "7" ,\
    "Resupt.txt" using 1:9 with lines lt 3 lw 2 title "8",\
    "Resupt.txt" using 1:10 with lines lt 3 lw 1 title "9" ,\
    "Resupt.txt" using 1:11 with lines lt 3 lw 2 title "10",\
    "Resupt.txt" using 1:12 with lines lt 3 lw 1 title "11" ,\
    "Resupt.txt" using 1:13 with lines lt 3 lw 2 title "12" 

pause -1 "Hit return to continue"


set terminal postscript eps enhanced color size 6,4 solid linewidth 2 'Arial' 24
set output "Resupt.eps"

replot