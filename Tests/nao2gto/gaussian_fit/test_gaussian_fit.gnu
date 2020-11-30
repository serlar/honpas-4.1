#!/usr/bin/env gnuplot

set style line 1 lt 7 lw 2 pt 2

plot "test_gaussian_fit.out" using 1:2 with linespoints title "Data", \
     "test_gaussian_fit.out" using 1:3 with lines title "G1 fit", \
     "test_gaussian_fit.out" using 1:4 with lines title "G2 fit", \
     "test_gaussian_fit.out" using 1:5 with lines title "G3 fit", \
     "test_gaussian_fit.out" using 1:6 with lines title "G4 fit", \
     "test_gaussian_fit.out" using 1:7 with lines title "G5 fit", \
     "test_gaussian_fit.out" using 1:8 with lines title "G6 fit", \
     "test_gaussian_fit.out" using 1:9 with lines title "G7 fit", \
     "test_gaussian_fit.out" using 1:10 with lines title "G8 fit"
pause -1
plot "test_gaussian_fit.chk" using 1:2 with lines title "Residuals"
pause -1
