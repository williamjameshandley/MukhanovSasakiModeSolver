set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb"white" behind
set logscale xy
plot "output/PPS.txt" u 1:2 w l
replot "output/PPS.txt" u 1:3 w l