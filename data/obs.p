set terminal postscript eps enhanced color dashed lw 1 "Helvetica" 16

set output "obs.eps"

plot "trace.txt" using 0:1 title ""
plot "trace.txt" using 0:2 title ""
plot "trace.txt" using 0:3 title ""
plot "trace.txt" using 0:4 title ""
plot "trace.txt" using 0:5 title ""
