set terminal postscript eps enhanced color dashed lw 1 "Helvetica" 16

set output "acc.eps"

plot "acceptance_rate.txt" using 0:1 title ""
