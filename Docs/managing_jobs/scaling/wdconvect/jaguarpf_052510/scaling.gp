set term post eps color
set output 'scaling.eps'

set key top right

set title 'Weak Scaling Behavior of Full-Star MAESTRO Simulations With 1 Level of Refinement on jaguarpf'

set logscale x

set xrange [10:110000];

set yrange [0:200];

set pointsize 2;

set size 1, 1;

set xlabel "Number of Processors";
set ylabel "Average Time per Time Step (seconds)";

set xtics nomirror ("12" 12, "768" 768, "12k" 12000, "96k" 96000)

plot 'scaling.txt'   using 1:2 with linespoints lw 1 pt 5 title "2-Level Hybrid MPI/OpenMP, 12 Threads"
