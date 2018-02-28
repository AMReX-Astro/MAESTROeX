set term post eps color
set output 'scaling.eps'

set key top right

set title 'Strong Scaling Behavior of 768^3 MAESTRO Scientific Production Runs on jaguarpf.ccs.ornl.gov'

set xrange [1000:50000];

set yrange [0:65];

set pointsize 2;

set size 1, 1;

set xlabel "Number of Processors";
set ylabel "Average Time per Time Step (seconds)";

set logscale x

set xtics nomirror ("1296" 1296, "5000" 5000, "10000" 10000, "24576" 24576, "32768" 32768, "49152" 49152)

set ytics nomirror 10

plot 'mpi.txt'   using 1:3 with linespoints lw 1 pt 5 title "MPI",\
     'omp6.txt'  using 1:3 with linespoints lw 1 pt 7 title "Hybrid MPI/OpenMP; 6 Threads", \
     'omp12.txt' using 1:3 with linespoints lw 1 pt 9 title "Hybrid MPI/OpenMP; 12 Threads"
