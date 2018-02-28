set term post eps color
set output 'scaling.eps'

set key top right

set title 'Strong Scaling Behavior of 768^3 MAESTRO Scientific Production Runs on intrepid.alcf.anl.gov'

set xrange [1500:60000];

set yrange [0:250];

set pointsize 2;

set size 1, 1;

set xlabel "Number of Processors";
set ylabel "Average Time per Time Step (seconds)";

set logscale x

set xtics nomirror ("1728" 1728, "5000" 5000, "10000" 10000, "16384" 16384, "55296" 55296)

set ytics nomirror 25

plot 'mpi.txt'   using 1:2 with linespoints lw 1 pt 5 title "MPI",\
     'omp4.txt'  using 1:2 with linespoints lw 1 pt 7 title "Hybrid MPI/OpenMP; 4 Threads"
