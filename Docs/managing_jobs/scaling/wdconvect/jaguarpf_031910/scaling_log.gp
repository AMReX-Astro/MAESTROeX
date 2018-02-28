set term post eps color
set output 'scaling_log.eps'

set key top right

set title 'Strong Scaling Behavior of 768^3 MAESTRO Scientific Production Runs on jaguarpf.ccs.ornl.gov'

set xrange [1000:55000];

set yrange [4:65];

set pointsize 2;

set size 1, 1;

set xlabel "Number of Processors";
set ylabel "Average Time per Time Step (seconds)";

set logscale x
set logscale y

set xtics nomirror ("2K" 2000, "4K" 4000, "14K" 14000, "33K" 33000, "49K" 49000)
set ytics nomirror ("5" 5, "" 6, "" 7, "" 8, "" 9, "10" 10, "20" 20, "" 30, "" 40, "50" 50)

plot 'mpi.txt'   using 1:3 w lp lt 1 lc 1 lw 1 pt 5 title "Pure MPI",\
     66182./x lt 2 lc 1 ti "Ideal Scaling", \
     'omp6.txt' using 1:3 w lp lt 1 lc 2 lw 1 pt 9 title "Hybrid MPI+OpenMP; 6 Threads", \
     76000./x lt 2 lc 2 ti "Ideal Scaling", \
     'omp12.txt' using 1:3 w lp lt 1 lc 3 lw 1 pt 9 title "Hybrid MPI+OpenMP; 12 Threads", \
     130118./x lt 2 lc 3 ti "Ideal Scaling"
