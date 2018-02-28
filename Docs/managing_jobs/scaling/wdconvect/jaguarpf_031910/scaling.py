import math
import numpy
import pylab

def scaling():

    procs_mpi = numpy.zeros(4, numpy.int)
    avg_mpi = numpy.zeros(4, numpy.float64)

    procs_omp6 = numpy.zeros(4, numpy.int)
    avg_omp6 = numpy.zeros(4, numpy.float64)

    procs_omp12 = numpy.zeros(4, numpy.int)
    avg_omp12 = numpy.zeros(4, numpy.float64)

    
    # all the following is 768^3 timings, with one block per MPI task

    # pure MPI
    procs_mpi[:] = [1728, 4096, 13824, 32768]
    avg_mpi[:]   = [39.6, 21.1, 16.3,  30.3]

    # 6 OpenMP threads per MPI task
    procs_omp6[:] = [1296, 3072, 10368, 24576]
    avg_omp6[:]   = [60.1, 26.7, 9.9,   6.5]

    # 12 OpenMP threads per MPI task
    procs_omp12[:] = [2592, 6144, 20736, 49152]
    avg_omp12[:]   = [57.3, 22.5, 8.3,   4.9]


    ax = pylab.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    pylab.plot(procs_mpi, avg_mpi, 
               ms=6, mew=0, marker="s", c="r", ls="solid", label=r"pure MPI")

    ptemp = numpy.arange(2)*procs_mpi[len(procs_mpi)-1] + float(procs_mpi[0])
    pylab.plot(ptemp, avg_mpi[0]*(procs_mpi[0]/ptemp), "r:", label="_nolegend_")


    pylab.plot(procs_omp6, avg_omp6, 
               ms=6, mew=0, marker="o", c="b", ls="solid", label=r"6 OpenMP threads / MPI task")

    ptemp = numpy.arange(2)*procs_omp6[len(procs_omp6)-1] + float(procs_omp6[0])
    pylab.plot(ptemp, avg_omp6[0]*(procs_omp6[0]/ptemp), "b:", label="_nolegend_")


    pylab.plot(procs_omp12, avg_omp12, 
               ms=6, mew=0, marker="^", c="g", ls="solid", label=r"12 OpenMP threads / MPI task")

    ptemp = numpy.arange(2)*procs_omp12[len(procs_omp12)-1] + float(procs_omp12[0])
    pylab.plot(ptemp, avg_omp12[0]*(procs_omp12[0]/ptemp), "g:", label="_nolegend_")

    pylab.xlabel("# of processors")
    pylab.ylabel("time to advance timestep")


    leg = pylab.legend(loc=3)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)



    pylab.axis([1000,100000,1,100])



    pylab.title(r"MAESTRO WD convection (768^3) scaling on Jaguar XT5",fontsize=11)

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("jaguar_xt5_scaling.eps")



if __name__== "__main__":
    scaling()

