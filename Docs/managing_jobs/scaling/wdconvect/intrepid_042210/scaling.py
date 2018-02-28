import math
import numpy
import pylab

def scaling():

    procs_mpi = numpy.zeros(3, numpy.int)
    avg_mpi = numpy.zeros(3, numpy.float64)

    procs_omp4 = numpy.zeros(4, numpy.int)
    avg_omp4 = numpy.zeros(4, numpy.float64)

    
    # all the following is 768^3 timings, with one block per MPI task

    # pure MPI
    procs_mpi[:] = [1728,  4096,  13824]
    avg_mpi[:]   = [227.9, 100.0, 68.9]

    # 4 OpenMP threads per MPI task
    procs_omp4[:] = [2048,  6912, 16384, 55296]
    avg_omp4[:]   = [211.5, 70.6, 35.1,  46.4]


    ax = pylab.subplot(111)
    ax.set_xscale('log')
    ax.set_yscale('log')

    pylab.plot(procs_mpi, avg_mpi, 
               ms=6, mew=0, marker="s", c="r", ls="solid", label=r"pure MPI")

    ptemp = numpy.arange(2)*procs_mpi[len(procs_mpi)-1] + float(procs_mpi[0])
    pylab.plot(ptemp, avg_mpi[0]*(procs_mpi[0]/ptemp), "r:", label="_nolegend_")


    pylab.plot(procs_omp4, avg_omp4, 
               ms=6, mew=0, marker="o", c="b", ls="solid", label=r"4 OpenMP threads / MPI task")

    ptemp = numpy.arange(2)*procs_omp4[len(procs_omp4)-1] + float(procs_omp4[0])
    pylab.plot(ptemp, avg_omp4[0]*(procs_omp4[0]/ptemp), "b:", label="_nolegend_")


    pylab.xlabel("# of processors")
    pylab.ylabel("time to advance timestep")


    leg = pylab.legend(loc=3)
    ltext = leg.get_texts()
    pylab.setp(ltext, fontsize='small')
    leg.draw_frame(0)



    pylab.axis([1000,100000,10,300])



    pylab.title(r"MAESTRO WD convection (768^3) scaling on Intrepid (BG/P)",fontsize=11)

    f = pylab.gcf()
    f.set_size_inches(6.0,6.0)

    pylab.savefig("intrepid_BGP_scaling.eps")



if __name__== "__main__":
    scaling()

