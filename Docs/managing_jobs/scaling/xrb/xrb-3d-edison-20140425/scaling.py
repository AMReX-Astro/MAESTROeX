import math
import numpy as np
import pylab 

def scaling():

    data = np.loadtxt("scaling.txt")

    problem_size = 384*384*768

    print data.shape

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]

    cores = N_mpi * N_omp

    # we ran in 3 batches, with 1, 6, and 12 threads -- separate them out
    idx_omp1 = N_omp[:] == 1
    idx_omp6 = N_omp[:] == 6
    idx_omp12 = N_omp[:] == 12

    pylab.loglog(cores[idx_omp1], t_total[idx_omp1], "o-", color="b", label="pure MPI")
    pylab.loglog(cores[idx_omp6], t_total[idx_omp6], "o-", color="r", label="MPI + 6 OpenMP threads")
    pylab.loglog(cores[idx_omp12], t_total[idx_omp12], "o-", color="g", label="MPI + 12 OpenMP threads")

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    pylab.loglog([cm, cM], t_total[0]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("NERSC Edison Scaling for 3-d XRB (384 x 384 x 768 zones)")

    pylab.ylim(1.,100.)

    pylab.savefig("xrb_edison_scaling_by_parallel.png")

    # we also ran with 3 different grid sizes: 32^3, 48^3, and 64^3.
    # this is essentially N_mpi
    pylab.clf()

    grids = np.unique(N_mpi)
    colors = ["g", "b", "r"]
    for g in grids:
        idx = N_mpi[:] == g
        gsize = int(round((problem_size/g)**(1./3.)))
        pylab.loglog(cores[idx], t_total[idx], "o-", color=colors.pop(),
                     label=r"${}^3$ grid".format(gsize))


    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    pylab.loglog([cm, cM], t_total[0]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("NERSC Edison Scaling for 3-d XRB (384 x 384 x 768 zones)")

    pylab.ylim(1.,100.)
    pylab.savefig("xrb_edison_scaling_by_grid.eps")



if __name__== "__main__":
    scaling()

