import math
import numpy as np
import pylab 

def scaling():

    data = np.loadtxt("scaling.txt")

    problem_size = 384*384*768

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    jmode = data[:,8]

    cores = N_mpi * N_omp

    
    # we ran in 4 batches, with 1, 4, 8, and 16 threads -- separate them out
    # we also ran in 2 modes on titan -- with 8 cores per compute node and 16
    idx_omp1_j1 = np.logical_and(N_omp[:] == 1, jmode[:] == 1)
    idx_omp4_j1 = np.logical_and(N_omp[:] == 4, jmode[:] == 1)
    idx_omp8_j1 = np.logical_and(N_omp[:] == 8, jmode[:] == 1)
    idx_omp16_j1 = np.logical_and(N_omp[:] == 16, jmode[:] == 1)

    idx_omp1_j2 = np.logical_and(N_omp[:] == 1, jmode[:] == 2)
    idx_omp4_j2 = np.logical_and(N_omp[:] == 4, jmode[:] == 2)
    idx_omp8_j2 = np.logical_and(N_omp[:] == 8, jmode[:] == 2)
    idx_omp16_j2 = np.logical_and(N_omp[:] == 16, jmode[:] == 2)

    pylab.loglog(cores[idx_omp1_j1], t_total[idx_omp1_j1], "o-", color="k", label="MPI")
    pylab.loglog(cores[idx_omp4_j1], t_total[idx_omp4_j1], "o-", color="b", label="MPI + 4 OpenMP threads")
    pylab.loglog(cores[idx_omp8_j1], t_total[idx_omp8_j1], "o-", color="r", label="MPI + 8 OpenMP threads")
    pylab.loglog(cores[idx_omp16_j1], t_total[idx_omp16_j1], "o-", color="g", label="MPI + 16 OpenMP threads")

    pylab.loglog(cores[idx_omp1_j2], t_total[idx_omp1_j2], "^-", color="k")
    pylab.loglog(cores[idx_omp4_j2], t_total[idx_omp4_j2], "^-", color="b")
    pylab.loglog(cores[idx_omp8_j2], t_total[idx_omp8_j2], "^-", color="r")
    pylab.loglog(cores[idx_omp16_j2], t_total[idx_omp16_j2], "^-", color="g")

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)
    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")


    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)
    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB (384 x 384 x 768 zones)")

    pylab.ylim(1.,200.)

    pylab.savefig("xrb_titan_scaling_by_parallel.png")
    pylab.savefig("xrb_titan_scaling_by_parallel.eps")


    #=========================================================================
    # we also ran with 3 different grid sizes: 32^3, 48^3, and 64^3.
    # this is essentially N_mpi
    pylab.clf()

    grids = np.unique(N_mpi)
    jtype = [1, 2]
    for j in jtype:
        colors = ["g", "b", "r"]
        for g in grids:
            idx = np.logical_and(N_mpi[:] == g, jmode[:] == j)
            gsize = int(round((problem_size/g)**(1./3.)))
            if j == 1:
                pylab.loglog(cores[idx], t_total[idx], "o-", color=colors.pop(),
                             label=r"${}^3$ grid".format(gsize))
            else:
                pylab.loglog(cores[idx], t_total[idx], "^-", color=colors.pop())

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)

    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

    # custom legend
    legs = []
    legnames = []
    colors = ["g", "b", "r"]
    for g in grids:
        legs.append(pylab.Line2D((0,1),(0,0), color=colors.pop(), marker=None))
        gsize = int(round((problem_size/g)**(1./3.)))
        legnames.append(r"${}^3$ grid".format(gsize))


    # now the wider (768) run
    data = np.loadtxt("scaling-768.txt")

    problem_size = 768.**3

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    jmode = data[:,8]

    cores = N_mpi * N_omp


    grids = np.unique(N_mpi)
    jtype = [2]
    for j in jtype:
        colors = ["b"]
        for g in grids:
            idx = np.logical_and(N_mpi[:] == g, jmode[:] == j)
            gsize = int(round((problem_size/g)**(1./3.)))
            if j == 1:
                pylab.loglog(cores[idx], t_total[idx], "o-", color=colors.pop(),
                             label=None)
            else:
                pylab.loglog(cores[idx], t_total[idx], "^--", color=colors.pop())

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)

    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")


                    
    legs.append(pylab.Line2D((0,1),(0,0), color="k", linestyle="-"))    
    legnames.append(r"$384\times384\times768$")

    legs.append(pylab.Line2D((0,1),(0,0), color="k", linestyle="--"))    
    legnames.append(r"$768^3$")

    legs.append(pylab.Line2D((0,1),(0,0), color="k", marker="o", linestyle=""))    
    legnames.append("1 CPU / compute unit (-j 1)")

    legs.append(pylab.Line2D((0,1),(0,0), color="k", marker="^", linestyle=""))    
    legnames.append("2 CPUs / compute unit (-j 2)")

    pylab.legend(legs, legnames, ncol=3, frameon=False, 
                 fontsize=13, numpoints=1)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB")

    pylab.ylim(1.,200.)

    pylab.tight_layout()
    pylab.savefig("xrb_titan_scaling_by_grid.png")
    pylab.savefig("xrb_titan_scaling_by_grid.eps", bbox_inches="tight")



if __name__== "__main__":
    scaling()

