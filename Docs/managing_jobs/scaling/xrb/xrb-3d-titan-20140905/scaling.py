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

    cores = N_mpi * N_omp

    
    gsize = 48
    pylab.loglog(cores, t_total, "o-", color="r",
                 label=r"${}^3$ grid".format(gsize))


    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)

    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")

 
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

    cores = N_mpi * N_omp


    gsize = 48
    pylab.loglog(cores, t_total, "o-", color="r",
                 label=None)

    # ideal
    cm = np.min(cores)
    cM = np.max(cores)
    id = np.argmin(cores)

    pylab.loglog([cm, cM], t_total[id]*cm/np.array([cm, cM]), ":", color="k", label="ideal scaling")


    legs = []
    legnames = []

    legs.append(pylab.Line2D((0,1),(0,0), color="k", linestyle="-"))    
    legnames.append(r"$384\times384\times768$")

    legs.append(pylab.Line2D((0,1),(0,0), color="k", linestyle="--"))    
    legnames.append(r"$768^3$")

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB")

    pylab.xlim(100,1.e5)
    pylab.ylim(1.,200.)

    pylab.tight_layout()
    pylab.savefig("xrb_titan_scaling_by_grid.png")
    pylab.savefig("xrb_titan_scaling_by_grid.eps", bbox_inches="tight")



if __name__== "__main__":
    scaling()

