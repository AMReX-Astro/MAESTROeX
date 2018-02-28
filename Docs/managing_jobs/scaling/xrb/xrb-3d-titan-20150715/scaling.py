import math
import numpy as np
import pylab 

def scaling():

    data = np.loadtxt("scaling.txt")

    N_mpi = data[:,0]
    N_omp = data[:,1]
    
    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    width = data[:,8]

    cores = N_mpi * N_omp

    
    # we ran 2 tests, one that was 384 zones wide, and one that was 768 zones wide
    idx_384 = width[:] == 384
    idx_768 = width[:] == 768

    pylab.loglog(cores[idx_384], t_total[idx_384], "o-", color="k", label=r"$384\times 384\times 768$")
    pylab.loglog(cores[idx_768], t_total[idx_768], "^-", color="r", label=r"$768\times 768\times 768$")

    # ideal
    cm = np.min(cores[idx_384])
    cM = np.max(cores[idx_384])
    id = np.argmin(cores[idx_384])
    pylab.loglog([cm, cM], t_total[idx_384][id]*cm/np.array([cm, cM]), ":", color="k")

    cm = np.min(cores[idx_768])
    cM = np.max(cores[idx_768])
    id = np.argmin(cores[idx_768])
    pylab.loglog([cm, cM], t_total[idx_768][id]*cm/np.array([cm, cM]), ":", color="k")

    pylab.text(600, 1.25, "Cray 8.3.9 compilers; 2015-07-15")
    pylab.xlim(512, 131072)

    pylab.legend(frameon=False)

    pylab.xlabel("number of cores")
    pylab.ylabel("average time to advance timestep")

    pylab.title("OLCF Titan Scaling for 3-d XRB")

    pylab.tight_layout()

    pylab.savefig("titan_xrb_scaling.png")


    pylab.loglog(cores[idx_384], t_react[idx_384], "o--", color="0.5", label=r"$384\times 384\times 768$ reactions")
    pylab.loglog(cores[idx_768], t_react[idx_768], "^--", color="m", label=r"$384\times 384\times 768$ reactions")

    pylab.savefig("titan_xrb_scaling_react.png")

if __name__== "__main__":
    scaling()

