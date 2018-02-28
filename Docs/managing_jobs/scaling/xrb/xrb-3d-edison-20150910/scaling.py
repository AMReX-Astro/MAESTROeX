import math
import numpy as np
import matplotlib.pyplot as plt

def scaling():

    plt.rcParams.update({'xtick.labelsize': 15,
                         'ytick.labelsize': 15,
                         'text.fontsize': 15})

    plt.rc("axes", linewidth=1.5)
    plt.rc("lines", markeredgewidth=1.5)

    data = np.loadtxt("scaling.txt")

    N_mpi = data[:,0]
    N_omp = data[:,1]

    t_adv = data[:,2]
    t_MAC = data[:,3]
    t_nodal = data[:,4]
    t_react = data[:,5]
    t_misc = data[:,6]
    t_total = data[:,7]
    std_total = data[:,8]
    comp = data[:,9]           # 1 == Intel; 2 == Cray

    cores = N_mpi * N_omp


    idx_intel = comp[:] == 1
    idx_MPI_intel = np.logical_and(comp[:] == 1, N_omp[:] == 1)
    idx_OMP_intel = np.logical_and(comp[:] == 1, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_intel], t_total[idx_MPI_intel], yerr=std_total[idx_MPI_intel],
                 marker="o", markeredgecolor="k", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="k", ls="none")

    plt.errorbar(cores[idx_OMP_intel], t_total[idx_OMP_intel], yerr=std_total[idx_OMP_intel],
                 marker="o", markeredgecolor="k", markersize=7, markeredgewidth=1.25,
                 color="k", ls="none")

    plt.plot(cores[idx_intel], t_total[idx_intel], color="k", ls="-", lw=1.5)


    idx_cray = comp[:] == 2
    idx_MPI_cray = np.logical_and(comp[:] == 2, N_omp[:] == 1)
    idx_OMP_cray = np.logical_and(comp[:] == 2, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_cray], t_total[idx_MPI_cray], yerr=std_total[idx_MPI_cray],
                 marker="^", markeredgecolor="r", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="r", ls="none")

    plt.errorbar(cores[idx_OMP_cray], t_total[idx_OMP_cray], yerr=std_total[idx_OMP_cray],
                 marker="^", markeredgecolor="r", markersize=7, markeredgewidth=1.25,
                 color="r", ls="none")

    plt.plot(cores[idx_cray], t_total[idx_cray], color="r", ls="-", lw=1.5)


    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")

    # ideal
    cm = np.min(cores[idx_intel])
    cM = np.max(cores[idx_intel])
    id = np.argmin(cores[idx_intel])
    plt.loglog([cm, cM], t_total[idx_intel][id]*cm/np.array([cm, cM]), ":", color="k")

    #cm = np.min(cores[idx_cray])
    #cM = np.max(cores[idx_cray])
    #id = np.argmin(cores[idx_cray])
    #plt.loglog([cm, cM], t_total[idx_cray][id]*cm/np.array([cm, cM]), ":", color="k")

    plt.text(300, 1.75, r"384$\times$384$\times$768 zones; 48$^3$ domain decomposition", fontsize=12)
    plt.text(300, 1.25, "2015-09-10", fontsize=12)
    plt.xlim(256, 32768)


    # custom legend
    legs = []
    legnames = []
    legs.append(plt.Line2D((0,1),(0,0), color="k", marker="o"))
    legnames.append(r"Intel compilers (15.0.1)")

    legs.append(plt.Line2D((0,1),(0,0), color="r", marker="^", markeredgecolor="r"))
    legnames.append(r"Cray compilers (8.4.0)")

    legs.append(plt.Line2D((0,1),(0,0), color="k",
                           marker="o", markeredgecolor="k", linestyle="none"))
    legnames.append("MPI + OpenMP")

    legs.append(plt.Line2D((0,1),(0,0), color="k",
                           marker="o", markeredgecolor="k", markerfacecolor="none", linestyle="none"))
    legnames.append("pure MPI")


    plt.legend(legs, legnames, ncol=2, frameon=False,
                 fontsize="small", numpoints=1)

    plt.xlabel("number of cores")
    plt.ylabel("average time to advance timestep")

    plt.title("NERSC Edison Scaling for Maestro 3-d XRB", fontsize="medium")

    plt.tight_layout()

    plt.savefig("edison_xrb_scaling.png")


if __name__== "__main__":
    scaling()
