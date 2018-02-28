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
    width = data[:,9]

    cores = N_mpi * N_omp


    # we ran 2 tests, one that was 384 zones wide, and one that was 768 zones wide
    # we'll plot these in different colors, with hollow symbols for pure MPI

    idx_384 = width[:] == 384
    idx_MPI_384 = np.logical_and(width[:] == 384, N_omp[:] == 1)
    idx_OMP_384 = np.logical_and(width[:] == 384, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_384], t_total[idx_MPI_384], yerr=std_total[idx_MPI_384],
                 marker="o", markeredgecolor="k", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="k", ls="none")

    plt.errorbar(cores[idx_OMP_384], t_total[idx_OMP_384], yerr=std_total[idx_OMP_384],
                 marker="o", markeredgecolor="k", markersize=7, markeredgewidth=1.25,
                 color="k", ls="none")

    plt.plot(cores[idx_384], t_total[idx_384], color="k", ls="-", lw=1.5)


    idx_768 = width[:] == 768
    idx_MPI_768 = np.logical_and(width[:] == 768, N_omp[:] == 1)
    idx_OMP_768 = np.logical_and(width[:] == 768, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_768], t_total[idx_MPI_768], yerr=std_total[idx_MPI_768],
                 marker="^", markeredgecolor="r", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="r", ls="none")

    plt.errorbar(cores[idx_OMP_768], t_total[idx_OMP_768], yerr=std_total[idx_OMP_768],
                 marker="^", markeredgecolor="r", markersize=7, markeredgewidth=1.25,
                 color="r", ls="none")

    plt.plot(cores[idx_768], t_total[idx_768], color="r", ls="-", lw=1.5)


    ax = plt.gca()
    ax.set_xscale("log")
    ax.set_yscale("log")

    # ideal
    cm = np.min(cores[idx_384])
    cM = np.max(cores[idx_384])
    id = np.argmin(cores[idx_384])
    plt.loglog([cm, cM], t_total[idx_384][id]*cm/np.array([cm, cM]), ":", color="k")

    cm = np.min(cores[idx_768])
    cM = np.max(cores[idx_768])
    id = np.argmin(cores[idx_768])
    plt.loglog([cm, cM], t_total[idx_768][id]*cm/np.array([cm, cM]), ":", color="k")

    plt.text(400, 1.75, "48$^3$ domain decomposition", fontsize=14)
    plt.text(400, 1.25, "Cray 8.4.0 compilers; 2015-08-31", fontsize=14)
    plt.xlim(256, 131072)


    # custom legend
    legs = []
    legnames = []
    legs.append(plt.Line2D((0,1),(0,0), color="k", marker="o"))
    legnames.append(r"384$\times$384$\times$768")

    legs.append(plt.Line2D((0,1),(0,0), color="r", marker="^", markeredgecolor="r"))
    legnames.append(r"768$\times$768$\times$768")

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

    plt.title("OLCF Titan Scaling for Maestro 3-d XRB", fontsize="medium")

    plt.tight_layout()

    plt.savefig("titan_xrb_scaling.pdf")


if __name__== "__main__":
    scaling()
