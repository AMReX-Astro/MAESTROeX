import math
import numpy as np
import matplotlib.pyplot as plt

def scaling():

    plt.rcParams.update({'xtick.labelsize': 14,
                         'ytick.labelsize': 14,
                         'font.size': 14})

    plt.rc("axes", linewidth=1.5)
    plt.rc("lines", markeredgewidth=1.5)


    # titan data

    data = np.loadtxt("scaling-titan.txt")

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


    # we ran 2 tests, one that was 384 zones wide, and one that was
    # 768 zones wide we'll plot these in different colors, with hollow
    # symbols for pure MPI

    idx_384 = width[:] == 384
    idx_MPI_384 = np.logical_and(width[:] == 384, N_omp[:] == 1)
    idx_OMP_384 = np.logical_and(width[:] == 384, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_384], t_total[idx_MPI_384], yerr=std_total[idx_MPI_384],
                 marker="o", markeredgecolor="C0", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="C0", ls="none")

    plt.errorbar(cores[idx_OMP_384], t_total[idx_OMP_384], yerr=std_total[idx_OMP_384],
                 marker="o", markeredgecolor="C0", markersize=7, markeredgewidth=1.25,
                 color="C0", ls="none")

    plt.plot(cores[idx_384], t_total[idx_384], color="C0", ls="-", lw=1.5)


    idx_768 = width[:] == 768
    idx_MPI_768 = np.logical_and(width[:] == 768, N_omp[:] == 1)
    idx_OMP_768 = np.logical_and(width[:] == 768, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_768], t_total[idx_MPI_768], yerr=std_total[idx_MPI_768],
                 marker="^", markeredgecolor="C0", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="C0", ls="none")

    plt.errorbar(cores[idx_OMP_768], t_total[idx_OMP_768], yerr=std_total[idx_OMP_768],
                 marker="^", markeredgecolor="C0", markersize=7, markeredgewidth=1.25,
                 color="C0", ls="none")

    plt.plot(cores[idx_768], t_total[idx_768], color="C0", ls="-", lw=1.5)


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


    # edison data

    data = np.loadtxt("scaling-edison.txt")

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


    idx_cray = comp[:] == 2
    idx_MPI_cray = np.logical_and(comp[:] == 2, N_omp[:] == 1)
    idx_OMP_cray = np.logical_and(comp[:] == 2, N_omp[:] != 1)

    plt.errorbar(cores[idx_MPI_cray], t_total[idx_MPI_cray], yerr=std_total[idx_MPI_cray],
                 marker="o", markeredgecolor="C1", markerfacecolor="none",
                 markersize=7, markeredgewidth=1.25,
                 color="C1", ls="none")

    plt.errorbar(cores[idx_OMP_cray], t_total[idx_OMP_cray], yerr=std_total[idx_OMP_cray],
                 marker="o", markeredgecolor="C1", markersize=7, markeredgewidth=1.25,
                 color="C1", ls="none")

    plt.plot(cores[idx_cray], t_total[idx_cray], color="C1", ls="-", lw=1.5)


    # ideal
    cm = np.min(cores[idx_cray])
    cM = np.max(cores[idx_cray])
    id = np.argmin(cores[idx_cray])
    plt.loglog([cm, cM], t_total[idx_cray][id]*cm/np.array([cm, cM]), ":", color="k")






    plt.xlim(256, 131072)
    plt.ylim(1, 100)

    # custom legend
    legs = []
    legnames = []
    legs.append(plt.Line2D((0,1),(0,0), color="k", marker="o"))
    legnames.append(r"384$\times$384$\times$768")

    legs.append(plt.Line2D((0,1),(0,0), color="k", marker="^", markeredgecolor="k"))
    legnames.append(r"768$\times$768$\times$768")

    legs.append(plt.Line2D((0,1),(0,0), color="k",
                           marker="o", markeredgecolor="k", linestyle="none"))
    legnames.append("MPI + OpenMP")

    legs.append(plt.Line2D((0,1),(0,0), color="k",
                           marker="o", markeredgecolor="k", markerfacecolor="none", linestyle="none"))
    legnames.append("pure MPI")

    legs.append(plt.Line2D((0,1),(0,0), color="C0"))
    legnames.append(r"OLCF titan")

    legs.append(plt.Line2D((0,1),(0,0), color="C1"))
    legnames.append(r"NERSC edison")


    plt.legend(legs, legnames, ncol=3, frameon=False,
                 fontsize="small", numpoints=1)

    plt.xlabel("number of cores")
    plt.ylabel("average time to advance timestep")

    plt.title("Maestro 3-d XRB", fontsize="medium")

    plt.tight_layout()

    f = plt.gcf()
    f.set_size_inches(8, 6)
    plt.savefig("titan_edison_maestro_scaling.png")
    plt.savefig("titan_edison_maestro_scaling.pdf")


if __name__== "__main__":
    scaling()
