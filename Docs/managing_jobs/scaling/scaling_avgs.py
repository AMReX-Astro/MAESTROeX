#!/usr/bin/env python

import argparse
import numpy as np

def get_time(line):
    idx = line.rfind("seconds")
    nline = line[:idx]
    return float(nline.split(":")[1])

def parser(N, ofile):

    f = open(ofile, "r")

    adv = []
    mac = []
    nodal = []
    react = []
    misc = []
    total = []

    for line in f:
        if "number of MPI processes" in line:
            n_mpi = int(line.split("=")[1])

        elif "number of threads" in line:
            n_omp = int(line.split("=")[1])

        elif "Advection       :" in line:
            adv.append(get_time(line))

        elif "MAC   Projection:" in line:
            mac.append(get_time(line))

        elif "Nodal Projection:" in line:
            nodal.append(get_time(line))

        elif "Reactions       :" in line:
            react.append(get_time(line))

        elif "Misc            :" in line:
            misc.append(get_time(line))

        elif "Time to advance timestep:" in line:
            total.append(get_time(line))

    if N == -1:
        N = len(adv)

    adv = np.array(adv[-N:])
    mac = np.array(mac[-N:])
    nodal = np.array(nodal[-N:])
    react = np.array(react[-N:])
    misc = np.array(misc[-N:])
    total = np.array(total[-N:])

    str = 9*"{:>10} "
    print "#  "+ str.format("MPI", "threads", "advection", "MAC",
                           "nodal", "reactions", "misc", "total", "total +/-")
    str = "{:10} {:10} {:>9.6f}  {:>9.6f}  {:>9.6f}  {:>9.6f}  {:>9.6f}  {:>9.6f}  {:>9.6f} "
    print "   "+str.format(n_mpi, n_omp, np.average(adv), np.average(mac), 
                           np.average(nodal), np.average(react), np.average(misc), np.average(total), np.std(total))

if __name__ == "__main__":
    p = argparse.ArgumentParser()

    p.add_argument("-N", help="number of timesteps (from the end) to include (-1 for all)",
                   type=int, default=10)

    p.add_argument("job_file", help="Maestro stdout file",
                   type=str, nargs=1)

    args = p.parse_args()

    parser(args.N, args.job_file[0])
