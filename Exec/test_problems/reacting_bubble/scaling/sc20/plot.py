import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

def plot():

    results_dir = './'

    results_files = [result for result in os.listdir(results_dir) if 'MAESTROeX' in result]

    n_gpus_per_node = 6

    throughput_list = []
    nnodes_list = []

    for results_file in results_files:

        nsteps = 0
        nzones = 0
        time = 0.0

        for line in open(results_dir + results_file):

            if len(line.split()) == 0:
                continue

            # Determine the number of MPI ranks and thus the number of nodes.

            if len(line.split()) == 6 and line.split()[0] == 'MPI' and line.split()[1] == 'initialized':
                n_ranks = int(line.split()[3])
                n_nodes = max(1, n_ranks / n_gpus_per_node)

            # For each step, add up the number of zones advanced and the walltime
            # for that step.

            if len(line.split()) == 4 and line.split()[0] == 'Level' and line.split()[1] == '0,' and line.split()[3] == 'cells':
                nsteps += 1
                nzones += int(line.split()[2])

            if len(line.split()) == 6 and line.split()[0] == 'Time' and line.split()[2] == 'advance':
                time += float(line.split()[5])

        nnodes_list.append(n_nodes)
        throughput_list.append(nzones / time / 1.e6)


    

    # Now we have all the results, so plot them.

    nnodes_arr = np.array(nnodes_list)
    throughput_arr = np.array(throughput_list)

    throughput_arr = np.array([x for _, x in sorted(zip(nnodes_arr, throughput_arr))])
    nnodes_arr = sorted(nnodes_arr)

    throughput_arr = throughput_arr / throughput_arr[0] / nnodes_arr

    plt.plot(nnodes_arr, throughput_arr, linestyle='-', lw=4, marker='o', markersize=14)

    plt.xlim([0.9 * min(nnodes_arr), 1.1 * max(nnodes_arr)])
    plt.ylim([0, 1.25 * max(throughput_arr)])

    plt.ylabel('Throughput (normalized)', fontsize=20)
    plt.xlabel('Number of nodes', fontsize=20)
    plt.title('Weak scaling of MAESTROeX reacting bubble', fontsize=20)
    plt.xscale('log', basex=2)
    ax = plt.gca()
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks([1, 2, 4, 8, 16, 32, 64, 128])
    plt.tick_params(labelsize=14)
    plt.tight_layout()

    plt.savefig('scaling.eps')
    plt.savefig('scaling.png')

def main():

    plot()

if __name__ == "__main__":

    main()
