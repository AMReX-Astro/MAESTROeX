#!/usr/bin/env python
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str,
                    help='Name of neutrino spectrum file to plot.')
parser.add_argument('-s', '--sumonly', action='store_true',
                    help='If supplied, only print the sum of the weights and quit.')
args = parser.parse_args()


# Read file
f = open(args.infile, 'r')
f.readline()
f.readline()
header = f.readline().strip().split()

data = {}
for h in header:
    data[h] = []

for line in f:
    ls = line.strip().split()
    for h, lsi in zip(header, ls):
        data[h].append(float(lsi))
f.close()

# Convert to numpy arrays
for h in header:
    data[h] = np.array(data[h])

if not args.sumonly:
    # Electron capture spectrum
    fig, ax = plt.subplots(constrained_layout=True)
    ax.semilogy(data['Energy'], data['Lambda_ecap23'], color='blue', label='ecap23')
    ax.set_xlim([0,1.5])
    ax.set_xlabel('$E_\\nu$ (MeV)')
    ax.set_ylabel('$f(E)$ (1/s)')
    plt.savefig('spectrum_ecap23.png', dpi=600)
    plt.clf()

    # Beta decay spectrum
    fig, ax = plt.subplots(constrained_layout=True)
    ax.semilogy(data['Energy'], data['Lambda_beta23'], color='green', label='beta23')
    ax.set_xlim([0,1.5])
    ax.set_xlabel('$E_\\nu$ (MeV)')
    ax.set_ylabel('$f(E)$ (1/s)')
    plt.savefig('spectrum_beta23.png', dpi=600)
    plt.clf()

    # Total spectrum
    fig, ax = plt.subplots(constrained_layout=True)
    ax.semilogy(data['Energy'], data['Lambda_ecap23'], color='blue', label='ecap23')
    ax.semilogy(data['Energy'], data['Lambda_beta23'], color='green', label='beta23')
    ax.semilogy(data['Energy'], data['Lambda_total'], color='black', label='total')
    ax.set_xlim([0,1.5])
    ax.set_xlabel('$E_\\nu$ (MeV)')
    ax.set_ylabel('$f(E)$ (1/s)')
    plt.savefig('spectrum.png', dpi=600)
    plt.clf()

number_rate_ecap23 = np.sum(data['Lambda_ecap23'])
number_rate_beta23 = np.sum(data['Lambda_beta23'])

print('Number rate for ecap23: {}'.format(number_rate_ecap23))
print('Number rate for beta23: {}'.format(number_rate_beta23))
print('Total number rate: {}'.format(number_rate_ecap23 + number_rate_beta23))

