#!/usr/bin/env python
"""
Given the output of fconv_slopes, plot the thermodynamic
gradients corresponding to an initial model.

Donald E. Willcox
"""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', type=str,
                    help='Name of file containing thermodynamic gradients to plot.')
parser.add_argument('-f', '--format', type=str, default='png',
                    help='Format of the desired output files. Can be, e.g. "png" or "eps". Defaults to "png".')
parser.add_argument('-rup', '--radius_upper', type=float,
                    help='Upper bound for the plotted radius.')
parser.add_argument('-o', '--outname', type=str, help='Base name of output file to use (w/o extension).')
args = parser.parse_args()

class ConvectiveGradients(object):
    def __init__(self, infile=None):
        if infile:
            self.r, self.actual, self.adiabatic, self.ledoux = np.loadtxt(infile, unpack=True)
            self.infile = infile
        else:
            self.r = []
            self.actual = []
            self.adiabatic = []
            self.ledoux = []
            self.infile = ''
        
    def plot(self, fmt=None, rup=None, outname=None, show=False):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        idxup = -1
        if rup:
            ax.set_xlim([0, rup])
            # Get the lowest index where radius > rup
            idxup = np.where(self.r > rup)[0][0]
        ax.set_xlabel('$\mathrm{r (cm)}$')
#        ax2 = ax.twinx()
        # ax.plot(self.r[:idxup], self.adiabatic[:idxup], color='blue', linestyle='-', label='adiabatic')
        # ax.plot(self.r[:idxup], self.actual[:idxup], color='green', linestyle='--', label='actual')
        # ax.plot(self.r[:idxup], self.ledoux[:idxup], color='red', linestyle=':', label='ledoux')
        dadiabatic = self.actual[:idxup]-self.adiabatic[:idxup]
        neg_idx, pos_idx = self.get_signed_indices(dadiabatic)
        # ax2.plot(self.r[:idxup][neg_idx], dadiabatic[neg_idx], color='black', marker='v', markersize=8,
        #          linestyle='-', label='actual-adiabatic (-)')
        # ax2.plot(self.r[:idxup][pos_idx], dadiabatic[pos_idx], color='black', marker='^', markersize=8,
        #          linestyle='-', label='actual-adiabatic (+)')
        dledoux = self.actual[:idxup]-self.ledoux[:idxup]        
        neg_idx, pos_idx = self.get_signed_indices(dadiabatic)
        # ax2.plot(self.r[:idxup][neg_idx], dledoux[neg_idx], color='magenta', marker='v', markersize=8,
        #          linestyle=':', label='actual-ledoux (-)')
        # ax2.plot(self.r[:idxup][pos_idx], dledoux[pos_idx], color='magenta', marker='^', markersize=8,
        #          linestyle=':', label='actual-ledoux (+)')
        ax.plot(self.r[:idxup], dadiabatic, color='blue', linestyle='-', label='adiabatic $\mathrm{\\nabla_{conv}}$')
        ax.plot(self.r[:idxup], dledoux, color='red', linestyle='-.', label='ledoux $\mathrm{\\nabla_{conv}}$')
        mx = max(np.amax(dadiabatic), np.amax(dledoux))
        mn = min(np.amin(dadiabatic), np.amin(dledoux))
        mlin = min(abs(mx), abs(mn))
        plt.yscale('symlog', linthreshy=0.5*mlin)
        ax.set_ylabel('$\mathrm{\\nabla_{actual} - \\nabla_{conv}}$')
        plt.legend()
        if fmt=='png':
            if not outname:
                outname = self.infile + '.png'
            plt.savefig(outname, dpi=300)
        else:
            if not outname:
                outname = self.infile + '.eps'
            plt.savefig(outname)
        if show:
            plt.show()
            plt.close(fig)

    def get_signed_indices(self, dvec):
        neg_idx = np.where(dvec < 0.0)
        pos_idx = np.where(dvec > 0.0)
        return neg_idx, pos_idx

if __name__=='__main__':
    cg = ConvectiveGradients(args.infile)
    cg.plot(args.format, args.radius_upper, args.outname)

    
