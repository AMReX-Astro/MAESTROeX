#!/usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from bisect import bisect_left
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infiles', type=str, nargs='+',
                    help='Name of diagnostic files to plot.')
parser.add_argument('-tmin', '--tmin', type=float, default=-1.0,
                    help='Minimum time to plot.')
parser.add_argument('-tmax', '--tmax', type=float, default=-1.0,
                    help='Maximum time to plot.')
args = parser.parse_args()



class Diagnostics(object):
    def __init__(self, file=None):
        self.file = file
        self.header = {}
        self.data = {}
        self.columns = []
        if self.file:
            self.read()

    def gendata(self):
        # Iterate over the data columns in the diagnostics file
        datacols = self.columns[:]
        datacols.remove('time')
        for dc in datacols:
            yield self.data['time'], self.data[dc], dc

    def plotall(self):
        # Plot all the data in the diagnostics file
        for time, data, col in self.gendata():
            self.plot(time, data, col)

    def plot(self, time, data, varlabel):
        # Plot variable 'varlabel' vs time
        fix, ax = plt.subplots(constrained_layout=True)
        print('Plotting {}'.format(varlabel))
        ax.plot(time, data)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel(varlabel)
        if args.tmin < 0.0:
            tmin = None
        if args.tmax < 0.0:
            tmax = None
        ax.set_xlim(left=tmin, right=tmax)
        file_label = ''.join(s for s in varlabel if s.isalnum())
        plt.savefig('diagnostics-{}.png'.format(file_label))
        plt.clf()

    def read(self):
        # Read lines of file
        f = open(self.file, 'r')
        assert(f)
        lines = []
        for l in f:
            ls = l.strip()
            if ls:
                lines.append(ls)
        f.close()

        # Read data from lines
        getheader = True        
        for l in lines:
            # Check for header (save the first one)
            if l[0]=='#':
                if getheader:
                    lh = l[2:]
                    if (not ':' in l) and ('time' in l):
                        # Get column headers
                        lhh = lh.split('  ')
                        for h in lhh:
                            hs = h.strip()
                            if hs:
                                self.data[hs] = []
                                self.columns.append(hs)
                    else:
                        # Get job metadata
                        h, v = lh.split(':', maxsplit=1)
                        h = h.strip()
                        v = v.strip()
                        self.header[h] = v
            # Read column data
            else:
                getheader = False
                timeindex = self.columns.index('time')
                ls = l.split()
                xvals = [float(lsi) for lsi in ls]
                xtime = xvals[timeindex]
                # print('----')
                # print('{}, time = {}'.format(len(self.data['time']), xtime))
                self.trimdata(xtime)
                # print('----')
                for c, x in zip(self.columns, xvals):
                    self.data[c].append(x)

        # Convert lists to arrays
        for k in self.data.keys():
            self.data[k] = np.array(self.data[k])

        # Print min/max times
        # print('min time: {}'.format(self.data['time'][0]))
        # print('max time: {}'.format(self.data['time'][-1]))

    def trimdata(self, trim_time):
        # Delete the data entries from trim_time onwards
        # print(self.data['time'])
        # print(trim_time)
        itrim = bisect_left(self.data['time'], trim_time)
        # print(itrim)
        for k in self.data.keys():
            self.data[k] = self.data[k][:itrim]
        # print(self.data['time'])

class SimulationDiagnostics(object):
    def __init__(self, files=[]):
        self.files = files
        self.diagnostics = [Diagnostics(f) for f in self.files]
        self.derived_diagnostics = Diagnostics()
        self.add_derived_fields()
        print('Added derived fields:')
        print(self.derived_diagnostics.columns)
        self.diagnostics.append(self.derived_diagnostics)

    def plotall(self):
        # Plot each of the columns against time
        for diag in self.diagnostics:
            diag.plotall()

    def add_derived_fields(self):
        # Add any derived fields, checking if their dependencies exist
        self.add_total_energy()

    def add_total_energy(self):
        xtime = np.array([])
        ntime = 'time'
        
        xgrav = np.array([])
        ngrav = 'grav pot energy'
        
        xeint = np.array([])
        neint = 'tot int energy'
        
        xekin = np.array([])
        nekin = 'tot kin energy'

        dependencies = {ntime:xtime,
                        ngrav:xgrav,
                        neint:xeint,
                        nekin:xekin}

        xtarget = np.array([])
        ntarget = 'tot energy'

        # Search for the dependency fields
        # in all the diagnostic files
        for kdep in dependencies.keys():
            if not dependencies[kdep].size:
                for diag in self.diagnostics:                
                    if kdep in diag.columns:
                        dependencies[kdep] = diag.data[kdep]
                        break

        # Check to make sure we found all the dependency fields
        for kdep in dependencies.keys():
            try:
                assert(dependencies[kdep].size)
            except:
                raise ValueError('Dependency field {} for derived field {} could not be found.'.format(kdep, ntarget))

        # Construct the derived field
        xtarget = dependencies[ngrav] + dependencies[neint] + dependencies[nekin]

        # Add time to derived field Diagnostics if not present
        if not ntime in self.derived_diagnostics.columns:
            self.derived_diagnostics.columns.append(ntime)
            self.derived_diagnostics.data[ntime] = dependencies[ntime]

        # Add derived field to derived field Diagnostics if not present
        if not ntarget in self.derived_diagnostics.columns:
            self.derived_diagnostics.columns.append(ntarget)
            self.derived_diagnostics.data[ntarget] = xtarget


if __name__ == "__main__":
    sdiag = SimulationDiagnostics(args.infiles)
    sdiag.plotall()
