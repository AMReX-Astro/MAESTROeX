import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# Load file
plotfile = sys.argv[1]
if plotfile[-1] != '/': plotfile+='/'
avgfile = plotfile + "averages"

header = np.loadtxt(avgfile, dtype="unicode", max_rows=1)
data = np.loadtxt(avgfile, skiprows=1).T
z = data[0]

# figure out subplots layout
Nvars = len(header)
Nrows = int(np.floor(np.sqrt(Nvars)))
Ncols = int(np.ceil(Nvars/Nrows))

# Plot
fig,axes = plt.subplots(Nrows,Ncols,figsize=(3*Ncols,2.5*Nrows))

for i,var in enumerate(header[1:]):
    ax = axes.ravel()[i] # ravel flattens the axes array
    ax.plot(z,data[i+1],'k-')
    ax.set_ylabel(var)

# remove the remaining subplots
i+=1
while i<Nrows*Ncols:
    ax = axes.ravel()[i].axis("off")
    i+=1

# x-label only on bottom right
axes[Nrows-1,0].set_xlabel("z (cm)")

if "/" in plotfile[:-1]: # it's a long path, strip the full path
    basename = os.path.basename(os.path.normpath(plotfile))
else:
    basename = plotfile[:-1]
#print(basename)

# title is model number
fig.suptitle(fr"Model {basename.split('_')[1]}")

plt.tight_layout()
#plt.show()

#if not os.path.exists("plt_avg/"): os.mkdir("plt_avg/")
fig.savefig("plt_avg/" + basename + ".png", bbox_inches="tight", dpi=300)
