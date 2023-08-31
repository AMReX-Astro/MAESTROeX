import numpy as np
import matplotlib.pyplot as plt

import mesa.py_mesa_reader as mr

fig,axes = plt.subplots(3,1,figsize=(6,8), sharex=True)
fig.subplots_adjust(hspace=0)

axes[0].set_ylabel(r"$T$ (K)")
axes[1].set_ylabel(r"$\rho$ (g cm$^{-3}$)")
axes[2].set_ylabel(r"$X$")
axes[2].set_xlabel(r"r (cm)")
axes[2].set_ylim([1e-3,1.1])

# MESA model
prof = mr.MesaData("mesa/profile37.data")

axes[0].plot(prof.R_cm-12e5, prof.T, 'k-', label="MESA")
axes[1].semilogy(prof.R_cm-12e5, prof.Rho, 'k-')

axes[2].semilogy(prof.R_cm-12e5, prof.h1, 'b-',label="h1")
axes[2].semilogy(prof.R_cm-12e5, prof.he4, 'r-',label="he4")
cno = prof.c12
for iso in ("c13","n13","n14","n15","o14","o15","o16","o17","o18","f17","f18","f19"):
    cno += prof.data(iso)
axes[2].semilogy(prof.R_cm-12e5, cno, 'g-',label="CNO")
axes[2].semilogy(prof.R_cm-12e5, prof.fe56, 'k-',label="fe56")

# toy atm model
all_data = np.loadtxt("toy_atm/toy_atm_3cm.hse", skiprows=26).T
names = ["density","temperature","pressure","hydrogen-1","helium-4","carbon-12","carbon-13","nitrogen-13","nitrogen-14","nitrogen-15","oxygen-14","oxygen-15","oxygen-16","oxygen-17","oxygen-18","fluorine-17","fluorine-18","fluorine-19","neon-18","neon-19","neon-20","magnesium-22","magnesium-24","iron-56"]
data = {name:all_data[i+1] for i,name in enumerate(names)}
z = all_data[0]

axes[0].plot(z, data["temperature"], 'k--', label="toy atm")
axes[1].semilogy(z, data["density"], 'k--')

axes[2].semilogy(z, data["hydrogen-1"], 'b--')
axes[2].semilogy(z, data["helium-4"], 'r--')
axes[2].semilogy(z, data["carbon-12"]+data["oxygen-14"]+data["oxygen-15"], 'g--')
axes[2].semilogy(z, data["iron-56"], 'k--')

axes[0].legend(frameon=False, ncol=2, bbox_to_anchor = (0.8,1.15), bbox_transform=axes[0].transAxes)
axes[2].legend(loc=3)

fig.savefig("toy_atm.png", bbox_inches="tight", dpi=500)
