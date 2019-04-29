import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("18m_500_s_rot_b_eq_1.hse.51200")

r = data[:,0]
rho = data[:,1]
T = data[:,2]

h1 = data[:,4]
he4 = data[:,6]
c12 = data[:,7]
o16 = data[:,9]
ne20 = data[:,10]
si28 = data[:,12]
s32 = data[:,13]
fe56 = data[:,21]
ni56 = data[:,22]

f = plt.figure()
ax_t = f.add_subplot(211)
ax_x = f.add_subplot(212)


ax_t.plot(r, rho, label=r"$\rho$")
ax_t.plot(r, T, label=r"$T$")

ax_t.set_xscale("log")
ax_t.set_yscale("log")

ax_t.legend(frameon=False, fontsize="small")


ax_x.plot(r, h1, label=r"${}^{1}\mathrm{H}$")
ax_x.plot(r, he4, label=r"${}^{4}\mathrm{He}$")
ax_x.plot(r, c12, label=r"${}^{12}\mathrm{C}$")
ax_x.plot(r, o16, label=r"${}^{16}\mathrm{O}$")
ax_x.plot(r, ne20, label=r"${}^{20}\mathrm{Ne}$")
ax_x.plot(r, si28, label=r"${}^{28}\mathrm{Si}$")
ax_x.plot(r, s32, label=r"${}^{32}\mathrm{S}$")
ax_x.plot(r, fe56, label=r"${}^{56}\mathrm{Fe}$")
ax_x.plot(r, ni56, label=r"${}^{56}\mathrm{Ni}$")

ax_x.set_xscale("log")
ax_x.set_yscale("log")

ax_x.legend(frameon=False, fontsize="small")

ax_x.set_ylim(1.e-5, 1.05)

ax_x.set_xlabel("r (cm)")

f.set_size_inches(7.0, 9.0)
f.savefig("model.png")
