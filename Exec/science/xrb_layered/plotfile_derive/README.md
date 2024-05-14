# Derived quantities

This tool evaluates the following quantities from the maestro plotfiles.
There is some hard code in here which assumes that the plotfiles are 2D cartesian, g=1.29e14.
This is based on [amrex-astro-diag](https://github.com/amrex-astro/amrex-astro-diag)

### Thermodynamic gradients

- The true temperature gradient in the model (`del`):
$$\nabla=\left(\frac{d\ln T}{d\ln P}\right)$$

- The adiabatic gradient (`del_ad`):
$$\nabla_{\rm ad}=\left(\frac{d\ln T}{d\ln P}\right)_S=\frac{\chi_T}{\Gamma_1c_V}\frac{P}{\rho T}$$

- The Ledoux gradient (`del_ledoux`)::
$$\nabla_{\rm L}=\nabla_{\rm ad}+B$$
For the composition term $B$, we use MESA's formulation (Paxton et al. 2013, Equation 8), but with a centered difference instead. At some grid point $k$ with pressure $P_k=P(\rho_k,T_k,X_k)$:
$$B_k=-\frac{1}{\chi_T}\frac{\ln P(\rho_k,T_k,X_{k+1})-\ln P(\rho_k,T_k,X_{k-1})}{\ln P_{k+1}-\ln P_{k-1}}$$
(The numerator terms are the local pressure evaluated with the composition of the above/below grid points.)

### Vertical Velocity squared
$v_y^2$, variable name `vy2`. This is to compare to the MLT prediction

### Fluxes

- The true convective heat flux from the temperature excess `Fconv_dT`
$$F_{\rm conv}=\rho c_p v_y \delta T$$
where $\delta T=T-\bar{T}$, which exists as a code variable in MAESTROeX.
We also output it with $T$ instead of $\Delta T$ as `Fconv_T`. These should be equivalent formulations if $<v>=0$ because $<vdT> = <v(T-T0)> = <vT> - <v>T0 = <vT>$

The following expressions are from the mixing-length theory of Kippenhahn & Wiegert. We assume efficient convection ($\nabla_{\rm element}\approx\nabla_{\rm ad}$),  and assume the mixing length $\ell$ is equal to the pressure scale height $H_p=p/\rho g$. We therefore indicate the normalization in terms of the mixing-length parameter $\alpha_{\rm MLT}=\ell/H_p$).

- The *mixing-length theory* (MLT) convective heat flux `Fconv_mlt`
$$F_{\rm MLT}=\frac{1}{2}\rho c_pv_yT(\nabla-\nabla_{\rm ad})$$
Missing a factor $\alpha_{\rm MLT}$

- The MLT flux can also be expressed as a function of the velocity `Fconv_mlt_v`, under the same assumptions as above:
$$F_{\rm MLT,v}=\frac{4\rho c_pT}{\delta gH_p}v^3$$
where $\delta=\left(\frac{d\ln\rho}{d\ln T}\right)_P=\frac{\chi_T}{\chi_\rho}$ comes from the EOS. This ignores composition gradients.
Missing a factor $1/\alpha_{\rm MLT}$.
We also calculate it using $|v|^3$, as `Fconv_mlt_vabs`.

- This is using the following form for the square of the velocity, which can be checked independently from the flux. We also output it here even though it's not a flux, `v2_mlt`.
$$v^2_{\rm MLT}=g\delta(\nabla-\nabla_{\rm ad})H_p/8$$
Missing a factor $\alpha_{\rm MLT}^2$.

- The kinetic energy flux `Fkin`:
$$F_{\rm kin}=\frac{1}{2}\rho v^3$$
This should be similar as the MLT flux, up to some constant scaling (and using $v^3$ instead of $|v|^3$.

- The radiative flux `Frad`:
$$F_{\rm rad}=-\frac{4ac T^3}{3\kappa\rho}\frac{dT}{dr}$$
This requires that you include `CONDUCTIVITY_DIR` in the makefile

- The chemical fluxes of species $i$
$$F_{\rm X_i}=\rho v X_i$$
Currently, only the hydrogen and carbon fluxes (`Fh1`,`Fc12`) are included in the code, but it is straightforward to include other species.


### Other numbers
- The pressure scale height `Hp` as $p/\rho g$
- The heat capacity `cp`

- ~~The Damkholer number, which is the ratio between mixing and reaction timescales. We write it as the ratio between the convective turnover time and the nuclear timescale:
$$Da=\frac{H_p/v_y}{c_pT/H_{\rm nuc}}$$~~

- ~~The Richardson number, which is the ratio between buoyancy and shear terms:
$$Ri=\frac{g}{\rho}\frac{\drho/dy}{(dv_x/dy)^2}$$~~

Unfortunately, I did not include Hnuc and vx in the small plotfiles so these cannot be evaluated.


To build, do:

```
make DIM=2
```

changing the `DIM` line to match the dimension of your plotfile.

It is also important that the network you build with matches
the one used for generating the plotfile.  This is set via
the `NETWORK_DIR` parameter in the `GNUmakefile`.

Runtime parameters are managed by AMReX's ParmParse.  To run,
you specify the plotfile via `diag.plotfile`, either in an inputs
file or on the command line, e.g.:

```
./derived.gnu.ex diag.plotfile=plt00000
```



