#!/bin/sh
slice-field.py $1 -res 2048 -w 1.5e8 -axis x
slice-field.py $1 -res 2048 -w 1.5e8 -axis y
slice-field.py $1 -res 2048 -w 1.5e8 -axis z

slice-field.py $1 -res 2048 -w 1.5e8 -axis x -f electron_fraction_asymmetry -min ' -3e-5' -max ' -0.5e-5' -cmap octarine
slice-field.py $1 -res 2048 -w 1.5e8 -axis y -f electron_fraction_asymmetry -min ' -3e-5' -max ' -0.5e-5' -cmap octarine
slice-field.py $1 -res 2048 -w 1.5e8 -axis z -f electron_fraction_asymmetry -min ' -3e-5' -max ' -0.5e-5' -cmap octarine

slice-enucdot.py $1 -w 1.5e8 -res 2048 -logmax 12 -symlog -rho 1.85e9 -axis x
slice-enucdot.py $1 -w 1.5e8 -res 2048 -logmax 12 -symlog -rho 1.85e9 -axis y
slice-enucdot.py $1 -w 1.5e8 -res 2048 -logmax 12 -symlog -rho 1.85e9 -axis z

slice-vel.py $1 -w 1.5e8 -rho 1.85e9 -res 2048 -dg -axis x
slice-vel.py $1 -w 1.5e8 -rho 1.85e9 -res 2048 -dg -axis y
slice-vel.py $1 -w 1.5e8 -rho 1.85e9 -res 2048 -dg -axis z

#plot-edot.py $1 -minr 2.5 -maxr 1500.0 -n 1000
