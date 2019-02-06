#!/usr/bin/env python

# Title: xin Generator (xinGen.py)
# Author: Adam Jacobs
# Creation Date: 04/20/2012

# Description: Very simple script for generating an input file of mass
# fractions for use by the test_react unit test.
#
# Revision History
# Programmer            Date                    Change
# ----------------------------------------------------------------------
# Adam Jacobs       04/20/2012              Code created.
# Adam Jacobs       05/16/2012              Updated interface.

# TODO:
# 1)Add more error checking.

###############
### Imports ###
###############
import sys

##################################
#### Global Data and Constants ###
##################################
# Solar values taken from Asplund, Grevesse, and Sauval (2005)
_X_SOL = 0.7393
_Y_SOL = 0.2485
_Z_SOL = 0.0122
_Z_C_SOL = 0.0031
_Z_O_SOL = 0.0081
_Z_N_SOL = 0.0010

################
#### Classes ###
################

##################
#### Functions ###
##################

#################
### Execution ###
#################

if __name__ == "__main__":
    # Get species count
    nspec = float(input('How many species are in the network? ' ))
    if(not (nspec % 1 == 0)):
        print('ERROR: Non-integer input for species count.')
        sys.exit(1)
    nspec = int(nspec)

    # Get grid size
    grid_size = float(input(
        'What is the grid size (cells on a side, e.g. 16, 32, or 64)? '))
    if(not (grid_size % 1 == 0)):
        print('ERROR: Non-integer input for grid size.')
        sys.exit(1)
    grid_size = int(grid_size)

    # Get mass fractions for each species
    xin = []
    for i in range(nspec):
        curx = []

        # Get user selection
        print('Select mass fraction type for species ', i, ': ')
        print('    A: Solar hydrogen value    [X=0.7393]')
        print('    B: Solar helium value      [Y=0.2485]')
        print('    C: Solar metallicity value [Z=0.0122]')
        print('    D: User-given uniform value')
        print('    E: User-given constant delta value')
        mft = input('Enter selection: ')

        # Solar X
        if mft.lower().startswith('a'):
            for j in range(grid_size):
                curx.append(_X_SOL)

        # Solar Y
        elif mft.lower().startswith('b'):
            for j in range(grid_size):
                curx.append(_Y_SOL)

        # Solar Z
        elif mft.lower().startswith('c'):
            for j in range(grid_size):
                curx.append(_Z_SOL)

        # Uniform
        elif mft.lower().startswith('d'):
            val = input('Enter uniform value: ')
            for j in range(grid_size):
                curx.append(val)

        # Delta
        elif mft.lower().startswith('e'):
            ini = input('Enter initial value (value at cell 0): ')
            delt = input(
                'Enter delta value (change in X from cell i to cell i+1): ')
            for j in range(grid_size):
                curx.append(ini + j * delt)

        # Invalid
        else:
            print('ERROR: Invalid selection.')
            sys.exit(1)

        xin.append(curx)

    # Write xin file
    label = 'xin'
    outfile = open(label, 'w')

    for i in range(nspec):
        for j in range(grid_size):
            outfile.write(str(xin[i][j]) + " ")
        outfile.write('\n')
