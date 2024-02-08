'''
Read MESA model (.mod) file, and create an equivalent .raw file which can be
read by initial model routines. The following columns in the mod file are used:
radius, density, temperature, pressure, mass fractions

Uses the external package py_mesa_reader (github.com/wmwolf/py_mesa_reader)
'''

import sys
import mesa_reader as mr

def get_periodic_table_dict():
    '''
    Load the periodic_table.list file and put it in a dictionary
    Each entry has the key:item format 'he':helium' (all lowercase)
    '''
    atom_dict = {}
    with open("periodic_table.list", encoding="UTF-8") as file:
        for line in file:
            long,short = line.split()
            atom_dict[short.lower()] = long.lower()
    return atom_dict

def main(filename):
    ''' Load with mesa_reader, write new file '''

    is_modfile = True if '.mod' in filename else False

    data = mr.MesaData(filename)
    npts = len(data.T)

    # Figure out the species
    species = []
    for name in data.bulk_names:
        if any(char.isdigit() for char in name) and '_' not in name:
            if len(name.rstrip('0123456789'))<=2:  # this gets rid of reaction columns like "pnhe4"
                species.append(name)
    # print(species)

    # Create file
    new_filename = filename.replace(".mod","").replace(".data","") + ".raw"
    with open(new_filename, 'w', encoding="UTF-8") as new_file:

        nvar0 = 2 if is_modfile else 3 # mod files don't have pressure

        # Write the header
        new_file.write(f"# npts = {npts}\n")
        new_file.write(f"# num of variables = {nvar0 + len(species)}\n")
        new_file.write("# density\n# temperature\n")
        if not is_modfile:
            new_file.write("# pressure\n")

        # Species
        atom_dict = get_periodic_table_dict()
        for spec in species:
            key = spec.rstrip('0123456789')
            anum = spec[len(key):]
            new_file.write(f"# {atom_dict[key]}-{anum}\n")

        # Write the data
        # Mesa data goes outside->in, we want the opposite (smallest radius first)

        # At the very bottom of the layer (in the context of NS atmosphere), the compression
        # is so large that the radius does not change to 6 decimal places, so save more

        for i in range(len(data.R)):

            if is_modfile:
                line = f"{data.R[-1-i]:.8e} {data.d[-1-i]:.6e} {data.T[-1-i]:.6e} "
            else:
                line = f"{data.R_cm[-1-i]:.8e} {data.Rho[-1-i]:.6e} {data.T[-1-i]:.6e} {data.P[-1-i]:.6e} "

            for spec in species:
                line += f"{data.bulk_data[spec][-1-i]:.6e} "
            new_file.write(line+"\n")

    print("Saved to ", new_filename)


if __name__ == "__main__":

    if len(sys.argv)<2:
        sys.exit("Give filename")

    main(sys.argv[1])
