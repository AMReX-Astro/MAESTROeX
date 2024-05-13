#!/bin/bash
#######################################################################
# Description: AMREX's faverage takes only 1 variable, and writes to
# plotfile/.slice.  This script loops through faverage for a given list
# of variables, and creates a column separated plotfile/averages
#
# Usage: averages_vars.sh plotfile var1 var2 var3 ...
#######################################################################

faverage=$AMREX_HOME/Tools/Plotfile/faverage.gnu.x86-milan.ex
plotfile=$1
curr_dir=$(pwd)

# Append "/" to plotfile in case it's not there
if [ "${plotfile: -1}" != "/" ]; then
    plotfile+="/"
fi

# Check if plotfile exists
if [ ! -d "$plotfile" ]; then
    echo "Error: Input file '$1' does not exist."
    exit 1
fi

# restart log file
rm -f $plotfile/.faverage_log

# Header line
header="height "

i=0
#for var in "$@"; do
IFS=' ' read -ra vars <<< "${@:2}"
for var in "${vars[@]}"; do
    # Skip the first argument (the input file)
    if [ "$var" = "$1" ]; then
        continue
    fi

    #echo $var
    header+="$var "

    # Average current variable which creates (or overwrites) plotfile/.slice
    $faverage -v "$var" "$plotfile" >> $plotfile/.faverage_log

    # go into there for file manipulation
    cd $plotfile

    # backup for validating this script
    #cp .slice .slice_"$var"

    # Delete header line
    sed -i '1d' .slice # linux
    #sed -i '' '1d' .slice # BSD (mac)

    # If first variable, copy to plotfile/averages
    # otherwise, paste the second column into it
    if [ $i -eq 0 ]; then
        cp .slice averages
    else
        # extract second column and write to temp file
        cut -b 30- .slice > _temp1

        # paste and write to other temp file
        paste -d ' ' averages _temp1 > _temp2

        # overwrite averages file
        cp _temp2 averages
    fi

    rm -f _temp*
    rm -f .slice
    cd $curr_dir
    ((i++))
done

# Recreate the header line
sed -i '1s/^/'"$header"' \n/' $plotfile/averages # linux
#sed -i '' '1s/^/'"$header"' \n/' $plotfile/averages # BSD (mac)

