"""
This script strips out preprocessor directives from C++ headers and Fortran
files and saves the results in source/preprocessed_files
"""

import os
import re

# directory of the source files
rootdir = "../Source"

outdir = "source/preprocessed_files"


def strip_directives(filename, filepath, outpath):
    """
    Read in file, remove all preprocessor directives and output
    """

    # r = re.compile(r"(^#.*$\n)")

    with open(os.path.join(filepath, filename)) as infile:
        txt = infile.read()

        outtxt = re.sub(r"(^#.*$\n)", '', txt, flags=re.M)

        with open(os.path.join(outpath, filename), 'w') as outfile:
            outfile.write(outtxt)


def delete_lines(filename, filepath):
    """
    For some reason sphinx-fortran does not like it when there are ' in the
    middle of strings, so we shall remove all of these.

    It also gets confused by 'double precision', so I'm going to replace all of
    these with 'real'
    """

    txt = ""
    with open(os.path.join(filepath, filename)) as infile:
        txt = infile.read()

    txt = re.sub(r"(call log\(\".*\")", "", txt, flags=re.M)
    txt = re.sub(r"double precision", "real", txt, flags=re.M)

    with open(os.path.join(filepath, filename), 'w') as outfile:
        outfile.write(txt)


if __name__ == "__main__":

    excl_files = []

    # make the output directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # loop over files in rootdir
    for f in sorted(os.listdir(rootdir)):
        if ((f[-2:] == ".H" and f[-4:] != "_F.H") or f[-4:] == ".F90"
                or f[-4:] == ".f90") and f not in excl_files:
            strip_directives(f, rootdir, outdir)
            delete_lines(f, outdir)

    # loop over source dir
    for subdir in sorted(os.listdir(rootdir)):
        if not os.path.isdir(os.path.join(rootdir, subdir)):
            continue

        # loop over files in subdirectories and run strip_directives on all
        # C++ header files and Fortran files
        for f in sorted(os.listdir(os.path.join(rootdir, subdir))):
            if ((f[-2:] == ".H" and f[-4:] != "_F.H") or f[-4:] == ".F90"
                    or f[-4:] == ".f90") and f not in excl_files:
                strip_directives(f, os.path.join(rootdir, subdir), outdir)
                delete_lines(f, outdir)
