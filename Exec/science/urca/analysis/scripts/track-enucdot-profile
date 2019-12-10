#!/usr/bin/env python
"""
Track the maximum enucdot location and FWHM in a 1D profile for the flame problem.

Donald E. Willcox
"""
import argparse
import re
import os
import yt

parser = argparse.ArgumentParser()
parser.add_argument('-re', '--regexp_file_names', type=str,
                    help='Python regular expression which files must match.')
parser.add_argument('-ifile', '--input_file_names', type=str,
                    help=('Name of file containing file names, one on each line, ' +
                          'to use as {file} entries in the template.'))
parser.add_argument('-o', '--out_task_file', type=str, default='enucdot.dat',
                    help='Name of the output file containing the profile data.')
parser.add_argument('-f', '--fields', type=str, nargs='+', default='enucdot',
                    help='Name of the fields to profile. Eg. "enucdot".')
parser.add_argument('-bf', '--binfields', type=str, nargs='+', default='radius',
                    help='Bin axis (field) along which to profile the desired field. Default is "x".')
parser.add_argument('-nbins', '--nbins', type=int, nargs='+', default=64,
                    help='Number of bins in each of the --binfields.')
parser.add_argument('-xmin', '--xmin', type=float, default=0.0,
                    help='Ignore cells with x < xmin. Default=0.0.')
parser.add_argument('-ymin', '--ymin', type=float, default=0.0,
                    help='Ignore cells with y < ymin. Default=0.0.')
parser.add_argument('-zmin', '--zmin', type=float, default=0.0,
                    help='Ignore cells with z < zmin. Default=0.0.')
args = parser.parse_args()

def doit(filename):
    ds = yt.load(filename)
    ad = ds.all_data()
    profile = yt.create_profile(ad, args.binfields, args.fields,
                                n_bins=tuple(args.nbins),
                                weight_field=None)

    print(profile.keys())
    print(profile["enucdot"].shape)
    print(profile["y"].shape)
    for y, e in zip(profile["y"], profile["enucdot"]):
        print('({}, {})'.format(y, e))
    exit()
    return ds.time, emax, elocmax, efwhm

if __name__ == '__main__':
    if not (args.regexp_file_names or args.input_file_names):
        print('ERROR: either -re or -ifile options must be supplied.')
        
    if args.regexp_file_names:
        # Sanity Check: compile the regular expression
        regexp = re.compile(args.regexp_file_names)

    subject_files = []

    # Get input file names if supplied from a file
    if args.input_file_names:
        fif = open(args.input_file_names, 'r')
        for l in fif:
            subject_files.append(l.strip())
        fif.close()

    # Get input file names if targeted by regular expression
    if args.regexp_file_names:
        # Get the list of subject files in the current directory
        ## Files in the current directory
        cwd_files = os.listdir()
        ## Filter files by regexp and skip if {file}.skip exists
        subject_files = []
        for f in cwd_files:
            if regexp.match(f):
                if not f + '.skip' in cwd_files:
                    subject_files.append(f)

    # Sort subject files
    try:
        subject_files = sorted(subject_files, key=lambda x: int(x[-6:]))
    except:
        subject_files = sorted(subject_files, key=lambda x: int(x[-5:]))

    # Loop over the subject files and compute needed values
    ds_time = []
    enucdot_max = []
    enucdot_loc_max = []
    enucdot_fwhm = []

    for f in subject_files:
        t, emax, elocmax, efwhm = doit(f)
        ds_time.append(t)
        enucdot_max.append(emax)
        enucdot_loc_max.append(elocmax)
        enucdot_fwhm.append(efwhm)

    # Save data in an output profile
    fo = open(args.out_task_file, 'w')
    fo.write('t     enucdot_max     enucdot_loc_max     enucdot_fwhm\n')
    for t, emax, elocmax, efwhm in zip(ds_time, enucdot_max,
                                       enucdot_loc_max, enucdot_fwhm):
        fo.write('{}  {}  {}  {}\n'.format(t, emax, elocmax, efwhm))
    fo.close()
