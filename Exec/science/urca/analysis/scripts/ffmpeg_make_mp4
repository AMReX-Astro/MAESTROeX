#!/usr/bin/env python
""" 
Use ffpmeg to transform a set of png files to an mp4 video.

Donald E. Willcox
"""
import os
import argparse
from FileToolkit import sort_input_filenames

parser = argparse.ArgumentParser()
parser.add_argument('infiles', type=str, nargs='+',
                    help='List of input files for encoding into a movie.')
parser.add_argument('-s', '--sort', action='store_true',
                    help='Sort input files sequentially by interpreting as an integer the first contiguous subset of the filename strings where the filenames differ. Does not strip leading paths.')
parser.add_argument('-name', '--moviename', type=str,
                    default='movie',
                    help='Name of output (not extension). Default is "movie.mp4".')
parser.add_argument('-ifps', '--input_fps', type=int, default=15,
                    help='Input rate at which to read images, i.e. number of pngs per second. Default is 15.')
parser.add_argument('-ofps', '--output_fps', type=int, default=25,
                    help='Output rate of movie frames (not the rate of plot images). Default is 25.')
args = parser.parse_args()


def symlink_inputs(base_link, filenames):
    """Create temporary symbolic links for the filenames."""

    symnames = []
    for i, f in enumerate(filenames):
        symn = '{}_{:010d}.png'.format(base_link, i)
        os.symlink(f, symn)
        symnames.append(symn)
    return symnames

def remove_symlinks(symnames):
    """Delete the temporary symbolic links."""

    for sf in symnames:
        os.remove(sf)

if __name__=='__main__':
    # Sort input files if needed
    if args.sort:
        input_files = sort_input_filenames(args.infiles)
    else:
        input_files = args.infiles

    # Symbolically link input files to a series of temporary files
    base_link = '__fftemp'
    symnames  = symlink_inputs(base_link, input_files)

    # Run ffmpeg on the files specified
    cmd = 'ffmpeg -r {} -i __fftemp_%010d.png -s 1440x1080 -vcodec libx264 -crf {} -pix_fmt yuv420p {}.mp4'.format(args.input_fps,
                                                                                                                   args.output_fps,
                                                                                                                   args.moviename)
    os.system(cmd)

    # Remove temporary files
    remove_symlinks(symnames)
