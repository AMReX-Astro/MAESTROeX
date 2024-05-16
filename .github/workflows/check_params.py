import argparse
import importlib
import os
import re
from pathlib import Path
import sys

# some parameters that are not defined in _cpp_parameters

whitelist = ["maestro.lo_bc",
             "maestro.hi_bc"]

# we don't have all of the radiation parametrs in the _cpp_parameters
# yet, so we won't check these namespaces

namespace_ignore = []

def doit(maestroex_dir):

    maestroex = Path(os.path.abspath(maestroex_dir))

    # import the module that defines the MAESTROeX runtime params
    sys.path.append(str(maestroex / "Source" / "param/"))
    import parse_maestro_params

    # read in the parameters defined in _cpp_parameters
    param_file = maestroex / "Source" / "param" / "_cpp_parameters"
    params = parse_maestro_params.read_param_file(str(param_file))

    namespaces = set(p.namespace for p in params)
    runtime_parameters = [f"{p.namespace}.{p.name}" for p in params]

    pattern = re.compile(r"[A-Za-z0-9_]+\.[A-Za-z0-9_]+", re.IGNORECASE)

    # loop over all the inputs files
    exec_path = maestroex / "Exec"
    for f in exec_path.glob("**/inputs*"):

        if os.path.isdir(f):
            continue

        # find all the params in each namespace
        with open(f) as infile:
            print(f"working on {f}")
            for line in infile:
                # remove comments
                idx = line.find("#")
                if idx > 0:
                    line = line[:idx]

                found_param = pattern.match(line)
                if not found_param:
                    continue

                p = found_param.group(0)
                nm = p.split(".")[0]
                if nm in namespaces and nm not in namespace_ignore:
                    if not (p in runtime_parameters or p in whitelist):
                        sys.exit(f"Error: {p} not valid")


if __name__ == "__main__":

    # we need the top-level Castro directory

    p = argparse.ArgumentParser()
    p.add_argument("maestroex_dir", type=str, nargs=1,
                   help="top level MAESTROeX directory")

    args = p.parse_args()

    doit(args.maestroex_dir[0])


