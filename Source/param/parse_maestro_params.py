#!/usr/bin/env python3

"""This script parses the list of C++ runtime parameters and writes
the necessary header and source files to make them available in
Maestro's C++ routines.

parameters have the format:

  name  type  default  ifdef

the first three (name, type, default) are mandatory:

  name: the name of the parameter.  This will be the same name as the
    variable in C++ unless a pair is specified as (name, cpp_name)

  type: the C++ data type (int, Real, bool, string)

  default: the default value.  If specified as a pair, (a, b), then
    the first value is the normal default and the second is for
    debug mode (#ifdef AMREX_DEBUG)

the next are optional:

   ifdef: only define this parameter if the name provided is #ifdef-ed

Any line beginning with a "#" is ignored

Commands begin with a "@":

   @namespace: sets the namespace that these will be under (see below)

     e.g. @namespace maestro

Note: categories listed in the input file aren't used for code generation
but are used for the documentation generation


For a namespace, name, we write out:

  -- name_params.H  (for maestro, included in Maestro.H):
     sets up the namespace and extern parameters

  -- name_queries.H  (for maestro, included in Maestro.cpp):
     does the parmparse query to override the default in C++

 -- name_job_info_tests.H
    this tests the current value against the default and outputs
    into a file

 -- runtime_params.cpp
    has the actual definition of the variables (without extern)

"""

import argparse
import re
import sys

import runtime_parameters as rp

CWARNING = """
// This file is automatically created by parse_maestro_params.py at build time.
// To update or add runtime parameters, please edit _cpp_parameters and rebuild.\n
"""

def read_param_file(infile):

    params = []

    namespace = None

    try:
        f = open(infile, encoding="UTF-8")
    except OSError:
        sys.exit("error opening the input file")

    for line in f:
        if line[0] == "#":
            continue

        if line.strip() == "":
            continue

        if line[0] == "@":
            # this is a command
            cmd, value = line.split(":")
            if cmd == "@namespace":
                fields = value.split()
                namespace = fields[0]

            else:
                sys.exit("invalid command")

            continue

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(r'".+"|[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

        name = fields[0]
        if name[0] == "(":
            name, cpp_var_name = re.findall(r"\w+", name)
        else:
            cpp_var_name = name

        dtype = fields[1].lower()

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try:
            ifdef = fields[3]
        except IndexError:
            ifdef = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(rp.Param(name, dtype, default,
                               cpp_var_name=cpp_var_name,
                               namespace=namespace,
                               debug_default=debug_default,
                               ifdef=ifdef))

    return params

def write_headers_and_source(params, out_directory):

    # output

    # find all the namespaces
    namespaces = sorted({q.namespace for q in params})

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]
        # sort by repr since None may be present
        ifdefs = sorted({q.ifdef for q in params_nm}, key=repr)

        # write name_params.H
        try:
            cp = open(f"{out_directory}/{nm}_params.H", "w", encoding="UTF-8")
        except OSError:
            sys.exit(f"unable to open {nm}_params.H for writing")

        cp.write(CWARNING)
        cp.write(f"#ifndef {nm.upper()}_PARAMS_H\n")
        cp.write(f"#define {nm.upper()}_PARAMS_H\n")

        cp.write("\n")
        cp.write(f"namespace {nm} {{\n")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cp.write(p.get_declare_string(with_extern=True))
            else:
                cp.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cp.write(p.get_declare_string(with_extern=True))
                cp.write("#endif\n")
        cp.write("}\n\n")
        cp.write("#endif\n")
        cp.close()

        # write name_queries.H
        try:
            cq = open(f"{out_directory}/{nm}_queries.H", "w", encoding="UTF-8")
        except OSError:
            sys.exit(f"unable to open {nm}_queries.H for writing")

        cq.write(CWARNING)

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cq.write(p.get_default_string())
                    cq.write(p.get_query_string())
                    cq.write("\n")
            else:
                cq.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cq.write(p.get_default_string())
                    cq.write(p.get_query_string())
                    cq.write("\n")
                cq.write("#endif\n")
            cq.write("\n")
        cq.close()

        # write the job info tests
        try:
            jo = open(f"{out_directory}/{nm}_job_info_tests.H", "w", encoding="UTF-8")
        except OSError:
            sys.exit(f"unable to open {nm}_job_info_tests.H")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    jo.write(p.get_job_info_test())
            else:
                jo.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    jo.write(p.get_job_info_test())
                jo.write("#endif\n")

        jo.close()

    # write a single C++ source file that actually defines the parameters
    # (one file for all namespaces)
    try:
        pf = open(f"{out_directory}/runtime_params.cpp", "w", encoding="UTF-8")
    except OSError:
        sys.exit(f"unable to open runtime_params.cpp")

    pf.write("#include <AMReX_REAL.H>\n")
    pf.write("#include <AMReX_Gpu.H>\n")

    for nm in namespaces:
        pf.write(f"#include <{nm}_params.H>\n")
    pf.write("\n")

    for nm in namespaces:
        params_nm = [q for q in params if q.namespace == nm]
        # sort by repr since None may be present
        ifdefs = sorted({q.ifdef for q in params_nm}, key=repr)

        pf.write(f"namespace {nm} {{\n")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    pf.write(p.get_declare_string())
            else:
                pf.write(f"#ifdef {ifdef}\n")
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    pf.write(p.get_declare_string())
                pf.write("#endif\n")
        pf.write("}\n\n")

    pf.close()


def main():
    """the main driver"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", type=str, default=None,
                        help="output directory for the generated files")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    p = read_param_file(args.input_file[0])
    write_headers_and_source(p, args.o)

if __name__ == "__main__":
    main()
