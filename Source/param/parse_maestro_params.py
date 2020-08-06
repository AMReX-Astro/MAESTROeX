#!/usr/bin/env python3

# This script parses the list of C++ runtime parameters and writes the
# necessary header files and Fortran routines to make them available
# in Maestro's C++ routines and (optionally) the Fortran routines
# through meth_params_module.
#
# parameters have the format:
#
#   name  type  default  need-in-fortran?  ifdef fortran-name  fortran-type
#
# the first three (name, type, default) are mandatory:
#
#   name: the name of the parameter.  This will be the same name as the
#     variable in C++ unless a pair is specified as (name, cpp_name)
#
#   type: the C++ data type (int, Real, bool, string)
#
#   default: the default value.  If specified as a pair, (a, b), then
#     the first value is the normal default and the second is for
#     debug mode (#ifdef AMREX_DEBUG)
#
# the next are optional:
#
#    need-in-fortran: if "y" then we do a pp.query() in meth_params.F90
#
#    ifdef: only define this parameter if the name provided is #ifdef-ed
#
#    fortran-name: if a different variable name in Fortran, specify here
#
#    fortran-type: if a different data type in Fortran, specify here
#
# Any line beginning with a "#" is ignored
#
# Commands begin with a "@":
#
#    @namespace: sets the namespace that these will be under (see below)
#      it also gives the C++ class name.
#
#      e.g. @namespace maestro Maestro
#
# Note: categories listed in the input file aren't used for code generation
# but are used for the documentation generation
#
#
# For a namespace, name, we write out:
#
#   -- name_params.H  (for maestro, included in Maestro.H):
#      sets up the namespace and extern parameters
#
#   -- name_declares.H  (for maestro, included in Maestro.cpp):
#      declares the runtime parameters
#
#   -- name_queries.H  (for maestro, included in Maestro.cpp):
#      does the parmparse query to override the default in C++
#
#  -- name_job_info_tests.H
#     this tests the current value against the default and outputs
#     into a file
#
# we write out a single copy of:
#
#   -- meth_params.F90
#      does the parmparse query to override the default in Fortran,
#      and sets a number of other parameters specific to the F90 routinse
#

import argparse
import re
import sys

FWARNING = """
! This file is automatically created by parse_maestro_params.py.  To update
! or add runtime parameters, please edit _cpp_parameters and then run
! mk_params.sh\n
"""

CWARNING = """
// This file is automatically created by parse_maestro_params.py.  To update
// or add runtime parameters, please edit _cpp_parameters and then run
// mk_params.sh\n
"""

param_include_dir = "../param_includes/"


class Param(object):
    """ the basic parameter class.  For each parameter, we hold the name,
        type, and default.  For some parameters, we also take a second
        value of the default, for use in debug mode (delimited via
        #ifdef AMREX_DEBUG)

    """

    def __init__(self, name, dtype, default,
                 cpp_var_name=None,
                 namespace=None, cpp_class=None,
                 debug_default=None,
                 in_fortran=0, f90_name=None, f90_dtype=None,
                 ifdef=None):

        self.name = name
        self.dtype = dtype
        self.default = default
        self.cpp_var_name = cpp_var_name

        self.namespace = namespace
        self.cpp_class = cpp_class

        self.debug_default = debug_default
        self.in_fortran = in_fortran

        if ifdef == "None":
            self.ifdef = None
        else:
            self.ifdef = ifdef

        if f90_name is None:
            self.f90_name = name
        else:
            self.f90_name = f90_name

        if f90_dtype is None:
            self.f90_dtype = dtype
        else:
            self.f90_dtype = f90_dtype

    def get_declare_string(self):
        # this is the line that goes into maestro_declares.H included
        # into Maestro.cpp

        if self.dtype == "int":
            tstr = "AMREX_GPU_MANAGED int {}::{}".format(
                self.namespace, self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "AMREX_GPU_MANAGED amrex::Real {}::{}".format(
                self.namespace, self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "AMREX_GPU_MANAGED bool {}::{}".format(
                self.namespace, self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "std::string {}::{}".format(
                self.namespace, self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        return "{};\n".format(tstr)

    def get_default_string(self):
        # this is the line that goes into maestro_declares.H included
        # into Maestro.cpp

        ostr = ""

        if not self.ifdef is None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += "{}::{} = {};\n".format(self.namespace,
                                            self.cpp_var_name, self.debug_default)
            ostr += "#else\n"
            ostr += "{}::{} = {};\n".format(self.namespace,
                                            self.cpp_var_name, self.default)
            ostr += "#endif\n"
        else:
            ostr += "{}::{} = {};\n".format(self.namespace,
                                            self.cpp_var_name, self.default)

        if not self.ifdef is None:
            ostr += "#endif\n"

        return ostr

    def get_f90_default_string(self):
        # this is the line that goes into read_method_params()
        # to set the default value of the variable

        ostr = ""

        # convert to the double precision notation Fortran knows
        # if the parameter is already of the form "#.e###" then
        # it is easy as swapping out "e" for "d"; if it is a number
        # like 0.1 without a format specifier, then add a d0 to it
        # because the C++ will read it in that way and we want to
        # give identical results (at least to within roundoff)

        if self.debug_default is not None:
            debug_default = self.debug_default
            if self.dtype == "Real":
                if "e" in debug_default:
                    debug_default = debug_default.replace("e", "d")
                else:
                    debug_default += "d0"

        default = self.default
        if self.dtype == "Real":
            if "e" in default:
                default = default.replace("e", "d")
            else:
                default += "d0"

        if self.dtype == "bool":
            if "true" in default:
                default = default.replace("true", ".true.")
            elif "false" in default:
                default = default.replace("false", ".false.")

        name = self.f90_name

        # for a character, we need to allocate its length.  We allocate
        # to 1, and the Fortran parmparse will resize
        if self.dtype == "string":
            ostr += "    allocate(character(len=1)::{})\n".format(name)
        else:
            ostr += "    allocate({})\n".format(name)

        if not self.debug_default is None:
            ostr += "#ifdef AMREX_DEBUG\n"
            ostr += "    {} = {};\n".format(name, debug_default)
            ostr += "#else\n"
            ostr += "    {} = {};\n".format(name, default)
            ostr += "#endif\n"
        else:
            ostr += "    {} = {};\n".format(name, default)

        return ostr

    def get_query_string(self, language):
        # this is the line that queries the ParmParse object to get
        # the value of the runtime parameter from the inputs file.
        # This goes into maestro_queries.H included into Maestro.cpp

        ostr = ""
        if not self.ifdef is None:
            ostr += "#ifdef {}\n".format(self.ifdef)

        if language == "C++":
            ostr += "pp.query(\"{}\", {}::{});\n".format(self.name,
                                                         self.namespace, self.cpp_var_name)
        elif language == "F90":
            ostr += "    call pp%query(\"{}\", {})\n".format(self.name,
                                                             self.f90_name)
        else:
            sys.exit("invalid language choice in get_query_string")

        if not self.ifdef is None:
            ostr += "#endif\n".format(self.ifdef)

        return ostr

    def default_format(self):
        """return the variable in a format that it can be recognized in C++ code"""
        if self.dtype == "string":
            return f'{self.default}'

        return self.default

    def get_job_info_test(self):
        # this is the output in C++ in the job_info writing

        ostr = f'jobInfoFile << ({self.namespace}::{self.cpp_var_name} == {self.default_format()} ? "    " : "[*] ") << "{self.namespace}.{self.cpp_var_name} = " << {self.namespace}::{self.cpp_var_name} << std::endl;\n'

        return ostr

    def get_decl_string(self):
        # this is the line that goes into maestro_params.H included
        # into Maestro.H

        if self.dtype == "int":
            tstr = "extern AMREX_GPU_MANAGED int {};\n".format(
                self.cpp_var_name)
        elif self.dtype == "Real":
            tstr = "extern AMREX_GPU_MANAGED amrex::Real {};\n".format(
                self.cpp_var_name)
        elif self.dtype == "bool":
            tstr = "extern AMREX_GPU_MANAGED bool {};\n".format(
                self.cpp_var_name)
        elif self.dtype == "string":
            tstr = "extern std::string {};\n".format(self.cpp_var_name)
        else:
            sys.exit("invalid data type for parameter {}".format(self.name))

        ostr = ""

        if not self.ifdef is None:
            ostr = "#ifdef {}\n".format(self.ifdef)

        ostr += tstr

        if not self.ifdef is None:
            ostr += "#endif\n"

        return ostr

    def get_f90_decl_string(self):
        # this is the line that goes into meth_params.f90

        if not self.in_fortran:
            return None

        if self.f90_dtype == "int":
            tstr = "integer          , allocatable, save :: {}\n".format(
                self.f90_name)
        elif self.f90_dtype == "Real":
            tstr = "double precision , allocatable, save :: {}\n".format(
                self.f90_name)
        elif self.f90_dtype == "bool":
            tstr = "logical          , allocatable, save :: {}\n".format(
                self.f90_name)
        elif self.f90_dtype == "string":
            tstr = "character (len=:), allocatable, save :: {}\n".format(
                self.f90_name)
            print("warning: string parameter {} will not be available on the GPU".format(
                self.f90_name))
        else:
            sys.exit("unsupported datatype for Fortran: {}".format(self.name))

        return tstr

    def get_cuda_managed_string(self):
        """this is the string that sets the variable as managed for CUDA"""
        if self.f90_dtype == "string":
            return "\n"
        else:
            cstr = ""
            if self.ifdef is not None:
                cstr += "#ifdef {}\n".format(self.ifdef)
            cstr += "  attributes(managed) :: {}\n".format(self.f90_name)
            if self.ifdef is not None:
                cstr += "#endif\n"
            return cstr


def write_meth_module(plist, meth_template, out_directory):
    """this writes the meth_params_module, starting with the meth_template
       and inserting the runtime parameter declaration in the correct
       place
    """

    try:
        mt = open(meth_template, "r")
    except:
        sys.exit("invalid template file")

    try:
        mo = open(f"{out_directory}/meth_params.F90", "w")
    except:
        sys.exit("unable to open meth_params.F90 for writing")

    mo.write(FWARNING)

    param_decls = [p.get_f90_decl_string() for p in plist if p.in_fortran == 1]
    params = [p for p in plist if p.in_fortran == 1]

    decls = ""

    for p in param_decls:
        decls += f"  {p}"

    cuda_managed_decls = [p.get_cuda_managed_string()
                          for p in plist if p.in_fortran == 1]

    cuda_managed_string = ""
    for p in cuda_managed_decls:
        cuda_managed_string += f"{p}"

    for line in mt:
        if line.find("@@f90_declarations@@") > 0:
            mo.write(decls)

            # Do the CUDA managed declarations

            mo.write("\n")
            mo.write("#ifdef AMREX_USE_CUDA\n")
            mo.write(cuda_managed_string)
            mo.write("#endif\n")

            # Now do the OpenACC declarations

#            mo.write("\n")
#            mo.write("  !$acc declare &\n")
#            mo.write("  !$acc create(")
#
#            for n, p in enumerate(params):
#                if p.f90_dtype == "string":
#                    print("warning: string parameter {} will not be available on the GPU".format(p.name),
#                          file=sys.stderr)
#                    continue
#
#                mo.write("{}".format(p.f90_name))
#
#                if n == len(params)-1:
#                    mo.write(")\n")
#                else:
#                    if n % 3 == 2:
#                        mo.write(") &\n  !$acc create(")
#                    else:
#                        mo.write(", ")

        elif line.find("@@set_maestro_params@@") >= 0:

            namespaces = list(set([q.namespace for q in params]))
            print("namespaces: ", namespaces)
            for nm in namespaces:
                params_nm = [q for q in params if q.namespace == nm]

                for p in params_nm:
                    mo.write(p.get_f90_default_string())

                mo.write("\n")

                mo.write('    call amrex_parmparse_build(pp, "{}")\n'.format(nm))

                for p in params_nm:
                    mo.write(p.get_query_string("F90"))

                mo.write('    call amrex_parmparse_destroy(pp)\n')

                mo.write("\n\n")

            # Now do the OpenACC device updates

#            mo.write("\n")
#            mo.write("    !$acc update &\n")
#            mo.write("    !$acc device(")
#
#            for n, p in enumerate(params):
#                if p.f90_dtype == "string": continue
#                mo.write("{}".format(p.f90_name))
#
#                if n == len(params)-1:
#                    mo.write(")\n")
#                else:
#                    if n % 3 == 2:
#                        mo.write(") &\n    !$acc device(")
#                    else:
#                        mo.write(", ")

        elif line.find("@@free_maestro_params@@") >= 0:

            params_free = [q for q in params if q.in_fortran == 1]

            for p in params_free:
                mo.write("    if (allocated({})) then\n".format(p.f90_name))
                mo.write("        deallocate({})\n".format(p.f90_name))
                mo.write("    end if\n")

            mo.write("\n\n")

        else:
            mo.write(line)

    mo.close()
    mt.close()


def parse_params(infile, meth_template, out_directory):

    params = []

    namespace = None
    cpp_class = None

    try:
        f = open(infile)
    except:
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
                cpp_class = fields[1]

            else:
                sys.exit("invalid command")

            continue

        # this splits the line into separate fields.  A field is a
        # single word or a pair in parentheses like "(a, b)"
        fields = re.findall(
            r'".+"|[\w\"\+\.-]+|\([\w+\.-]+\s*,\s*[\w\+\.-]+\)', line)

        name = fields[0]
        if name[0] == "(":
            name, cpp_var_name = re.findall(r"\w+", name)
        else:
            cpp_var_name = name

        dtype = fields[1]

        default = fields[2]
        if default[0] == "(":
            default, debug_default = re.findall(r"\w+", default)
        else:
            debug_default = None

        try:
            in_fortran_string = fields[3]
        except:
            in_fortran = 0
        else:
            if in_fortran_string.lower().strip() == "y":
                in_fortran = 1
            else:
                in_fortran = 0

        try:
            ifdef = fields[4]
        except:
            ifdef = None

        try:
            f90_name = fields[5]
        except:
            f90_name = None

        try:
            f90_dtype = fields[6]
        except:
            f90_dtype = None

        if namespace is None:
            sys.exit("namespace not set")

        params.append(Param(name, dtype, default,
                            cpp_var_name=cpp_var_name,
                            namespace=namespace,
                            cpp_class=cpp_class,
                            debug_default=debug_default,
                            in_fortran=in_fortran, f90_name=f90_name, f90_dtype=f90_dtype,
                            ifdef=ifdef))

    # output

    # find all the namespaces
    namespaces = list(set([q.namespace for q in params]))

    for nm in namespaces:

        params_nm = [q for q in params if q.namespace == nm]
        ifdefs = {q.ifdef for q in params_nm}

        # write name_declares.H
        try:
            cd = open(f"{out_directory}/{nm}_declares.H", "w")
        except:
            sys.exit(f"unable to open {nm}_declares.H for writing")

        cd.write(CWARNING)
        cd.write(f"#ifndef _{nm.upper()}_DECLARES_H_\n")
        cd.write(f"#define _{nm.upper()}_DECLARES_H_\n")

        for ifdef in ifdefs:
            if ifdef is None:
                for p in [q for q in params_nm if q.ifdef is None]:
                    cd.write(p.get_declare_string())
            else:
                cd.write("#ifdef {}\n".format(ifdef))
                for p in [q for q in params_nm if q.ifdef == ifdef]:
                    cd.write(p.get_declare_string())
                cd.write("#endif\n")

        cd.write("#endif\n")
        cd.close()

        # write name_params.H
        try:
            cp = open(f"{out_directory}/{nm}_params.H", "w")
        except:
            sys.exit(f"unable to open {nm}_params.H for writing")

        cp.write(CWARNING)
        cp.write(f"#ifndef _{nm.upper()}_PARAMS_H_\n")
        cp.write(f"#define _{nm.upper()}_PARAMS_H_\n")

        cp.write("\n")
        cp.write("namespace {} {{\n".format(nm))

        for p in params_nm:
            cp.write(p.get_decl_string())

        cp.write("};\n\n")
        cp.write("#endif\n")
        cp.close()

        # write maestro_queries.H
        try:
            cq = open(f"{out_directory}/{nm}_queries.H", "w")
        except:
            sys.exit(f"unable to open {nm}_queries.H for writing")

        cq.write(CWARNING)

        for p in params_nm:
            cq.write(p.get_default_string())
            cq.write(p.get_query_string("C++"))
            cq.write("\n")

        cq.close()

        # write the job info tests
        try:
            jo = open(f"{out_directory}/{nm}_job_info_tests.H", "w")
        except IOError:
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

    # write the Fortran module
    write_meth_module(params, meth_template, out_directory)


def main():
    """the main driver"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", type=str, default=None,
                        help="template for the meth_params module")
    parser.add_argument("-o", type=str, default=None,
                        help="output directory for the generated files")
    parser.add_argument("input_file", type=str, nargs=1,
                        help="input file containing the list of parameters we will define")

    args = parser.parse_args()

    parse_params(args.input_file[0], args.m, args.o)


if __name__ == "__main__":
    main()
