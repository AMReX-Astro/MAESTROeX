#!/bin/sh

./parse_maestro_params.py -m meth_params.template ./_cpp_parameters

./set_variables.py _variables
