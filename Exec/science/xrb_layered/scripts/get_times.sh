#!/bin/bash
RUN_DIR=$1
PLOT_DIR=PLOTS
PLOT_BASENAME="xrb_0"

if [ "${RUN_DIR: -1}" == "/" ]; then
    RUN_DIR=${RUN_DIR:0:-1}
fi
#echo $RUN_DIR

ftime=$AMREX_HOME/Tools/Plotfile/ftime.gnu.x86-milan.ex
$ftime $RUN_DIR/$PLOT_DIR/$PLOT_BASENAME* > $RUN_DIR/"times.txt"
