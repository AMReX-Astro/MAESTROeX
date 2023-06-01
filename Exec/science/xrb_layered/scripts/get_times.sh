#!/bin/bash
RUN_DIR=$1
PLOT_DIR=PLOTS
PLOT_BASENAME="xrb_0"

ftime=$AMREX_HOME/Tools/Plotfile/ftime.gnu.x86-milan.ex
$ftime $RUN_DIR/$PLOT_DIR/$PLOT_BASENAME* > $RUN_DIR/"times.txt"
