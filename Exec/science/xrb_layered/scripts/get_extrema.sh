#!/bin/bash
RUN_DIR=$1
PLOT_DIR=PLOTS
PLOT_BASENAME="xrb_0"

VARS="tfromp MachNumber Hnuc"

fextrema=$AMREX_HOME/Tools/Plotfile/fextrema.gnu.x86-milan.ex
$fextrema -v "$VARS" $RUN_DIR/$PLOT_DIR/$PLOT_BASENAME* > $RUN_DIR/"extrema.txt"
