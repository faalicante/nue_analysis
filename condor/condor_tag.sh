#!/bin/bash

CELL=$3

SNDBUILD_DIR=/afs/cern.ch/work/f/falicant/public/SNDBUILD/sw
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
source /afs/cern.ch/work/f/falicant/public/fedra/setup_new.sh	

root -l -b -q /eos/user/f/falicant/nue_search/R1B121/gen1/hist/hist_XYP_b121_$CELL.root /afs/cern.ch/work/f/falicant/public/nue_analysis/tag_basetrack.C\($CELL\)

