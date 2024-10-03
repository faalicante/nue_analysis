#!/bin/bash

BRICK=$3

SNDBUILD_DIR=/afs/cern.ch/work/f/falicant/public/SNDBUILD/sw
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
source /afs/cern.ch/work/f/falicant/public/fedra/setup_new.sh	

root -l -b -q /eos/user/f/falicant/Simulations_sndlhc/nuecc_withcrisfiles_25_July_2022/b0000$BRICK/hist_XYP_nue.root /afs/cern.ch/work/f/falicant/public/nue_analysis/tag_basetrack_mc.C\($BRICK\)

