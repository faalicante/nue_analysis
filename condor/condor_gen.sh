#!/bin/bash

CELL=$3

SNDBUILD_DIR=/afs/cern.ch/work/f/falicant/public/SNDBUILD/sw
source /cvmfs/sndlhc.cern.ch/SNDLHC-2023/Aug30/setUp.sh
eval `alienv load -w $SNDBUILD_DIR --no-refresh sndsw/latest`
source /afs/cern.ch/work/f/falicant/public/fedra/setup_new.sh	

cd /eos/user/f/falicant/nue_search/R1B121/gen1

python /afs/cern.ch/work/f/falicant/public/nue_analysis/genXYP.py --cell $CELL

