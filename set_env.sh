#!/bin/bash
# Script that sets up enviornment on UCSB servers

# Environment variable for different Linux versions
export RUN_KERNEL=$(uname -r | cut -d '-' -f1)

# Setting up ROOT using cmssw at ucsb servers. 
# Please setup ROOT yourself if you are not working on a ucsb server.
. /cvmfs/cms.cern.ch/cmsset_default.sh
if [ "$RUN_KERNEL" == "3.10.0" ]; then
  export SCRAM_ARCH=slc7_amd64_gcc700
  cd /net/cms11/data/pico/cc7/CMSSW_10_2_13/src
elif [ "$RUN_KERNEL" == "2.6.32" ]; then
  cd /net/cms29/cms29r0/pico/CMSSW_10_2_11_patch1/src
fi
. /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -

## Scons setup
# Multicore build
export SCONSFLAGS="-j $(nproc --all)"
export SET_ENV_PATH=set_env.sh # environment to use for build
# For shared built libraries
export LD_LIBRARY_PATH=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/lib/$RUN_KERNEL${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
#for python bindings
#ldconfig -n $(dirname $(readlink -e "$BASH_SOURCE"))/lib
#export LD_LIBRARY_PATH=$(dirname $(readlink -e "$BASH_SOURCE"))/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$(dirname $(readlink -e "$BASH_SOURCE"))/bindings:$PYTHONPATH

# Setup job environment, sometimes used for running combine
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/jb_utils/set_env.sh
source $(dirname $(readlink -e "$BASH_SOURCE"))/modules/queue_system/set_env.sh

# Setup for automatically using different folders in scripts for certain users.
usernames=("jbkim oshiro")
if [[ "${usernames[@]}" =~ "$USER" ]]; then
  if [ "$HOSTNAME" = cms37.physics.ucsb.edu ]; then
    export LOCAL_PICO_DIR="/data/"
    echo "LOCAL_PICO_DIR is set to $LOCAL_PICO_DIR. Set it so blank if you want to use net files."
  fi
else
  export LOCAL_PICO_DIR=""
fi
