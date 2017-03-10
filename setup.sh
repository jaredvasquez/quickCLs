#!/bin/bash

source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/gcc49/setup.sh
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.16-HiggsComb-x86_64-slc6-gcc49-opt/bin/thisroot.sh

which gcc
which root

# dynamic libraries
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`/lib:/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/boost/boost-1.60.0-python2.7-x86_64-slc6-gcc49/boost-1.60.0-python2.7-x86_64-slc6-gcc49/lib/

# Greet the user
echo "Setting up environment for compiling/running SFrame"

if [ $_DIRFIT ]; then
    echo _DIRFIT is already defined, use a clean shell
    return 1
fi

# speficy the SFRAME base directory, i.e. the directory in which this file lives
export _DIRFIT=${PWD}

# Modify to describe your directory structure. Default is to use the a structure where
# all directories are below the SFrame base directory specified above
export _BIN_PATH=${_DIRFIT}/bin
export _LIB_PATH=${_DIRFIT}/lib

# Check if bin/lib/include directories exist, if not create them
if [ ! -d ${_BIN_PATH} ]; then
    echo Directory ${_BIN_PATH} does not exist ... creating it
    mkdir ${_BIN_PATH}
fi
if [ ! -d ${_LIB_PATH} ]; then
    echo Directory ${_LIB_PATH} does not exist ... creating it
    mkdir ${_LIB_PATH}
fi

# The Makefiles depend only on the root-config script to use ROOT,
# so make sure that is available
if [[ `which root-config` == "" ]]; then
    echo "Error: ROOT environment doesn't seem to be configured!"
fi

if [[ `root-config --platform` == "macosx" ]]; then

    # With Fink ROOT installations, DYLD_LIBRARY_PATH doesn't have
    # to be defined for ROOT to work. So let's leave the test for it...
    export DYLD_LIBRARY_PATH=${_LIB_PATH}:${DYLD_LIBRARY_PATH}

else

    if [ ! $LD_LIBRARY_PATH ]; then
        echo "Warning: so far you haven't setup your ROOT enviroment properly (no LD_LIBRARY_PATH): SFrame will not work"
    fi

    export LD_LIBRARY_PATH=${_LIB_PATH}:${LD_LIBRARY_PATH}

fi

export PATH=${_BIN_PATH}:${PATH}
export PYTHONPATH=${_DIRFIT}/python:${PYTHONPATH}

export PAR_PATH=./:${_LIB_PATH}
