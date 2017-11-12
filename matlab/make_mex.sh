#!/bin/bash
# Generate MEX files by compiling from Matlab
# This scripts intends to automate the compilation process
# by following the instructions provided in
# http://laurentkneip.github.io/opengv/page_installation.html#sec_installation_5
# We assume the installation configuration is standard
# so that every dependency can be found in an automated way

# We assume a launcher command is available: matlab
# Add OpenGV library directory to the path
export LD_LIBRARY_PATH=../build/lib:$LD_LIBRARY_PATH

# Find path to Eigen library (assumes CMake has cached EIGEN_INCLUDE_DIRS)
# See https://stackoverflow.com/questions/8474753/how-to-get-a-cmake-variable-from-the-command-line
EigenPath=$(cmake -L ../build | grep EIGEN_INCLUDE_DIRS | cut -d "=" -f2)

# Call Matlab with the compilation command
matlab -nodisplay -nosplash -nodesktop \
-r "mex -I../include -I${EigenPath} -L../build/lib -lopengv opengv.cpp -cxx, mex -I../include -I${EigenPath} -L../build/lib -lopengv opengv_donotuse.cpp -cxx, exit"


