#!/bin/sh

# DEBUGGIN
set -x

# get the script's directory
# (snippet courtesy stackoverflow: http://bit.ly/MmvSIz)
DIR=$( cd "$( dirname "$0" )" && pwd )
cd $DIR

# start the killing
rm ./CIL/ext/librastersystem.so

rm -rf ./build
rm -rf ./cython_imaging.egg-info
rm -rf ./dist


# TODO: leave no pyc file behind
