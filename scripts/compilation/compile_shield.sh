#!/bin/bash
# Script to compile SHiELD

if [ "$#" -ne 3 ]; then
  echo "Usage: ./compile_shield.sh COMP"
  echo "choices for COMP :  debug, repro, prod"
  echo "choices for clean:  noclean, clean, cleanall"
  echo "choices for mode :  32bit, 64bit"
  exit 1
fi
COMP="$1"
clean="$2"
mode="$3"
compiler="intel"

model=nh
cd ../../SHiELD_build/Build
./COMPILE shield $COMP $compiler $mode $clean

# Copy the log file to test directory
cp -r build_shield_$model.$COMP.$mode.$compiler.out ../../scripts/compilation/build_shield_$model.$COMP.$mode.$compiler.out
