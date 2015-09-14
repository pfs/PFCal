#!/bin/bash

source g4env.sh 

echo "Making userlib"
cd userlib/
mkdir obj
mkdir lib
mkdir bin
make dictionary
cd -

echo "Making main"
make

echo "Making analysis"
cd analysis/
mkdir obj
mkdir lib
mkdir bin
make
cd -

echo "...all done"

