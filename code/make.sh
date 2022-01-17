#!/bin/bash
rm -r -f lib
mkdir "lib"
gcc "classicFW.c" -o "lib/classicFW"
mpicc "blockedFW.c" -o "lib/blocked"
mpic++ "improvement.cpp" -o "lib/improvement"  
g++ "graph_gen.cpp" -o "lib/graph_gen"
mpicc "2dsparseFW.c" -o "lib/2dsparseFW" -lm 
cp "test.py" "lib/test.py"
