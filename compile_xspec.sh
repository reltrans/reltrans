#!/bin/bash

echo "initpackage crossen lmodel.dat ." | xspec

rm -f *~ *.o
rm -f *FunctionMap.* lpack_* 
rm -f *.mod Makefile
