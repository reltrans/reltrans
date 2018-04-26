#!/bin/bash

echo "initpackage reltrans lmodel_relxill.dat ." | xspec

rm -f *~ *.o
rm -f *FunctionMap.* lpack_* 
rm -f *.mod Makefile
