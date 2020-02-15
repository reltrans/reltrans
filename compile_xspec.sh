#!/bin/bash

mkdir makeshield         #make the directory to shield the correct makefile
mv Makefile makeshield/  #move the correct makefile in the shield 

echo "initpackage reltrans lmodel_reltrans.dat `pwd` \nexit" | xspec

rm -f Makefile        #remove the incorrect makefile created by initpackage
mv makeshield/Makefile . # restore the correct makefile from the shield
rm -vf *~ *.o         #prepare for the next compilation        
rm -vrf makeshield    # remove the shield 

echo "hmake \nexit" | xspec   #run the hmake in xspec with the makefile

echo " lmod reltrans ." | xspec  

#cleaning up
rm -vf *~ *.o
rm -vf *FunctionMap.* lpack_* *tcl
rm -vf *.mod
