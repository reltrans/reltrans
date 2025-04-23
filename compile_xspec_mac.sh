#!/bin/bash

echo "initpackage reltrans lmodel_reltrans.dat `pwd` \nexit" | xspec

#Include in the failing Makefile created by xspec the lib that it needs to compile the fftw
#For Mac OS
sed -i '' '1s/^/libs = -L fftw\/fftw_comp\/lib\/ -lfftw3 -lm\'$'\n/' Makefile
sed -i '' '1s/^/incs = -I fftw\/fftw_comp\/include\/ \'$'\n/' Makefile
sed -i '' '1s/^/optimization = -O3 \'$'\n/' Makefile
sed -i '' 's/HD_FFLAGS		=/HD_FFLAGS = ${optimization} ${incs}/g' Makefile
sed -i '' 's/HD_SHLIB_LIBS           =/HD_SHLIB_LIBS = ${optimization} ${libs}/g' Makefile

#For Linux OS (it needs to be tested!!)
# sed -i  '1s/^/libs = -L fftw\/fftw_comp\/lib\/ -lfftw3 -lm \n/' Makefile
# sed -i  '1s/^/incs = -I fftw\/fftw_comp\/include\/ \n/' Makefile
# sed -i  '1s/^/optimization = -O3 \n/' Makefile
# sed -i  's/HD_FFLAGS		=/HD_FFLAGS = ${optimization} ${incs}/g' Makefile
# sed -i  's/HD_SHLIB_LIBS           =/HD_SHLIB_LIBS = ${optimization} ${libs}/g' Makefile


rm -vf *~ *.o         #prepare for the next compilation        

echo "hmake \nexit" | xspec   #run the hmake in xspec with the makefile

echo " load libreltrans.dylib" | xspec  

#cleaning up
rm -vf *~ *.o
rm -vf *FunctionMap.* lpack_* *tcl
rm -vf *.mod
rm -vf Makefile


# mkdir makeshield         #make the directory to shield the correct makefile
# mv Makefile makeshield/  #move the correct makefile in the shield 
# rm -f Makefile        #remove the incorrect makefile created by initpackage
# mv makeshield/Makefile . # restore the correct makefile from the shield
# rm -vrf makeshield    # remove the shield 
