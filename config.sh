#!/bin/bash

pwdPATH=$(pwd)
flib=fftw_comp 
fftw=fftw-3.3.9

mkdir fftw
cp $fftw.tar.gz  fftw/
cd fftw
tar -zxvf $fftw.tar.gz

mkdir $flib

cd $fftw
./configure --prefix=$pwdPATH/fftw/$flib/
make CFLAGS="-fPIC"
make install 


