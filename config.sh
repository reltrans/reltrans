#!/bin/bash

pwdPATH=$(pwd)
flib=fftw_comp 
fftw=fftw-3.3.10

mkdir fftw
cp $fftw.tar  fftw/
cd fftw
tar -zxvf $fftw.tar

mkdir $flib

cd $fftw
./configure --prefix=$pwdPATH/fftw/$flib/
make CFLAGS="-fPIC"
make install 


