#!/bin/bash

pwdPATH=$(pwd)
flib=fftw_comp 
fftw=fftw-3.3.8

cd fftw
tar -zxvf $fftw.tar.gz

mkdir $flib

cd $fftw
./configure --prefix=$pwdPATH/fftw/$flib/
make
make install 


