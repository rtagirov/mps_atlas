#!/bin/bash
export CXX=''
rm -rf build lib
mkdir build
mkdir lib
cd build
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/old/netcdf-3.6.1.tar.gz
tar xfz netcdf-3.6.1.tar.gz
cd netcdf-3.6.1/src/
./configure --disable-cxx --prefix="$PWD"/../../../lib
make
make check
make install
