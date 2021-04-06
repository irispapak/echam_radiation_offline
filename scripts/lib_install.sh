#!/bin/bash 
#-ev
#
# Just run this script from the directory you want ECHAM libraries to be installed, sit back and wait :)
# if you have any questions ask Marc 
linker(){
cd $1/lib
if [[ -d ../lib64 ]]; then
 infiles=../lib64/*
 for file in ${infiles[*]}; do
  ln -s $file
 done
fi
}
addlib(){
if [[ x$LD_LIBRARY_PATH != x ]]; then
 export LD_LIBRARY_PATH=${1}/lib:$LD_LIBRARY_PATH
else
 export LD_LIBRARY_PATH=${1}/lib
fi
}

# This has to be changed
opt=

export SW_DIR=${opt}/src

set -ex
which gcc
which gfortran

which gcc
which gfortran

export FC=gcc
export F77=gfortran
export F95=gfortran
export zlib_VERSION=1.2.11
export szip_VERSION=2.1.1
export hdf5_VERSION=1.10.6
export netcdf4_VERSION=c-4.7.4
export netcdf_fortran_VERSION=4.5.3
export mpich_VERSION=3.3


mkdir -p ${SW_DIR}/source

# zlib
cd ${SW_DIR}/source
if [ ! -e zlib-${zlib_VERSION}.tar.gz ]; then
  wget http://www.zlib.net/zlib-${zlib_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/zlib-${zlib_VERSION} ]; then
  tar xzf zlib-${zlib_VERSION}.tar.gz
fi
dir_zlib=${opt}/zlib/${zlib_VERSION}
  echo processing ${zlib_VERSION}
  cd ${SW_DIR}/source/zlib-${zlib_VERSION}
  ./configure --prefix=${dir_zlib}
  make
  make check install
  rm -f ${opt}/zlib/current
  ln -s ${dir_zlib} ${opt}/zlib/current
addlib ${opt}/zlib/current


#szip
cd ${SW_DIR}/source
if [ ! -e szip-${szip_VERSION}.tar.gz ]; then
  wget http://www.hdfgroup.org/ftp/lib-external/szip/${szip_VERSION}/src/szip-${szip_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/szip-${szip_VERSION} ]; then
  tar xzf szip-${szip_VERSION}.tar.gz
fi
dir_szip=${opt}/szip/${szip_VERSION}
echo processing szip-${szip_VERSION}
cd szip-${szip_VERSION}
./configure --prefix=${dir_szip}
make
make check install
#linker ${dir_szip}
rm -f ${opt}/szip/current
ln -s ${dir_szip} ${opt}/szip/current
addlib ${opt}/szip/current


# hdf5
cd ${SW_DIR}/source
if [ ! -e hdf5-${hdf5_VERSION}.tar.gz ]; then
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.6/src/hdf5-${hdf5_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/hdf5-${hdf5_VERSION} ]; then
  tar xzf hdf5-${hdf5_VERSION}.tar.gz
fi
dir_hdf5=${opt}/hdf5/${hdf5_VERSION}
echo processing hdf5-${hdf5_VERSION}
cd hdf5-${hdf5_VERSION}
./configure --prefix=${dir_hdf5} --with-szlib=${dir_szip}/lib --with-zlib=${dir_zlib}/include,${dir_zlib}/lib
make -j 4
make install
#linker ${dir_hdf5}
rm -f ${opt}/hdf5/current
ln -s ${dir_hdf5} ${opt}/hdf5/current
addlib ${opt}/hdf5/current


# netcdf4
cd ${SW_DIR}/source
if [ ! -e netcdf-${netcdf4_VERSION}.tar.gz ]; then
  wget http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-${netcdf4_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/netcdf-${netcdf4_VERSION} ]; then
  tar xzf netcdf-${netcdf4_VERSION}.tar.gz
fi
dir_netcdf=${opt}/netcdf/${netcdf4_VERSION}_${F77}
echo processing netcdf-${netcdf4_VERSION}
cd netcdf-${netcdf4_VERSION}
CPPFLAGS="-I${dir_zlib}/include -I${dir_szip}/include -I${dir_hdf5}/include"  LDFLAGS="-L${dir_zlib}/lib -L${dir_szip}/lib -L${dir_hdf5}/lib" ./configure --prefix=${dir_netcdf} --disable-dap
make -j 4
make check install
rm -f ${opt}/netcdf/current
ln -s ${dir_netcdf} ${opt}/netcdf/current
addlib ${opt}/netcdf/current


# netcdf4_fortran
cd ${SW_DIR}/source
if [ ! -e netcdf-fortran-${netcdf_fortran_VERSION}.tar.gz ]; then
  wget http://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-fortran-${netcdf_fortran_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/netcdf-fortran-${netcdf_fortran_VERSION} ]; then
  tar xzf netcdf-fortran-${netcdf_fortran_VERSION}.tar.gz
fi
dir_netcdf_fortran=${opt}/netcdf-fortran/${netcdf_fortran_VERSION}
echo processing netcdf-fortran-${netcdf_fortran_VERSION}
cd netcdf-fortran-${netcdf_fortran_VERSION}
export FC=gfortran 
export FFLAGS='-O3' 
CPPFLAGS="-I${dir_zlib}/include -I${dir_szip}/include -I${dir_hdf5}/include -I${dir_netcdf}/include" LD_LIBRARY_PATH=${dir_netcdf}/lib:${LD_LIBRARY_PATH}  LDFLAGS="-L${dir_zlib}/lib -L${dir_szip}/lib -L${dir_hdf5}/lib -L${dir_netcdf}/lib"  ./configure --prefix=${dir_netcdf} 
make -j 4
make install


#mpich
cd ${SW_DIR}/source
if [ ! -e mpich-${mpich_VERSION}.tar.gz ]; then
  wget http://www.mpich.org/static/downloads/${mpich_VERSION}/mpich-${mpich_VERSION}.tar.gz
fi
if [ ! -d ${SW_DIR}/source/mpich-${mpich_VERSION} ];then
  tar xzf mpich-${mpich_VERSION}.tar.gz
fi
dir_mpich=${opt}/mpich/${mpich_VERSION}_${F77}
echo processing mpich-${mpich_VERSION}
export CC=gcc
export CXX=gcc
export FC=gfortran
export F77=gfortran
cd mpich-${mpich_VERSION}
./configure --prefix=${dir_mpich} 
make -j 4
make install
rm -f ${opt}/mpich/current
ln -s ${dir_mpich} ${opt}/mpich/current
addlib ${opt}/mpich/current

