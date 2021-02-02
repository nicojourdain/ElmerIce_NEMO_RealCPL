#!/bin/bash

#------------------------------------------------#
#	Script for Elmer on Froggy [Dec2015]
#	       INTEL 15 compilers
#Module to load:
#source /applis/site/env.bash
#module load intel-devel/2015
#module load intel-mkl/11.3.0.109
#module load netcdf (auto loaded: hdf5/1.8.15-patch1_intel-15.0.3 && netcdf/4.3.3.1_intel-15.0.3)
#module load cmake (only for installation > 2.8.9)
#------------------------------------------------#

CMAKE=cmake

# Installation directory (set these!)
ELMERSRC="/home/ltavard/MyElmer/Elmer_RELEASE/elmer_release"
BUILDDIR="/home/ltavard/MyElmer/Elmer_RELEASE/build"
IDIR=/home/ltavard/MyElmer/Elmer_RELEASE/install_$ELMER_REV


echo "------------------------------------------"
echo "-- Building Elmer from within " ${BUILDDIR}
echo "-- Installation into " ${IDIR}
echo "-- CMake info"
cmake --version
echo "------------------------------------------"
cd ${BUILDDIR}
pwd
ls -ltr
echo "------------------------------------------"
echo $CMAKE $ELMERSRC 
echo "------------------------------------------"
# you can add a -DCMAKE_TOOLCHAIN_FILE=$TOOLCHAIN,
# if you have a toolchain file declared
$CMAKE -Wno-dev $ELMERSRC  \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DCMAKE_COLOR_MAKEFILE:BOOL=ON \
    -DCMAKE_INSTALL_PREFIX=$IDIR \
    -DCMAKE_C_COMPILER=/opt/intel/2015/composer_xe_2015.3.187/bin/intel64/icc\
    -DCMAKE_CXX_COMPILER=/opt/intel/2015/composer_xe_2015.3.187/bin/intel64/icpc \
    -DCMAKE_Fortran_COMPILER=/opt/intel/2015/composer_xe_2015.3.187/bin/intel64/ifort \
    -DWITH_ELMERGUI:BOOL=FALSE \
    -DWITH_ELMERGUILOGGER:BOOL=FALSE \
    -DWITH_ELMERGUITESTER:BOOL=FALSE \
    -DWITH_ElmerIce:BOOL=TRUE \
    -DWITH_Hypre:BOOL=FALSE \
    -DWITH_MKL:BOOL=TRUE \
    -DWITH_MPI:BOOL=TRUE \
    -DWITH_Mumps:BOOL=TRUE \
    -DWITH_OpenMP:BOOL=TRUE \
    -DWITH_Trilinos:BOOL=FALSE \
    -DMPI_CXX_COMPILER=/opt/intel/2015_mpi/compilers_and_libraries_2016.0.079/linux/mpi/intel64/bin/mpiicpc \
    -DMPI_C_COMPILER=/opt/intel/2015_mpi/compilers_and_libraries_2016.0.079/linux/mpi/intel64/bin/mpiicc \
    -DMPI_Fortran_COMPILER=/opt/intel/2015_mpi/compilers_and_libraries_2016.0.079/linux/mpi/intel64/bin/mpiifort \
    -DMUMPS_D_LIB=/home/ltavard/CODES/MUMPS_5.0.0/lib/libdmumps.a \
    -DMUMPS_COMMON_LIB=/home/ltavard/CODES/MUMPS_5.0.0/lib/libmumps_common.a \
    -DMUMPS_PORD_LIB=/home/ltavard/CODES/MUMPS_5.0.0/lib/libpord.a \
    -DMumps_INCLUDE_DIR=/home/ltavard/CODES/MUMPS_5.0.0/include/

echo "------------------------------------------"
echo "-- Compile & install with"
echo "-- make -j4 && make install "
echo "-- change the -j4 to the number of available cores on your system"
echo "------------------------------------------"

