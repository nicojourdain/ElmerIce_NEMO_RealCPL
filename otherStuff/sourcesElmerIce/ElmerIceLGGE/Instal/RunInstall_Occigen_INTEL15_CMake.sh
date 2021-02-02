#!/bin/bash

#------------------------------------------------#
#	Script for Elmer on Occigen [Dec2015]
#	       INTEL 15 compilers
#Module to load:
#module load intel/15.0.3.187
#module load intelmpi
#module load netcdf
#module load cmake (only for installation > 2.8.9)
#------------------------------------------------#

CMAKE=cmake

# Installation directory (set these!)
ELMERDIR="/home/ltavard/MyElmer/Elmer_RELEASE"
ELMERSRC="$ELMERDIR/elmer_release"
BUILDDIR="$ELMERDIR/build"
# Set ELMER_HOME in your ~/.bashrc !
IDIR=$ELMER_HOME

#Set path to mumps (optional)
MUMPS_DIR="/home/ltavard/CODES/MUMPS_5.0.0"
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
    -DCMAKE_C_COMPILER=/opt/software/intel/composer_xe_2015.3.187/bin/intel64/icc \
    -DCMAKE_CXX_COMPILER=/opt/software/intel/composer_xe_2015.3.187/bin/intel64/icpc\
    -DCMAKE_Fortran_COMPILER=/opt/software/intel/composer_xe_2015.3.187/bin/intel64/ifort\
    -DCMAKE_Fortran_FLAGS_DEBUG=-g\
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
    -DMPI_CXX_COMPILER=/opt/software/intel/impi/5.0.1.035/bin64/mpiicc \
    -DMPI_C_COMPILER=/opt/software/intel/impi/5.0.1.035/bin64/mpiicpc \
    -DMPI_Fortran_COMPILER=/opt/software/intel/impi/5.0.1.035/bin64/mpiifort \
    -DMUMPS_D_LIB=$MUMPS_DIR/lib/libdmumps.a \
    -DMUMPS_COMMON_LIB=$MUMPS_DIR/lib/libmumps_common.a \
    -DMUMPS_PORD_LIB=$MUMPS_DIR/lib/libpord.a \
    -DMumps_INCLUDE_DIR=$MUMPS_DIR/include/ \

echo "------------------------------------------"
echo "-- Compile & install with"
echo "-- make -j4 && make install "
echo "-- change the -j4 to the number of available cores on your system"
echo "------------------------------------------"

