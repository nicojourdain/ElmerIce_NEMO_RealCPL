#!/bin/sh -f
#the compiler wrapper scripts
export CC=mpicc
export CXX=mpic++
export FC=mpif90
export F77=mpif90

# Version of Mac OS X 
export MACOSX_DEPLOYMENT_TARGET=10.6

#the compiler flags
#export CFLAGS="-I/opt/hypre/include -O5 -ftree-vectorize"
#export CXXFLAGS="-I/opt/hypre/include -O5 -ftree-vectorize"
#export FCFLAGS="-I/opt/hypre/include -O5 -m64 -ftree-vectorize -funroll-loops"
#export F77FLAGS="-I/opt/hypre/include -O5 -m64 -ftree-vectorize -funroll-loops"
#export FFLAGS="-I/opt/hypre/include -O5 -m64 -ftree-vectorize -funroll-loops"

export CFLAGS="-I/opt/hypre/include -I/opt/include -I/usr/X11/include -O2 -funroll-loops"
export CXXFLAGS="-I/opt/hypre/include -I/opt/include -I/usr/X11/include -O2 -funroll-loops"
export FCFLAGS="-I/opt/hypre/include -O2 -funroll-loops"
export F77FLAGS="-I/opt/hypre/include -O2 -funroll-loops"
export FFLAGS="-I/opt/hypre/include -O2 -funroll-loops"

#linking to Hypre library and Mac OS X's Blas/Lapack
export LDFLAGS="-L/opt/hypre/lib -lHYPRE -L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A -lBLAS -L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A -lLAPACK -L/opt/lib/ -L/usr/X11/lib/"

##paths to elmer and MPI
export ELMER_HOME="/Applications/elmerdev"
export MPI_HOME="/opt/openmpi"

#we want to use Mac OS X's super fast Blas and Lapack
export BLAS_HOME="/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A"
export LAPACK_HOME="/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A"

#general paths
export PATH="/opt/openmpi/bin:$PATH"
export LD_LIBRARY_PATH="/opt/openmpi/lib:/opt/hypre/lib:$LD_LIBRARY_PATH"

# modules
modules="matc umfpack mathlibs elmergrid meshgen2d eio hutiter fem post"

# configure and build
for m in $modules; do
    cd $m ; make distclean ; ./configure --prefix=$ELMER_HOME --with-mpi-dir=$MPI_HOME --with-blas=$BLAS_HOME/libBLAS.dylib --with-lapack=$LAPACK_HOME/libLAPACK.dylib LIBS="-L/opt/lib/-ltcl8.5 -ltk8.5 -L/usr/X11/lib/-lGL -lGLU -lX11" && make && make install && cd ..
done
