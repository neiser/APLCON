#!/bin/bash

set -e

CACHE=$HOME/cache
CACHE_REV=$CACHE/.rev1
NCPU=2

pushd $PWD

if [ -f $CACHE_REV ]; then
    export PATH=$CACHE/cmake/bin:$PATH
    export PATH=$CACHE/gcc/bin:$PATH
    export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
    export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH    
else
    # Clean the cache
    echo "Cleaning the cache..."
    rm -rf $CACHE
    mkdir -p $CACHE
    
    # Get CMake 3.1
    wget https://github.com/Viq111/travis-container-packets/releases/download/cmake-3.1.2/cmake.tar.bz2 -O $CACHE/cmake.tar.bz2
    tar -xjf $CACHE/cmake.tar.bz2 -C $CACHE
    rm $CACHE/cmake.tar.bz2
    export PATH=$CACHE/cmake/bin:$PATH
    
    # Get GCC 4.8
    wget https://github.com/Viq111/travis-container-packets/releases/download/gcc-4.8.2/gcc.tar.bz2 -O $CACHE/gcc.tar.bz2
    tar -xjf $CACHE/gcc.tar.bz2 -C $CACHE
    rm $CACHE/gcc.tar.bz2
    export PATH=$CACHE/gcc/bin:$PATH
    export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
    export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
    export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH

    touch $CACHE_REV
fi

popd

mkdir build && cd build
cmake -DCTEST_PARALLEL_JOBS=$NCPU -DCOVERAGE=On ..
make -j$NCPU
make -j$NCPU build_and_test
