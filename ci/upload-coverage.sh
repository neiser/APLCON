#!/bin/bash

set -e

CACHE=$HOME/cache

export PATH=$CACHE/gcc/bin:$PATH
export LIBRARY_PATH=$CACHE/gcc/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$CACHE/gcc/lib64:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$CACHE/gcc/include/c++/4.8.2:$CPLUS_INCLUDE_PATH

bash <(curl -s https://codecov.io/bash)
