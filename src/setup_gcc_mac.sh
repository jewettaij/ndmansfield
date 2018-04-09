#!/bin/sh

export ANSI_C="gcc"
export ANSI_CPP="g++"
export L_COMP="ar rs"

export LFLAGS=""
# Note: Static linking where avaiable is recommended. If you're using linux use:
# export LFLAGS="-static"          

export MY_FLAGS="-std=c++11 -O3 -DNDEBUG -ffast-math -finline-functions -DSTRTOLD_UNSUPPORTED"
export CFLAGS="-c $MY_FLAGS"
export CPP_PRELINKER_COMMAND="echo"
export COMPILER_TEMP_FILES=""
export LINKER_TEMP_FILES=""
