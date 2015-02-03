#!/bin/bash

# Determine system architecture and determine compiler
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX = mex
endif

ifeq ($(UNAME), Darwin)
	MEX = /Applications/MATLAB_R2014a.app/bin/mex
endif

# Compiler Flags
#MEX_FLAG_REGULAR = CFLAGS="\$CFLAGS -std=c99"
#-Ofast is available only after gcc4.6
MEX_FLAG_REGULAR = COPTIMFLAGS="-O3" CXXOPTIMFLAGS="-O3"


# Files need to be compiled
#FD_SRC = fwdTimeCpmlFor2dAw_mex.c rvsTimeCpmlFor2dAw_mex.c

# Lib file
LIB = finiteDifference.o

all: fd

fd: finiteDifference.o
	$(MEX) $(MEX_FLAG_REGULAR) diffOperator_mex.c
	$(MEX) $(MEX_FLAG_REGULAR) fwdTimeCpmlFor2dAw_mex.c ${LIB}
	$(MEX) $(MEX_FLAG_REGULAR) rvsTimeCpmlFor2dAw_mex.c ${LIB}


finiteDifference.o: finiteDifference.c finiteDifference.h
	$(MEX) -c $(MEX_FLAG_REGULAR) finiteDifference.c

