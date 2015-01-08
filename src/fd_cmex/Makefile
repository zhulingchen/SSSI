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
#MEX_FLAG = -O CFLAGS="\$CFLAGS -std=c99"

# Files need to be compiled
#FD_SRC = fwdTimeCpmlFor2dAw_mex.c rvsTimeCpmlFor2dAw_mex.c

# Lib file
LIB = finiteDifference.o

all: fd

fd: finiteDifference.o
	$(MEX) diffOperator_mex.c
	$(MEX) fwdTimeCpmlFor2dAw_mex.c ${LIB}
	$(MEX) rvsTimeCpmlFor2dAw_mex.c ${LIB}

finiteDifference.o: finiteDifference.c finiteDifference.h
	$(MEX) -c finiteDifference.c

