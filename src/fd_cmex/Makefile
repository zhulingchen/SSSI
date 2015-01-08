#!/bin/bash

# Determine system architecture
UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX = mex
endif

ifeq ($(UNAME), Darwin)
	MEX = /Applications/MATLAB_R2014a.app/bin/mex
endif

# Default compiler


# Compiler Flags
# CFLAG = -lm

# Files need to be compiled
FD_SRC = fwdTimeCpmlFor2dAw_mex.c rvsTimeCpmlFor2dAw_mex.c

# Lib file
LIB = finiteDifference.o

all: fd

fd: finiteDifference.o
	$(MEX) fwdTimeCpmlFor2dAw_mex.c ${LIB}
	$(MEX) rvsTimeCpmlFor2dAw_mex.c ${LIB}

finiteDifference.o: finiteDifference.c finiteDifference.h
	$(MEX) -c finiteDifference.c

