close all;
clear;
clc;

mex CFLAGS='\$CFLAGS -std=c99' estimate_steering_parameters.c 