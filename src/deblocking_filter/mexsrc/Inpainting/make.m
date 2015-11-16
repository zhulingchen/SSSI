close all;
clear;
clc;

mex CFLAGS='\$CFLAGS -std=c99' inpainting_order_zero.c;
mex CFLAGS='\$CFLAGS -std=c99' calculate_prediction_cost.c;
mex CFLAGS='\$CFLAGS -std=c99' create_image_alignment_template_c.c;
mex CFLAGS='\$CFLAGS -std=c99' template_matching.c;
mex CFLAGS='\$CFLAGS -std=c99' template_matching_single_feature_point.c;