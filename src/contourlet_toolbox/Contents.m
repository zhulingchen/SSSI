%Contourlet Toolbox (Version 2.0)
%
%Demos
% DECDEMO    Demonstrates contourlet decomposition and reconstruction 
% NLADEMO    Demo for contourlet nonlinear approximation.
% NLADEMO2   Nonlinear approximation demo using only the finest scale
% DENOISEDEMO  Denoise demo
%
%Main functions: Contourlets pyramidal directional filter bank
% PDFBDEC    Pyramidal Directional Filter Bank (or Contourlet) Decomposition
% PDFBREC    Pyramid Directional Filterbank Reconstruction
% 
%Retrieve filters by names  
% PFILTERS   Generate filters for the Laplacian pyramid
% DFILTERS   Generate directional 2D filters
% LDFILTER   Generate filter for the ladder structure
%
%Utility functions of the contourlet transform
% SHOWPDFB   Show contourlet or PDFB coefficients. 
% PDFB2VEC   Convert the output of the PDFB into a vector form
% VEC2PDFB   Convert the vector form to the output structure of the PDFB
% PDFB_TR    Retain the most significant coefficients at certain subbands
% PDFB_NEST  Estimate the noise standard deviation in the PDFB domain
%
%Laplacian pyramid 
% LPDEC      Laplacian Pyramid Decomposition
% LPDEC      Laplacian Pyramid Reconstruction
%
%Wavelet filter bank 
% WFB2DEC    2-D Wavelet Filter Bank Decomposition
% WFB2REC    2-D Wavelet Filter Bank Decomposition
%
%Directional filter bank 
% DFBDEC     Directional Filterbank Decomposition
% DFBREC     Directional Filterbank Reconstruction
% DFBDEC_L   Directional Filterbank Decomposition using Ladder Structure
% DFBREC_L   Directional Filterbank Reconstruction using Ladder Structure
% DFBIMAGE   Produce an image from the result subbands of DFB
% 
%Two-channel 2D filter banks (used in the DFB)  
% FBDEC      Two-channel 2D Filterbank Decomposition
% FBDEC_L    Two-channel 2D Filterbank Decomposition using Ladder Structure
% FBREC      Two-channel 2D Filterbank Reconstruction
% FBREC_L    Two-channel 2D Filterbank Reconstruction using Ladder Structure  
% 
%Multidimensional filtering (used in building block filter banks)
% SEFILTER2  2D seperable filtering with extension handling 
% EFILTER2   2D Filtering with edge handling (via extension)
% EXTEND2    2D extension
 
%Multidimensional sampling (used in building block filter banks)
% PDOWN      Parallelogram Downsampling
% PUP        Parallelogram Upsampling
% QDOWN      Quincunx Downsampling
% QUP        Quincunx Upsampling
% QUPZ       Quincunx Upsampling (with zero-pad and matrix extending)
% DUP        Diagonal Upsampling
% RESAMP     Resampling in 2D filterbank
% RESAMPZ    Resampling of matrix
% RESAMPC    Mex file used in RESAMP
%
%Polyphase decomposition (used in the ladder structure implementation)
% QPDEC      Quincunx Polyphase Decomposition
% QPREC      Quincunx Polyphase Reconstruction
% PPDEC      Parallelogram Polyphase Decomposition
% PPREC      Parallelogram Polyphase Reconstruction
%
%Support functions to avoid visual distortion (used in DFB)
% BACKSAMP   Backsampling the subband images of the directional filter bank 
% REBACKSAMP Re-backsampling the subband images of the DFB
%
%Support functions for generating filters 
% FFILTERS   Fan filters from diamond shape filters
% LD2QUIN    Quincunx filters from the ladder network structure
% MCTRANS    McClellan transformation
% MODULATE2  2D modulation
% 
%Other support functions
% COMPUTESCALE   Comupute display scale for PDFB coefficients
% SMTHBORDER Smooth the borders of a signal or image
% SNR        Compute the signal-to-noise ratio