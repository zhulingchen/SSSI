% CREWES Seismic toolbox
% Seismic processing tools 
%
%
% Transforms
%  FFTRL: forward Fourier transform for real vectors.
%  IFFTRL: inverse Fourier transform for real time series
%  IFKTRAN: inverse fk transform
%  FKTRAN_MC: forward fk transform using less memory than FKTRAN
%  FKTRAN: Forward fk transform
%  IFKTRAN_MC: memory conserving inverse fk transform
%  TPTRAN: tau-p transform for linear trajectories (slant stacking)
%  ITPTRAN: inverse tau-p transform for linear trajectories (filtered back projection)
%  TEST_TPTRAN: demo script for forward and inverse tau-p transforms
%
% Spectra
%  BURG: compute the Burg (maximum entropy) spectrum
%  MULTITAPER: estimate the spectrum of a short seismic trace using multitaper method
%
% Filters, convolution
%  CONVM: convolution followed by truncation for min phase filters
%  CONVZ: convolution then truncation for non-min phase filters
%  BURGPR: compute a Burg (maximum entropy) prediction error filter
%  FILTF: apply a bandass filter to a trace
%  FILTSPEC: designs the filter spectrum for FILTF
%  SECTCONV: runs convz on each trace in a matrix
%  SECTFILT: runs FILTF on each trace in a section
%  TVSW: time variant spectral whitening
%
% Amplitude adjustment
%  AEC: automatic envelope correction, a better AGC.
%  BALANS: match the rms power of one trace to another
%  CLIP: clips the amplitudes on a trace
%  GAINMUTE: Apply gain and top mute to a shot record
%  BANDWIDTH_XFER: transfer the bandwidth of one signal to another
%  
% Interpolation, resampling
%  RESAMP: resample a signal using sinc function interpolation
%  SINC: sinc function evaluation
%  SINCI: sinc function interpolation for time series without nan's
%  INTERPSINC: identical to SINCI except for the order of the input arguments.
%  SINCINAN: sinc function interpolation for signals with embedded nan's 
%  SINQUE: sinc function evaluation
%  SECTRESAMP: runs resamp on each trace in a seismic section
%
% Attributes
%  PICKER: make picks in a seismic matrix
%  IPICK: interactive interface to PICKER
%  FIND_ZERO_CROSSINGS: as the name says
%  INS_PHASE: Compute the instantaneous phase useing complex trace theory
%  INS_AMP: Compute the magnitude of the complex trace.
%  INS_FREQ: Compute the instantaneous frequency useing complex trace theory
%  FOMELFREQ: Compute a local frequency by Fomel's method
%  GABORFREQ: Compute a local frequency by a Gabor Method
%  DOM_FREQ: Compute the dominant frequency of a signal (centroid method)
%  TEST_INSTANTANEOUS_FREQ: compare various methods of local frequency
%
% Deconvolution, inversion
%  LEVREC: solve Tx=b using Levinson's recursion
%  DECONF: frequency domain stationary spiking deconvolution
%  DECONW: time domain stationary spiking decon (Wiener)
%  DECONB: time domain stationary spiking decon using Burg spectra
%
% Utilities
%  CONSTPHASE: estimate constant phase rotation between two signals by lsq
%  CONSTPHASE2: estimate cons. phs. rotation by systematic search over angles
%  TVCONSTPHASE: time variant version of constphase
%  TODB: convert to decibels
%  TOMIN: compute the minimum phase equivalent of a signal or wavelet
%  TOINV: compute the causal (min. phs.) inverse of a signal
%  TOINVF: compute freq. domain non-causal inverse of a signal
%  TOALL: compute the all pass equivalent of a signal
%  TOZERO: compute the zero phase equivalent of a signal
%  
% Auto and cross correlation
%  AUTO: single-sided autocorrelation
%  AUTO2: returns the two-sided autocorrelation 
%  MAXCORR: given two signals, find max crosscorrelation coef. and its lag
%
% Moveout and Traveltime adjustment
%  NMOR: normal moveout removal (forward and inverse)
%  NMOR_SRM: Normal movout for surface related multiples (forward and reverse)
%  STAT: static shift a trace
%
% Stacking
%  NMOR_CMP: Remove NMO from a shot gather and map traces to CMP locations
%  CMPSTACK: Common midpoint stack and gather creation.
%
% Seismic line geometry
%  NOMINAL_LINE_GEOM: create source/receiver positions for a 2D model
%  STACKING_CHART: create plots of source/receiver and midpoint/offset coords
%
