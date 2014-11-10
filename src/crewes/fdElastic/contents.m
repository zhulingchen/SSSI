% CREWES elastic finite-difference seismic modelling
% The top boundary is a free surface suitable for propagating surface
% (Rayleigh) waves.
%
% Script drivers
% fdModDown ..... build a geological model from top to bottom
% fdModUp ....... build a geological model from bottom to top
% fdLayers ...... run a forward model, select parameters by modifying code
%                    The control menu code is in this module
%
% Modelling subroutines
% fdModUD ....... Make Vp, Vs and density vectors, X/Z boundary matrices
% fdLayDoc ...... Carry f-d model through several time steps
% plotPack ...... plot a snapshot (state) with a choice of plot types
% plotcolourk ... Plot vector data with direction coded around the spectrum
% plotcoloura ... Plot vector data with direction coded around the spectrum,
%                   with absolute amplitude (consistent scaling)
% plotptwist2 ... Plot pressure and twist functions of the displacement field
%
% Private subroutines
% fdModUD
% fdInitGen
% fdInitMod2
% fdLoadModel
% fdInitPad
% fdInitTrace
% fdInitArray1
% fdInitArray
% fdInitBound
% fdInitForce
% fdLoadPrior
% fdSEGY2
% plotTraces
% padright6
% zbounds
%
% Dot mat files
% The two files wncEdge.mat and wncUnity.mat (wnc: wave number correction)
% are generated to work with surface waves. The generating software, which
% accounts for body waves, is not yet ready to be released.
%
% The main script, fdLayers, writes ASCII files and a parameter file for
% command ascii2segy which is also part of the 2008/9 CREWES software release.
