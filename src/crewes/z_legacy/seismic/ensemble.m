function [x,tx,a,ta,w,tw]=ensemble(ntraces,nsamps,dt,fdom)
% [x,tx,a,ta,w,tw]=ensemble(ntraces,nsamps)
%
% ENSEMBLE generates an ensemble (profile) of synthetic traces
% for use with ensemble averaging experiments
%
% ntraces= number of rows in the matricies x and a
% nsamps= number of columns in x, keep this a power of 2!! 
% dt= time sample rate desired on output
% fdom= dominant frequency of the waveform 
% x= an ntraces by nsamps matrix of gaussian random noise where
%    each row has been convolved with the waveform w.
% tx= time coordinate vector for x
% a= the two-sided autocorrelation matrix of x, ntraces by 2*nsamps
%    (There is an extra sample added to keep things a power of two)
% ta= time coordinate vector for a
% w= minimum phase waveform with a dominant frequency of fdom
% tw= time coordinate vector for w
%
% by G.F. Margrave, July, 1991
%
% NOTE: It is illegal for you to use this software for a purpose other
% than non-profit education or research UNLESS you are employed by a CREWES
% Project sponsor. By using this software, you are agreeing to the terms
% detailed in this software's Matlab source file.
 
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by 
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
% 
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code
%
% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the 
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may 
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers 
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any 
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE
 if nargin<4, fdom=15;end
 if nargin<3, dt=.002;end
 if nargin<2, nsamps=256;end
 if nargin<1, ntraces=20;end
% generate the random matrix
 mtx=randn(ntraces,nsamps);
 tx=xcoord(0.,dt,mtx(1,:));
% generate the wavelet
 nw=nsamps;
 [w,tw]=wavemin(dt,fdom,(nw-1)*dt);
% convolve them to get x
 x=convm(mtx,w);
% autocorrelate x
 a=auto2(x);
 ta=xcoord(-max(tx),dt,a(1,:));
