function [tvs,trow,fcol]=tvspec(trin,t,win_length,ovlpct)
% [tvs,trow,fcol]=tvspec(trin,t,win_length,ovlpct)
% [tvs,trow,fcol]=tvspec(trin,t,win_length)
% 
% TVSPEC computes a complex valued time variant Fourier spectrum. Most
% SPIN applications will want only the amplitude spectrum obtained by
%  tvs=abs(tvs); 
% The raised cosine style MWINDOW is used with a hard-wired 30% taper.
%
% trin= input trace
% t= time coordinate vector
% win_length= length of time window (sec)
% ovlpct= overlap percentage of the windows
%    ********* default= 80.0 **********
% tvs= ouput time variant spectal matrix (complex). There will be
%      one spectrum per row of tvs.
% trow= vector specifying the time of each row
% fcol= vector specifying the frequency of each column
%
% by G.F. Margrave, May 1991
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
% set defaults
 if nargin<=3
  ovlpct=80.;
 end
% compute the window function 
  dt=t(2)-t(1);
  ntwind=2.^nextpow2(win_length/dt+1); % pad window to next power of 2
  window=mwindow(ntwind,30);
% main loop setup
  nend=ntwind;
  nt=length(trin);
  factor=(100.-ovlpct)/100.; % window moveup factor
  numspec= floor((nt-ntwind)/(factor*ntwind))+1;% number of spectra to becomputed
  nspec=1;
% do first spectrum and pre-allocate arrays
% note: all spectra are computed at length nt even though they have
%     only ntwind live samples. This achieves interpolation to the 
%     frequency sample rate needed for application via the SPIN 
%     algorithm.
    temp= [trin(1:nend).*window zeros(1,nt-nend)];
    [spc,fcol]=fftrl(temp,t);
  tvs= zeros(numspec,length(spc));
  tvs(nspec,:)=spc;
  nend=floor(nend+factor*ntwind);
  nspec=nspec+1;
% now enter main loop
  while nend<=nt
    nbeg=nend-ntwind+1; % first sample in window
    trow(nspec)=dt*((nend+nbeg)/2.-1);
% apply window and compute spectrum 
    temp= [zeros(1,nbeg-1) trin(nbeg:nend).*window zeros(1,nt-nend)];
    [spc,fcol]=fftrl(temp,t);
    tvs(nspec,:)=abs(spc);
    nend=floor(nend+ntwind*factor); % last sample in next window
    nspec=nspec+1;
  end
      
