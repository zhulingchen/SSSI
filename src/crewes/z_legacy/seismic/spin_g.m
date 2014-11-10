function trout=spin_g(trin,t,smooth_f,smooth_t,phase,delt,win_length)
% trout=spin_g(trin,t,smooth_f,smooth_t,phase,delt,win_length)
% trout=spin_g(trin,t,smooth_f,smooth_t,phase,delt)
% trout=spin_g(trin,t,smooth_f,smooth_t,phase)
% trout=spin_g(trin,t,smooth_f,smooth_t)
% trout=spin_g(trin,t,smooth_f)
% trout=spin_g(trin,t)
% 
% Time variant spectral inversion 'generic' method.
% Algorithm is described in: 
%               "A New Method of Time Variant Spectral Inversion"
%                 by G.F. Margrave, 1989 ICGC
%
% trin= input trace
% t= time coordinate for trin
% smooth_f= frequency smoother length in Hz
%   *********** default = 10 hz *************
% smooth_t= temporal smoother length in seconds
%   *********** default = .4 secs **********
% phase= phase flag, 0= zero phase, 1= minimum phase
%   *********** default = 0 ********** 
% delt= temporal increment (seconds) for computing the tvs of trin
% *************** default = .1 ************
% win_length= window length (seconds) for computing the tvs of trin
% ************** default = .3 *************
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
 if nargin < 7, win_length=.3; end
 if nargin < 6, delt=.1; end
 if nargin < 5, phase=0; end
 if nargin < 4, smooth_t=.4; end
 if nargin < 3, smooth_f=10; end
 
% compute the time variant amplitude spectrum of the input trace
 [tvs,trow,fcol]=tvspec(trin,t,delt,win_length);
 tvs=abs(tvs)+i*eps; % this peculiar statement computes the amplitude spectrum
%                    while keeping the tvs complex.
 nt=length(trow);nf=length(fcol);
% smooth the tvs in frequency and time using a 2-d convolve
% The time coordinate is constant for each row of the tvs and
% likewise, frequency is constant for each column.
% 
 df=fcol(2)-fcol(1);
 dt=trow(2)-trow(1);
 nfs=ceil(smooth_f/df);
 nts=ceil(smooth_t/dt);
 smooth_2d= ones(nts,nfs);
 tvs=conv2(tvs,smooth_2d);
% now extract the "central part" of the smoothed tvs
 nf2=round(nfs/2);nt2=round(nts/2);
 tvs=tvs(nt2:nt+nt2-1,nf2:nf+nf2-1); 
% now compute the minimum phase spectrum of each row
 if phase==1,
   tvs = -log(real(tvs))+i*eps;
 end
trout= tvs;
