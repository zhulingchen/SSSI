function [mfilt,tm]=matchf(trin,trdsign,t,fsmoth,flag)
% [mfilt,tm]=matchf(trin,trdsign,t,fsmoth,flag)
%
% MATCHF designs a frequency domain match filter.
% The match filter will always have time zero in the middle and
% should be applied with convz.
%
% trin= input trace to be matched to trdsign
% trdsign= input trace which is to be matched
% t= time coordinate vector for trin
% ***** note: trin and trdsign must be the same length
% fsmoth= length of frequency smoother in Hz.
% flag=0 ... smooth the filter spectrum prior to inverse transform
%     =1 ... smooth the spectra of trin and trdsign before division
%     =2 ... both 0 and 2   
%  ********* default= 0 ********      
%
% mfilt= output mlength match filter
% tm= time coordinate for the match filter
%
% by G.F. Margrave, June 1991
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
 if nargin<5, flag=0; end
% forward fft's
 [Trin,f]=fftrl(trin,t);
 [Trdsign,f]=fftrl(trdsign,t);
% smooth
 if flag ~=0,
   Trin=convz(Trin,boxf(fsmoth,f));
   Trdsign=convz(Trdsign,boxf(fsmoth,f));
 end
% compute operator spectrum
 Mf=Trdsign./Trin;
% smooth
 if flag~=1,
   Mf=convz(Mf,boxf(fsmoth,f));
 end
% inverse transform
 mfilt=ifftrl(Mf,f);
 mfilt=fftshift(mfilt);
 tm=xcoord(-max(t)/2.,.002,mfilt);
   
 
 
