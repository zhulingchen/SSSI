function imp=seisinv1(trin,t,znot,flow,fhigh,n)
% imp=seisinv1(trin,t,znot,flow,fhigh,n)
%
% SEISINV1 estimates acoustic impedence from a seismic trace
% using a Burg prediction filter to provide the low frequency
% component. 
%
% trin ... input seismic trace
% t ... time coordinate vector for trin
% znot ... first impedance value
% flow ... lowest frequency in trin to keep
% fhigh ... highest signal frequency in trin
% n ... number of points in Burg prediction filter
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
% check for row vector and transpose if needed
aa=size(trin);
bb=size(t);
if aa(1)==1
	trin=trin';
end
if bb(1)==1
	t=t';
end
%pad if needed
trpad=padpow2(trin);
tpad=xcoord(t(1),t(2)-t(1),trpad);
%compute spectrum
[Trin,f]=fftrl(trpad,tpad);
%predict the lows
Trout=predlowf(Trin,f,flow,fhigh,n);
%ifft and remove pad
rcs=ifftrl(Trout,f);
if(length(rcs)~=length(trin))
	rcs=rcs(1:length(trin));
end
%integrate
imp=rcs2imp(rcs,znot);
