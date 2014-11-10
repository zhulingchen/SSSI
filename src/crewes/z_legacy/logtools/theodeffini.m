function [sonicnum,densitynum,waveletnum,densityopt,holealg,a,m,theotype,name]=...
				theodeffini
% [sonicnum,densitynum,waveletnum,densityopt,holealg,a,m,theotype,name]=theodeffini
%
% Call this to complete the theogram definition dialog. Return values mean:
%	sonicnum = number of the sonic section selected for theogram
%	densitynum = number of the density section selected for theogram
%	waveletnum = number of the wavelet selected for theogram
%	densityopt = density option selection 
%                    1: use Gardners exclusively, 
%                    2: use density logs but fill holes with Gardners 
%                    3: use density logs exclusively
%	holealg = hole filling algorithm selection 
%                 1: constant 2: linear 3: mean
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
	params=get(gca,'userdata');
	if( params==-1 ) %test for a cancel
		sonicnum=-1;
		return;
	end
	sonicnum=params(1);
	densitynum=params(2);
	waveletnum=params(3);
	densityopt=params(4);
	holealg=params(5);
%		4: layer mean 5: layer trend
%   domain = 1: compute rcs in time  or 2: compute rcs in depth
%	a= scalar coefficient in gardners relation. Will be -1 if densopt=3
%	m= scalar exponent in gardners relation. Will be -1 if densopt=3
%	name = name chosen for the rc section
	a=params(6);
	m=params(7);
	theotype=params(8);
	name=setstr(params(8:length(params)));
