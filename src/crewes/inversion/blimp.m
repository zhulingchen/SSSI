function imp=blimp(trin,implog,t,flow,fhigh,delf)
% imp=blimp(trin,implog,t,flow,fhigh,delf)
%
% BLIMP estimates acoustic impedence from a seismic trace
% using a well log to provide the low frequency component.
% The algorithm is described in Ferguson and Margrave (1996 CREWES
% annual report).  This used to be called SEISINV2. Blimp is an 
% acronym for band limited impedance. 
%	
% trin ... input seismic trace
% implog ... input impedance log (in time)
% t ... time coordinate vector for trin
% flow ... lowest frequency in trin to keep
% fhigh ... highest signal frequency in trin
% delf ... width of Gaussian rolloff filter to be applied to
%	log at flow and trin at flow+delf
%	****** default min([5 flow/5]) *******
%
% G.F. Margrave, CREWES Project, U of Calgary, 1995-96
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
% check for row vectors and transpose if needed
aa=size(trin);
bb=size(t);
cc=size(implog);
if aa(1)==1
	trin=trin';
end
if bb(1)==1
	t=t';
end
if cc(1)==1
	implog=implog';
end
%integrate
if(nargin<6)
	delf=min([5 flow/5]); %gaussian on low end
end
impbl=rcs2impbl(trin,t,flow+delf,fhigh,delf);
%zero pad to impbl
impbl=padpow2(impbl);
%impbl=padpow2(impbl,1); %pad twice for better f sampling
%remove trend of log
p=polyfit(t,implog,1);
implog=implog-polyval(p,t);
implog=pad(implog,impbl);%zero pad
t2=xcoord(t(1),t(2)-t(1),implog);
%merge log and bandlimited impedance
imp=mergetrcs(implog,impbl,t2,flow,delf,fhigh);
imp=imp(1:length(trin));
%restore trand
imp=imp+polyval(p,t);
