function [pm,p]=synseis_press(rin)
% [pm,p]=synseis_press(r)
%
% SYNSEIS_PRESS computes a complete 1-D seismogram assuming a pressure 
% recording. Based of the algorithm in 
%	Waters, "Reflection Seismology", 1981, John Wiley, pp 128-135
% and as discussed in
%   Margrave, "Numerical Methods of Reflection Seismology", 2001
%
% r ... input reflection coeficients, regularly sampled in TIME 
%			(2-way time)
% pm ... output 1-D impulse response, primaries plus multiples,
%	sampled at the same times as r.
% p ... primaries only showing the effects of transmission losses
% 
% G.F. Margrave, CCR, Oct 1994
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
%initialize a few things
r=rin(2:length(rin));
r0=rin(1);
d=zeros(size(rin));
d(1)=1;
pm=zeros(size(rin));
p=zeros(size(rin));
pm(1)=1; %direct wave
p(1)=1;
%loop over output times
for k=1:length(r)
	
	%zero upgoing wave at each step
	u=0.0;%upgoing primaries plus multiples
	up=0.0;%upgoing, primaries only wave
	%step from r(k) to the surface
	for j=k:-1:1
		%update downgoing wave
		d(j+1)=(1-r(j))*d(j) -r(j)*u;
		%update upgoing wave
		u=r(j)*d(j) + (1+r(j))*u;
		if( j==k )
			up = u;
		else
			up = (1+r(j))*up;
		end
	end
	%step to surface
    d(1)= -r0*u;
	pm(k+1)=u+d(1);
	p(k+1)=up;
end
