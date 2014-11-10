function amat=event_dip(amat,t,x,tlims,xlims,amp)
% EVENT_DIP: inserts a dipping (linear) event in a matrix
%
% amat=event_dip(amat,t,x,tlims,xlims,amp)
%
% EVENT_DIP inserts a dipping (linear) event in a matrix.
%
% amat ... the matrix of size nrows by ncols
% t ... vector of length nrows giving the matrix t coordinates
% x ... vector of length ncols giving the matrix x coordinates
% tlims ... vector of length 2 giving the t limits of the event
% xlims ... vector of length 2 giving the x limits of the event
% amp ... vector of length 2 giving the amplitudes at either end
%	of the event. Intermediate amplitudes are interpolated in x.
%   *************** default [1 1] *****************
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

if(nargin<6)
    amp=[1 1];
end

%loop over columns
nc= between(xlims(1),xlims(2),x,2);

spike=zeros(length(t),1);
spike(1)=1;

if(length(amp)==1) amp=[amp amp]; end

if(nc~=0)
	tmin=t(1);
	tmax=t(length(t));
	dt=t(2)-t(1);
	for k=nc
		tk = interp1(xlims,tlims,x(k));
		a = interp1(xlims,amp,x(k));
		if( between(tmin,tmax,tk) )
			amat(:,k)=amat(:,k) + a*stat(spike,t,tk-tmin);
		end
	end
end
