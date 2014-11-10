function vrms=vint2vrms(vint,t,tout)
% VINT2VRMS: convert interval to rms velocity
% 
% vrms=vint2vrms(vint,t,tout)
% vrms=vint2vrms(vint,t)
%
% VINS2VRMS computes rms velocity as a function of time
% given interval velocity as a function of time.
%
% vint = input interval velocity vector
% t = input time vector to go with vint
% tout = vector of output times at which vrms estimates are
%	desired. Requirement tout >= t(1)
%*********** default tout=t *********
%
% NOTE: Times can be either one-way or two-way.
%
%
% G.F. Margrave June 1995, CREWES Project
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

%

%test input arguments
if(length(vint)~=length(t))
	error('vint and t must have same lengths')
end

if nargin<3
	tout=t;
end

%force to column vectors
vint=vint(:);t=t(:);tout=tout(:);

if(nargin < 3)
	%integrate vint

	vrms=zeros(size(vint));
	nt=length(t);
    dt=diff(t);
	i1=1:nt-1;
    
	vrms(i1)=cumsum(dt.*(vint(i1).^2));
	vrms(i1)=sqrt(vrms(i1)./(t(i1+1)-t(1)));
	vrms(end)=vrms(end-1);
% 	dt=[t(1); diff(t)];
%     tmp=cumsum(dt.*(vint.^2));
%     vrms=sqrt(tmp./t);

else
   ind=find(tout<t(1));
   if(~isempty(ind))
		error('tout must be greater than t(1)');
   end
	
	nt=length(t);
	i1=1:nt-1;
	vrms2=zeros(size(t));
	dt=diff(t);
	vrms2(i1)=cumsum(dt.*(vint(i1).^2));

	vrms2(i1)=vrms2(i1)./(t(i1+1)-t(1));
	vrms2(end)=vrms2(end-1); %densly sample rms^2 from surface

   vrms=pwlint(t,sqrt(vrms2),tout);
   %check for nans at the end
   ind=find(isnan(vrms));
   if(~isempty(ind))
       vrms(ind)=vrms(ind(1)-1);
   end
	
%   vrms22=pcint(t,vrms2,tout); %interpolate surface rms^2 values at tout

%   vrms=zeros(size(tout));
%   i1=1:length(tout)-1;
%   vrms(i1)=sqrt((vrms22(i1+1).*tout(i1+1)-vrms22(i1).*tout(i1))...
%      ./(tout(i1+1)-tout(i1)));
%   vrms(end)=vrms(end-1);
   
end
