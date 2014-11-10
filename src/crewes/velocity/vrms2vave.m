function vave=vrms2vave(vrms,t)
% VRMS2VAVE: convert rms velocity to average
%
% vave=vrms2vave(vrms,t)
%
% VRMS2VAVE computes average velocity as a function of time
% given rms velocity as a function of time. The method is simply to call
% vrms2vint followed by vint2vave. Non-physical vint values resulting from
% the vrms function are thrown out and new values are interpolated from
% neighbors.
%
% vrms = input rms velocity vector
% t = input time vector to go with vrms
% tout = vector of output times at which vave estimates are
%	desired. Requirement tout >= t(1)
%*********** default tout=t *********
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

%test input arguments
if(length(vrms)~=length(t))
	error('vrms and t must have same lengths')
end

%force column vectors
vrms=vrms(:);
t=t(:);

vint=vrms2vint(vrms,t,1);
vave=vint2vave(vint,t);