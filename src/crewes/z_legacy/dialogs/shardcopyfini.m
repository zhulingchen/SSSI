function [ps_scale,xlength,ylength,xscale,yscale]=shardcopyfini
% [ps_scale,xlength,ylength,xscale,yscale]=shardcopyfini
%
% shardcopyfini terminates the hardcopy dialog and returns the
% user supplied scales for x and y as well as the final plot 
% size after scaling (xlength and ylength) and the scale factor
% needed for CHVSUB (ps_scale). The calling program will probably
% want to report these to the user somehow as it is up to her to
% actually submit the plotfile to CHVSUB and type the ps_scale
% value into CHVSUB
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
a=get(gca,'userdata');
set(gca,'userdata',[]);
if(a(1)==-999.)
	ps_scale=-999.;
	return
end
ps_scale=a(1);
xlength=a(2);
ylength=a(3);
xscale=a(4);
yscale=a(5);
