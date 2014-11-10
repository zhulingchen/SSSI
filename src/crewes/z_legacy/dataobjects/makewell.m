function well=makewell(name,location,wellog)
%
% well=makewell(name,location,wellog)
%
% Make a earthobject representing a well with name, location, and 
% (optionally) a log
%
% name = string giving the well name
% location = 3 or 4 element vector giving: [x y elevation inline_distance]
% wellog = random Earth Object containing the log. If the log is in depth 
%          then it should have two data fields: 'depth' & 'samples' and 
%          be of datatype 'zlog', else it should have data fields 
%          'time' and 'samples' and be of type 'tlog'
%
% by G.F. Margrave
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

% Make an empty well object
well=contobj(name,'well');

% put the location information in
well=objset(well,'location',location);

% put the log in
if( nargin > 2)
	nam=objget(wellog,'name');
	well=objset(well,nam(:)',wellog);
end
