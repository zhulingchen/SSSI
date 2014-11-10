function tracehead=AFD_makeSEGYtraceheaders(rx,ry,dt,sx,sy,ep)
% tracehead=AFD_makeSEGYtraceheaders(rx,ry,dt,sx,sy,ep)
% tracehead=AFD_makeSEGYtraceheaders(rx,ry,dt,sx,sy)
%
% AFD_makeSEGYtraceheaders will create a set of SEG-Y Revision 1
%  Trace headers wich can be used with the SEGY_Toolbox
%
% Inputs:
%   rx= must be a vector of x coordinates for the receivers.  There must be
%     one value per trace.
%   ry= is the y coordinate of the receivers and can either be a vector the 
%     same length as rx or a scalar.
%   dt= is the time sampling rate and can either be a vector the same 
%     length as rx or a scalar.
%   sx= is the x coordinate of the shot and can either be a vector the 
%     same length as rx or a scalar.
%   sy= is the y coordinate of the shot and can either be a vector the 
%     same length as rx or a scalar.
%   ep= is the source point number and can either be a vector the same 
%     length as rx or a scalar,
%
%
% Heather Lloyd 2010, Kevin Hall 2009, Chad Hogan 2004
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

if length(rx)<1
    me=MException('SEGY_maketraceheaders:InputNotCorrectSize',...
        'rx must be the same length as the traces these headers correspond to.');
    throw(me)
end
% create standard traceheader vector
tracehead=zeros(91,length(rx));
% set trace numbers to some standards
tracehead(1,:)=1:length(rx);
tracehead(2,:)=1:length(rx);
tracehead(3,:)=1;
tracehead(4,:)=1:length(rx);
tracehead(5,:)=1;
tracehead(6,:)=1:length(rx);
tracehead(7,:)=1:length(rx);
tracehead(8,:)=1;
tracehead(9,:)=1;
tracehead(10,:)=1;
tracehead(11,:)=1;
% set scale fields =1

tracehead([20,21,26,65,77,83,91],:)=1;
% check for decimals and scale appropriately
Dsx=0;for k=1:length(sx),ms=num2str(sx(k),6);m=strfind(ms,'.');dsx=length(ms)-m;if dsx>Dsx,Dsx=dsx;end;end;
Dsy=0;for k=1:length(sy),ms=num2str(sy(k),6);m=strfind(ms,'.');dsy=length(ms)-m;if dsy>Dsy,Dsy=dsy;end;end;
Drx=0;for k=1:length(rx),ms=num2str(rx(k),6);m=strfind(ms,'.');drx=length(ms)-m;if drx>Drx,Drx=drx;end;end;
Dry=0;for k=1:length(ry),ms=num2str(ry(k),6);m=strfind(ms,'.');dry=length(ms)-m;if dry>Dry,Dry=dry;end;end;
if ~isempty(Dsx) || ~isempty(Dsy) || ~isempty(Drx) || ~isempty(Dry)
    scl=10^(max([Dsx, Dsy, Drx, Dry]));
    sx=sx*scl;
    sy=sy*scl;
    rx=rx*scl;
    ry=ry*scl;
    tracehead(21,:)=scl;
end
% put the values in the header 
tracehead(22,:)=sx;
tracehead(23,:)=sy;
tracehead(24,:)=rx;
tracehead(25,:)=ry;

if nargin>5
    tracehead(5,:)=ep;
end

% convert dt to appropriate value and scale
if length(dt)>1
    dt=dt(2)-dt(1);
end
if dt<1
    dt=dt*10^6;
elseif dt<1000
    dt=dt*10^3;
end
dt=dt*ones(size(rx))/1000;
tracehead(40,:)=dt;



