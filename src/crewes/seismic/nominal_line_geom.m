function [xshots,xrecs]=nominal_line_geom(x,nshots,spreadlength,offedge)
% NOMINAL_LINE_GEOM: compute source and receiver positions for 2D synthetic
%
% [xshots,xrecs]=nominal_line_geom(x,nshots,offedge)
%
% Function distributes shots evenly along a 2D line while keeping them a
% reasonable distance from the edge. Receiver interval is taken as twice the
% model grid spacing so that the cmp spacing is equal to the model grid.
% Receivers span the model and the receiver spread is identical for all
% shots. The shot positions are chosen to be between adjacent receivers.
%
% x ... x coordinate vector for the numerical model. That is x(1) is the
%       beginning of the model, x(end) is the end, x(2)-x(1) is the numerical
%       grid spacing, and x(k+1)>x(k) for any k in [1 length(x)].
% nshots ... number of shots desired.
% spreadlength ... length of geophone spread (equivalent to maximum offset)
% ************* default = x(end)-x(1) **************
% offedge ... offset of first and last shot from the model edge
% ************** default: .05*(x(end)-x(1)) **********
%
% xshots ... vector of shot locations
% xrecs ... cell array of receiver locations (identical for all shots). The
%           length(xrec) will equal the length(xshots)
%
% G.F. Margrave, CREWES Project, August 2013
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
xmax=max(x);
dx=x(2)-x(1);
if(dx<0)
    error('x coordinates must be increasing')
end

if(nargin<4)
    offedge=.05*(x(end)-x(1));
end
if(nargin<3)
    spreadlength=x(end)-x(1);
end

%make spread length an even number of geophone spacings
dg=2*dx;
spreadlength=round(spreadlength/(2*dg))*2*dg;
nspr=spreadlength/(dg);
xrnom=(-nspr:nspr)*dx;%nominal geometry

%shot locations first
if(nshots==1)
    %special case
    xshots=xmax/2;
else
    dxshot=(xmax-2*offedge)/(nshots-1);
    xshots=offedge:dxshot:xmax-offedge;
    %round to the grid
    xshots=round(xshots/dx)*dx;
    %shift shot locations to the nearest odd numbered grid point
    xshots=ceil(xshots/(2*dx))*2*dx;
    %note that the grid points corresponsding to the shots are
    %ishots=round(xshots/dx)+1. The +1 means that when x shot is an even
    %number the corresponding grid point is an odd number.
end
%receiver locations
xrecs=cell(1,nshots);
for k=1:nshots
    %The basic idea is that we want receivers at every other grid point (so
    %that the cmp sampling is equal to the grid spacing). We will also keep
    %the receivers the same for each shot by placing them at even numbered
    %grid points. Note that the shots were placed at odd numbered grid
    %points.
    xr=xshots(k)+xrnom;
    if(xr(1)<x(1))
        del=ceil((x(1)-xr(1))/dg)*dg;
        xr=xr+del;
    end
    if(xr(end)>x(end));
        del=ceil((xr(end)-x(end))/dg)*dg;
        xr=xr-del;
    end
    xrecs{k}=xr;
end
