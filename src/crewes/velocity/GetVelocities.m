function data = GetVelocities(fig)

% data = GetVelocities(fig)
%
% Returns velocity-analysis data in a structure.
% fig .... figure number containing velocity analyis accomplished by the 
% VelocityAnalysis() function. The default is the current figure. If you
% have the wrong figure active, you'll get an empty result. 
%
% data.vrms : RMS Velocity in two-way time
% data.vave : Average velocity in two-way time
% data.vint : Interval velocity
% data.z    : Depth coordinate for interval velocity
% data.t    : Time coordinate for all velocities
%
% See also: VelocityAnalysis
%
% Chad Hogan, 2008
%
% $Id: GetVelocities.m,v 1.1 2009/05/25 21:00:56 cmhogan Exp $

if (nargin < 1)
    fig = gcf;
end

kids = get(fig, 'Children');
velax = kids(1);

ud = get(velax, 'UserData');

data.vrms = ud.vrmst;
data.vint = ud.vintt;
data.vave = ud.vavet;
data.z     = ud.zint;
data.t     = ud.t;