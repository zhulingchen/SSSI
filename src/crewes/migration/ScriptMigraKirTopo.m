%  Nov-2011

% Script Migra Kirchhoff from topography - Converted wave
% Method shot profile
% An example

% Requirements:
%   * The time tables for Sources and Receivers
%   * The seismic data of the horizontal component. Better if data Zero 
%       Phase to obtain the correct depths.
%   * Aditional data required: field record, its coordinates (t,x),
%       velocity model of the source, velocity model of the 
%       receiver, coordinates of these models (xv, xv), coordinates
%       of the energy source (xshot, zshot), surface topography (zt)
%       
%
% It can also be used for pure PP reflections, using the seismic data of
% the vertical component as an input, and the P-wave time tables for the
% receivers.


%% Migration Kirchhoff - Three shots


[smigt41x,smigz41x,tmig41x,xmig41x,vmod41x]=kirk_shotcz(topo41xs,t,x,zt,xshot1,vp2,vs2,xv,zv,tr4s19(:,:,1),tr4rall,[]);
[smigt45x,smigz45x,tmig45x,xmig45x,vmod45x]=kirk_shotcz(topo45xs,t,x,zt,xshot5,vp2,vs2,xv,zv,tr4s19(:,:,5),tr4rall,[]);
[smigt49x,smigz49x,tmig49x,xmig49x,vmod49x]=kirk_shotcz(topo49xs,t,x,zt,xshot9,vp2,vs2,xv,zv,tr4s19(:,:,9),tr4rall,[]);

save migraKirkcTopoResult

%% Stack Migra
stkMig19=smigz41x+smigz45x+smigz49x;
plotimage(stkMig19,zv,xmig45x)
title('Stack C-wave Mig Kir 1,5,9 ')
grid


