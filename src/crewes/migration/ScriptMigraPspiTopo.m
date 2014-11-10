
%  Nov-2011

% Script Migra Pspi from topography - Converted wave
% Method shot profile
% An example

% Requirements:
%   * The velocity models for P and S waves
%   * The seismic data of the horizontal component. Better if data Zero 
%       Phase to obtain the correct depths.
%   * Aditional data required: field record, its coordinates (t,x),
%       velocity model of the source (velp), velocity model of the 
%       receiver (vels), coordinates of these models (xv, xv), coordinates
%       of the energy source (xshot, zshot), surface topography (zt)
%       
% It can also be used for pure PP reflections, using the seismic data of
% the vertical component as an input, and the P-wave time tables for the
% receivers.
%


%% Migration PSPI - Three shots


[smdc41x,smcc41x,illum41x]=pspi_shot_cwavez(topo41xs,t,x,velp,vels,xv,zv,xshot1,zshot1,zt,[0 60],0.001);
[smdc45x,smcc45x,illum45x]=pspi_shot_cwavez(topo45xs,t,x,velp,vels,xv,zv,xshot5,zshot5,zt,[0 60],0.001);
[smdc49x,smcc49x,illum49x]=pspi_shot_cwavez(topo49xs,t,x,velp,vels,xv,zv,xshot9,zshot9,zt,[0 60],0.001);

save migraPspicTopoResult
%% Migration Stack
% stkmigcc=smcc41x+smcc45x+smcc49x;
% plotimage(stkmigcc,zv,xv)


%%

