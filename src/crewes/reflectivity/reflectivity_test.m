%************************************************************************
%
%*************  Main Program for Reflectivity Modeling       ****************
%
%************************************************************************
clear all

% Input parameter files

model = load('model.dat');
parameters = load('parameters.par');
geometry = load('geometry.dat');

%Calculate wave fields with the reflectivity codes

[uR,uZ,t] = reflectivity(model,parameters,geometry);

%Plot the final modeling results
figure;
plotseis(1,uR,t,geometry);ylabel('Time (second)');xlabel('Offset (km)');title('radial component')
figure;
plotseis(1,uZ,t,geometry);ylabel('Time (second)');xlabel('Offset (km)');title('vertical component')
