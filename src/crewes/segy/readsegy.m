function [traces,sampint,texthead,binaryhead,extendedhead] = readsegy(sgyfile)

% function [traces,sampint,texthead,binaryhead,extendedhead] = readsegy(sgyfile)
%
% Reads segy data into a structured array format in Matlab workspace.
%
% traces    = a 2-D array containing the data of interest
% sampint   = sample interval in seconds (optional)
%

[traces,texthead,binaryhead,extendedhead]=SEGY_read(sgyfile);
sampint = SEGY_getHeader(binaryhead,'hdt')/1000;







