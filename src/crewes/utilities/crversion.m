function [version,versiontext] = crversion()
%
% function [version,versiontext] = crversion()
%
% CREWES Matlab toolbox zip or tgz files from the website contain
% a file called svninfo.txt. crversion returns the svn version number
% from svninfo.txt if it exists, or [] if it does not
%

version     = []; %assume failure
versiontext = ''; %assume failure
SVNINFO = which('crversion.txt');

if isempty(SVNINFO)
   %disp('CREWES toolbox version number was not found');
else
    try
        fid   = fopen(SVNINFO);
        tline = fgetl(fid);
        versiontext = tline;
        version = str2double(tline(strfind(tline,':')+2:end));
        if isnan(version) version =0.0; end
            
        fclose(fid);       
    catch ex
        disp(ex.message);
    end
end