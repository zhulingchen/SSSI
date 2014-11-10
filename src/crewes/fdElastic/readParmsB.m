function [Dt,Dxz,xMin,lengthX,lengthZ,nstep,nframes,shotDepth,shotX,....
        centreFreq,energySrc,cycleSrc,timeSrc,trX,trZ,iLbnd,iRbnd,iTbnd,iBbnd,....
        mvXmax,mvZmax,mvPlot,mvAmp,mvClip,gfdFile,wncvar,mvTif,contBlk] = ....
            readParmsB(parmFile)
%function [parms] = readParms(parmFile)
%
%parmFile ... name of file containing finite-difference parameters, ending in .parm
%
%Dt      .... FD sample rate in seconds
%Dxz     .... FD sample rate in (metres)
%xMin    .... Left edge of the FD model extracted from the geological model
%lengthX .... X-length of the FD model
%lengthZ .... Z-height of the FD model
%nstep   .... Number of FD time steps
%nframes .... Number of movie frames to display or plot
%shotDepth  . Z (depth) of the FD initializing source
%shotX   .... X (from the FD model left) of the initializing source
%centreFreq . The centre frequency of the FD source wavelet in time
%energySrc  . The directivity code of the energy source
%cycleSrc.... No. of cycles (2-windowed cosine only)
%timeSrc .... Code to specify wavelet: 1-Ricker, 2-windowed cosine
%trX     .... X-position of the 'well' where trace (time) data will be collected
%trZ     .... Z-depth of the 'line' where trace (time) data will be collected
%iLbnd   .... Boundary code left
%iRbnd   .... Boundary code right
%iTbnd   .... Boundary code top
%iBbnd   .... Boundary code bottom
%mvXmax  .... X-length of movie frames to display or plot
%mvZmax  .... Z-depth of movie frames to display or plot
%mvPlot  .... The code number of the movie plot
%mvAmp   .... The amplitude for mvPlot=4 (larger is higher amplitude)
%mvClip  .... The movie amplitude clip level (1 is unclipped)
%gfdFile .... The geological definition file (within quotes)
%wncvar  .... The wavenumber correction file (within quotes), ('' is none) 
%mvTif   .... The root name of the movie tiff file (within quotes), ('' is none)
%contBlk .... The indicator of a 'cont' (continuous) or 'block' (blocked)
%               geological model
%
%To pick up an additional parameter, see **** comments,
%   and add it to the return list
%
% P.M. Manning, Dec 2011
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
fid = fopen(parmFile,'rt');


parmStr(1,:) = 'Dt        ';
parmStr(2,:) = 'Dxz       ';
parmStr(3,:) = 'xMin      ';
parmStr(4,:) = 'lengthX   ';
parmStr(5,:) = 'lengthZ   ';
parmStr(6,:) = 'nstep     ';
parmStr(7,:) = 'nframes   ';
parmStr(8,:) = 'shotDepth ';
parmStr(9,:) = 'shotX     ';
parmStr(10,:) = 'centreFreq';
parmStr(11,:) = 'energySrc ';
parmStr(12,:) = 'cycleSrc  ';
parmStr(13,:) = 'timeSrc   ';
parmStr(14,:) = 'trX       ';
parmStr(15,:) = 'trZ       ';
parmStr(16,:) = 'iLbnd     ';
parmStr(17,:) = 'iRbnd     ';
parmStr(18,:) = 'iTbnd     ';
parmStr(19,:) = 'iBbnd     ';
parmStr(20,:) = 'mvXmax    ';
parmStr(21,:) = 'mvZmax    ';
parmStr(22,:) = 'mvPlot    ';
parmStr(23,:) = 'mvAmp     ';
parmStr(24,:) = 'mvClip    ';
parmStr(25,:) = 'gfdFile   ';
parmStr(26,:) = 'wncvar    ';
parmStr(27,:) = 'mvTif     ';
parmStr(28,:) = 'contBlk   ';
%Add another file variable name to this parmStr array *****
nStr = size(parmStr);
parms = zeros(nStr(1),1);

indChar = [25 26 27 28];    %Increment the indices if new parm is earlier *****
%If the added variable is a character string, include its number in indChar *****
nChParms = length(indChar);
while 1
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    k = strfind(tline,'%');
    if isempty(k); k = length(tline)+1; end
    if k>5
        dline = tline(1:k-1);
        [parmName] = strread(dline, '%s',1);
        %[parmName Cstring] = sscanf(dline, '%s %s');
        %disp(parmName)
        ind = strmatch(parmName,parmStr);
        if ind > 0                          %Found parameter
            %parmsC(ind,:) = Cstring;
            compCh = find(indChar-ind);
            nCh = length(compCh);
            if nCh==nChParms                %Retrieve numbers
                [parmName value] = strread(dline, '%s %f');
                %[parmName value] = sscanf(dline, '%s %f');
                %[value] = strread(dline, '%f');
                parms(ind) = value;
            else                            %Retrieve strings
                %[valueC] = strread(dline, '%s %s,');
                %[valueC] = strread(dline, '%s,');
                %[parmName Cstring] = strread(dline, '%s');
                [Cstring] = strread(dline, '%c');
                %parms(ind) = Cstring;
%                 cmd = ([parmName '= disp('''  Cstring ')''']);
%                 disp(cmd)
%                 eval(cmd)
                %disp(Cstring)
                eval([Cstring;';'])
            end
        else
            disp([parmName',' not found'])
        end
    end
end
fclose(fid);
%stop
%disp([parms(1),parms(2),parms(3),parms(4)])
Dt         = parms(1);
Dxz        = parms(2);
xMin       = parms(3);
lengthX    = parms(4);
lengthZ    = parms(5);
nstep      = parms(6);
nframes    = parms(7);
shotDepth  = parms(8);
shotX      = parms(9);
centreFreq = parms(10);
energySrc  = parms(11);
cycleSrc   = parms(12);
timeSrc    = parms(13);
trX        = parms(14);
trZ        = parms(15);
iLbnd      = parms(16);
iRbnd      = parms(17);
iTbnd      = parms(18);
iBbnd      = parms(19);
mvXmax     = parms(20);
mvZmax     = parms(21);
mvPlot     = parms(22);
mvAmp      = parms(23);
mvClip     = parms(24);
%Add new numeric input names to the copy from parms *****
if cycleSrc<1; cycleSrc = 1; end
if mvXmax>lengthX; mvXmax = lengthX; end
if mvZmax>lengthZ; mvZmax = lengthZ; end
%Then add name to the parameter list                *****
%   and add name to the calling programs            *****

% itrAcq     = parms(21);
% mvTif      = parms(22);

%         for iParm = 1:nStr
%             if strcmpi(parmName,parmStr(iParm,:))
%                 %parms(iParm) = sscanf(rest,'%f');
%                 parms(iParm) = value;
%                 disp('hit')
%                 break
%             end
%         end

