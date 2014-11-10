function [Dt,Dxz,fractx,nWide,iPlTerm,wncvar,iZlvl,pvel,svel] = ....
            readParmsCr4(parmFile)
%Read parameters for design of Wave Number Correction filters (.parc)
fid = fopen(parmFile,'rt');

parmStr(1,:) = 'Dt        ';
parmStr(2,:) = 'Dxz       ';
parmStr(3,:) = 'fractx    ';
parmStr(4,:) = 'nWide     ';
parmStr(5,:) = 'iPlTerm   ';
parmStr(6,:) = 'wncvar    ';
parmStr(7,:) = 'zLvl      ';
parmStr(8,:) = 'pvel      ';
parmStr(9,:) = 'svel      ';
nStr = size(parmStr);
parms = zeros(nStr(1),1);
indChar = 6;                %Indices of parameters in text format
nChParms = length(indChar);
iZlvl = zeros(5);
pvel = zeros(5);
svel = zeros(5);
while 1
    tline = fgetl(fid);
    if ~ischar(tline); break; end
    k = strfind(tline,'%');
    if isempty(k); k = length(tline)+1; end
    if k>5
        dline = tline(1:k-1);
        [parmName, dline2] = strtok(dline);
        ind = strmatch(parmName,parmStr);
        if ind > 0                          %Found parameter
            compCh = find(indChar-ind);
            nCh = length(compCh);
            if nCh==nChParms                %Retrieve numbers
                if ind < 6
                    value = strread(dline2, '%f');
                    parms(ind) = value;
                    %disp(value)
                else
                    if ind == 7
                        zLvl = strread(dline2, '%f');
                        disp(zLvl')
                    end
                    if ind == 8
                        pvel = strread(dline2, '%f');
                        disp(pvel')
                    end
                    if ind == 9
                        svel = strread(dline2, '%f');
                        disp(svel')
                    end
                end
            else                            %Retrieve strings
                wncvar = (strread(dline2, '%c'))';
                disp(wncvar)
            end
        else
            disp([parmName',' not found'])
        end
    end
end
fclose(fid);
%stop
%disp([parms(1),parms(2),parms(3),parms(4)])
Dt          = parms(1);
Dxz         = parms(2);
fractx      = parms(3);
nWide       = parms(4);
iPlTerm     = parms(5);
iZlvl = floor(zLvl/Dxz);
% pvel        = parms(3);
% svel        = parms(4);
% p2vel       = parms(5);
% s2vel       = parms(6);
%disp(parms(1:7))
% if mvXmax>lengthX; mvXmax = lengthX; end
% if mvZmax>lengthZ; mvZmax = lengthZ; end
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

