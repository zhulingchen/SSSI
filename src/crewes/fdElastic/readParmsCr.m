function [Dt,Dxz,pvel,svel,fractx,nWide,iPlTerm,wncvar] = ....
            readParmsCr(parmFile)
%Read parameters for design of a Wave Number Correction filter (.parc)
%function [parms] = readParms(parmFile)
%Key code 3 Read parameters from a text file *****************************
fid = fopen(parmFile,'rt');


parmStr(1,:) = 'Dt        ';
parmStr(2,:) = 'Dxz       ';
parmStr(3,:) = 'pvel      ';
parmStr(4,:) = 'svel      ';
parmStr(5,:) = 'fractx    ';
parmStr(6,:) = 'nWide     ';
parmStr(7,:) = 'iPlTerm   ';
parmStr(8,:) = 'wncvar    ';
nStr = size(parmStr);
parms = zeros(nStr(1),1);

indChar = (8);
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
        if ind > 0
            %parmsC(ind,:) = Cstring;
            compCh = find(indChar-ind);
            nCh = length(compCh);
            if nCh==nChParms
                [parmName value] = strread(dline, '%s %f');
                %[parmName value] = sscanf(dline, '%s %f');
                %[value] = strread(dline, '%f');
                parms(ind) = value;
            else
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
                %eval(Cstring)
            end
        else
            disp([parmName',' not found'])
        end
    end
end
fclose(fid);
%disp([parms(1),parms(2),parms(3),parms(4)])
Dt          = parms(1);
Dxz         = parms(2);
pvel        = parms(3);
svel        = parms(4);
fractx      = parms(5);
nWide       = parms(6);
iPlTerm     = parms(7);
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

