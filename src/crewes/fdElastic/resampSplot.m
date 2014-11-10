function [factt,indX,Xarr,tifFile] = ....
        resampSplot(fact,spMin,nSpf,Dxz,kOut)
% function [factt,indX,Xarr,tifFile] = ....
%         resampSplot(fact,spMin,nSpf,Dxz,kOut)
    %Change trace selection parameters from a menu here,
        % set resample vectors
    % Mostly copied from resampSEGY
%The input parameters are
    %fact .... The default amplitude factor
    %spMin ... The offset or depth of the first available trace
    %nSpf .... The number of traces available
    %Dxz ..... Spatial sample rate
    %kOut .... The switch that gives the type of trace output required
%The output parameters are
    %factt ... The chosen amplitude factor
    %indX .... The index array of the chosen traces
    %Xarr .... The offset or depth array of the chosen traces
    %tifFile . The chosen unique part of a tif file name
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
%DtPhR = Dt*1000;
%tZero = tMin; phSp = Dxz; nPh = nSpf; DtPh = DtPhR; %phMax = spMin+nSpf*Dxz;
phSp = Dxz; nPh = nSpf; factt = fact; %phMax = spMin+nSpf*Dxz;
phMin = spMin; spMax = spMin+Dxz*nSpf; nPhMax = nSpf; phMax = spMax;
%factt = fact;
if kOut==1
    nPh = 50;   %Default number of WVA traces
    phSp = (phMax-phMin)/nPh;
end
  tifFile = '';
  str1 = 'Change seisplot';
  %str6 = ['First time = ',num2str(tZero)];
  %str7 = ['Delta t (ms) = ',num2str(DtPhR)];
  str8 = ['File name = ',tifFile];
  str99 = 'All OK';
kOptn = 0;
while kOptn < 7
    str3 = ['First trace = ',num2str(phMin)];
    str4 = ['Phone space. = ',num2str(phSp)];
    str5 = ['No. of traces. = ',num2str(nPh)];
    str6 = ['Last trace = ',num2str(phMax)];
    str7 = ['Amplitude = [',num2str(factt),']'];
    kOptn = menu(str1,str3,str4,str5,str6,str7,str8,str99);
    if kOptn==1
        phMin = input('Type first position ');              %First
        %str3 = ['First phone = ',num2str(phMin)];
        if (phMin<spMin)
            phMin = spMin;
        end
        %nPhMax = floor((phMax-phMin)/phSp)+1;
        nPh = floor((phMax-phMin)/phSp)+1;
    end
    if kOptn==2
        phSp = input('Type trace spacing ');                %Spacing
         %nPhMax = floor((phMax-phMin)/phSp)+1;           %- limits number
        nPh = floor((phMax-phMin)/phSp)+1;                  %- changes number
    end
    if kOptn==3
        nPh = input('Type no. of traces ');                 %No of traces
        nPhMax = floor((phMax-phMin)/Dxz)+1;                %- sets spacing
        if nPh > nPhMax
            nPh = nPhMax;
        end
        phSp = (phMax-phMin)/nPh;
    end
    if kOptn==4
        phMax = input('Type last position ');               %Last
        if (phMax>spMax)
            phMax = spMax;
        end
        %nPhMax = floor((phMax-phMin)/phSp)+1;
        nPh = floor((phMax-phMin)/phSp)+1;
    end
    if kOptn==5
        factt = input('Type amplitude ');                   %Amp
    end
    if kOptn==6
        tifFile = input('Name a tiff file ','s');           %tiff
        str8 = ['File name = ',tifFile];
    end
    if nPh>nPhMax 
        nPh = nPhMax;
    end
end
    %disp([phMin nPh phSp])
indX = zeros(nPh,1);    %Indices to select in X
Xarr = zeros(nPh,1);
for iPh = 1:nPh
    pos = phMin+(iPh-1)*phSp;
    indX(iPh) = round((pos-spMin)/Dxz)+1;
    %disp(indX(iPh))
    Xarr(iPh) = pos;
end
% tMax = tMin+nStep*DtPhR;
% nTsamp = round((tMax-tZero)/DtPh);
% indT = zeros(nTsamp,1);    %Indices to select in time
% Tarr = zeros(nTsamp,1);
% for iTime = 1:nTsamp
%     time = tZero+(iTime-1)*DtPh;
%     indT(iTime) = round((time-tMin)/DtPhR)+1;
%     Tarr(iTime) = time;
% end
%disp(Xarr)
