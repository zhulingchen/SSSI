function fdSaveRetr2(modelDir,SorR,permName)
%fdSaveRetr(SorR,permName)
%Save trace files to longer term names, or retrieve
%
%SorR is a switch to indicate whether to save or retrieve trace data
    %The valid values are 's' or 'r'
%permName will become the root name of the (permenant) trace files
%         or the permenant name of the trace data to retrieve
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
if SorR == 'r'
    disp('Retrieve')
else
    disp('Save')
end
% pref = 'Name';
% source(1,:) = 'wellFileX.mat'; dest(1,:) = 'wFX.fdt';
% source(2,:) = 'wellFileZ.mat'; dest(2,:) = 'wFZ.fdt';
% source(3,:) = 'surfFileX.mat'; dest(3,:) = 'sFX.fdt';
% source(4,:) = 'surfFileZ.mat'; dest(4,:) = 'sFZ.fdt';
% source(5,:) = 'headrFile.mat'; dest(5,:) = 'hed.fdt';
source(1,:) = [modelDir,'\wellFileX.mat'];
disp(source(1,:))
source(2,:) = [modelDir,'\wellFileZ.mat'];
source(3,:) = [modelDir,'\surfFileX.mat'];
source(4,:) = [modelDir,'\surfFileZ.mat'];
source(5,:) = [modelDir,'\headrFile.mat'];
% dest(1,:) = [modelDir,'\',permName,'_wFX.fdt'];
% disp(dest(1,:))
% dest(2,:) = [modelDir,'\',permName,'_wFZ.fdt'];
% dest(3,:) = [modelDir,'\',permName,'_sFX.fdt'];
% dest(4,:) = [modelDir,'\',permName,'_sFZ.fdt'];
% dest(5,:) = [modelDir,'\',permName,'_hed.fdt'];
dest(1,:) = [permName,'_wFX.fdt'];
disp(dest(1,:))
dest(2,:) = [permName,'_wFZ.fdt'];
dest(3,:) = [permName,'_sFX.fdt'];
dest(4,:) = [permName,'_sFZ.fdt'];
dest(5,:) = [permName,'_hed.fdt'];

for ii=1:5
    %%destA = [pref,dest(ii,:)];
    %destA = [permName,dest(ii,:)];
    %disp(destA)
    if(strcmpi('r',SorR))
        copyfile(dest(ii,:),source(ii,:));
    else
        copyfile(source(ii,:),dest(ii,:));
    end
end

