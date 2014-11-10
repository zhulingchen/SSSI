function [tops,ztops]=readtops(filename)

% [tops,ztops]=readtops(filename)
%
% Read in tops information from a standard ascii file. The file is
% defined as one with an arbitrary number of header lines containing any
% information followed by two columns of formation top name and depth.
% The header lines must begin with two slashes (//)
% but are otherwise arbitrary and are ignored by this function. 
% The top names may be arbitrary character strings (including embedded
% whitespace) while the top depths are taken to be the first numerical value
% surrounded by while space. 
%
% Example tops file:
%
% // this is a comment
% top bun     1392.3
% ketchup     1512.0
% mustard     1512.0
% relish      1565.0
% burger 1    1572.0
% burger 2    1682.0
% bottom bun  1.792e3
%
% G.F. Margrave November 1994, updated 2003.
% Improved by Henry Bland at some later date.
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

fid=fopen(filename);

if(fid==-1)
	tops=-1;
	ztops=[];
	return;
end

%allocate some matricies
maxtops=100;
tops=setstr(32*ones(maxtops,30));
ztops=zeros(maxtops,1);

maxname=0;
ntops=0;

%read through the file skipping any // lines
buf=fgetl(fid);
while ischar(buf)
    if isempty(regexp(buf,'^\s*//'))
        % Find the start end end of the last number on the line
        % Number may contain digits minuses, the letter 'e' for exponents.
        % The trailing dollar sign ensures that this is found at the end of
        % the line (possible whitespace is ok between number and end of
        % line.
%        [startOfNum, endOfNum] = regexp(buf,'[0-9.e---]+\s*$');
        [startOfNum, endOfNum] = regexp(buf,'[-+]?([0-9]*\.)?[0-9]+([eE][-+]?[0-9]+)?\s*$');
        if (~isempty(startOfNum)) && (startOfNum(end) >= 3)
            depth = str2num(buf(startOfNum(end):endOfNum(end)));
            topName = trimWhitespace(buf(1:startOfNum(end)-1));
            if ~isempty(topName)
                ntops=ntops+1;
                if(ntops>maxtops)
                    tops=[tops;setstr(32*ones(10,30))];
                    ztops=[ztops;zeros(10,1)]
                    maxtops=maxtops+10;
                end
                ztops(ntops)=depth;
                if(length(topName)>maxname)
                    maxname=length(topName);
                end
                tops(ntops,1:length(topName))=topName;
            end
        end
    end
    buf=fgetl(fid);
end

tops=tops(1:ntops,1:maxname);
ztops=ztops(1:ntops);

fclose(fid);

function noWhitespace=trimWhitespace(t)
    t = regexprep(t,'^\s+','');
    t = regexprep(t,'\s+$','');
    noWhitespace = t;    
    
