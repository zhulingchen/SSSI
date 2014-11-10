function [traces,texthead,binaryhead,extendedhead]=SEGY_read(filein, varargin)
% [traces,texthead,binaryhead,extendedhead]=SEGY_read(filein)
%
% SEGY_read reads in the entire SEG-Y file.
%
% Inputs:
%    filein= the name of the SEG-Y file that is to be read in.
%
% Outputs:
%    traces= a Trace object.  To get the trace header values use
%       SEGY_getHeader.  To get the trace data values use SEGY_getData.
%    texthead= a TextHeader object.  To get the character array use
%       SEGY_getHeader.
%    binaryhead= a BinaryHeader object. To get the binary header values use
%       SEGY_getHeader.  
%    extendedhead= a cell array of TextHeader Objects. To get the character array use
%       SEGY_getHeader.
%
% Optional Arguments 
% To only read in some traces:
%   [...]=SEGY_read(filein,'traces',{search parameters});
%      search parameters can be a series of header words and values as
%          defined in SEGY_getHeader
%
%
% Heather Lloyd 2010, Kevin Hall 2009, Chad Hogan 2004
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
%
if nargin<1 
    [filename, pathname, filterindex] = uigetfile( {'*.sgy;*.SEGY;*.segy;*.SGY', 'All SEG-Y Files (*.sgy, *.SGY, *.segy, *.SEGY)'; ...
        '*.*','All Files (*.*)'},'Select a SEG-Y File');
    if  filename~=0
    filein=[pathname,filename];
    else
        me=MException('SEGY_READ:InvalidFilename','Please select a valid Filename');
        throw(me)
    end
end
name={};
 for i = 1:2:length(varargin)
                name = varargin{i};
                value = varargin{i+1};
 end
if any(strcmpi(name,'traces'))
    segyobj=SegyFile(filein,'traces','1');
    segyobj.trace=segyobj.trace.getTraces(value{:}); 
else
    segyobj=SegyFile(filein);
end
 
texthead=segyobj.textheader;
binaryhead=segyobj.binaryheader;
traces=segyobj.trace;
extendedhead=segyobj.extendedheader;
end