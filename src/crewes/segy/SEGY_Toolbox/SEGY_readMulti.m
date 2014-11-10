function [traces,texthead,binaryhead,extendedhead]=SEGY_readMulti(file, varargin)
% [traces,texthead,binaryhead,extendedhead]=SEGY_readMulti({})
% [texthead,binaryhead,traces,extendedhead]=SEGY_readMulti({file1,file2,file3})
% [texthead,binaryhead,traces,extendedhead]=SEGY_readMulti({path})
% SEGY_readMulti reads in the entire first SEG-Y file and then only reads
%   in the traces from the other SEG-Y files and concatenates these traces
%   to the first trace.
%
%
% Inputs:
%    file= a cell array containing the names of the SEG-Y files that are  
%        to be read in.  file can also be the path to a directory
%        containing *.sgy files.  If file is set to {} then a gui prompt
%        will ask the user to select files.
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
%   [...]=SEGY_read({filein},'traces',{search parameters});
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
if nargin==0
    file=[];
end

if isempty(file) || ~iscell(file);
    [file,path]=uigetfile('*.sgy','Select SEG-Y Files to Concatenate','multiselect','on');
    if file==0
        disp('No Files Were Selected')
    end
end


if length(file)==1
    
if isdir(file{:})
    path=file{:};
    file={};
    list=ls(path);
    m=1;
    for k=1:size(list)
        dec=strfind(lower(list(k,:)),'.');
        if strcmpi('.sgy',strunpad(list(k,dec:end))) ||  strcmpi('.segy',strunpad(list(k,dec:end)))
            file{m}=strunpad(list(k,:));
            m=m+1;
        end
    end
end
end

file=sort(file);

[texthead,binaryhead,traces,extendedhead]=SEGY_read([path,file{1}], varargin{:});
if isempty(extendedhead)
    varargin{end+1}='extendedheader';
    varargin{end+1}='no';
else
    varargin{end+1}='extendedheader';
    varargin{end+1}='no';
    varargin{end+1}='hdroffset';
    varargin{end+1}=num2str(3600+length(extendedheader)*3200);
end
for k=2:length(file)
    sgytmp=SegyFile([path,file{k}], 'all','no',varargin{:});
    bintmp=sgytmp.binaryheader;trchtmp=sgytmp.traceheader;
    if SEGY_getHeader(bintmp,'jobid')==SEGY_getHeader(binaryhead,'jobid')
       if SEGY_getHeader(bintmp,'format')==SEGY_getHeader(binaryhead,'format') 
           trchtmp.filefmt=traces.traceheader.filefmt;
           trchtmp.machineformat=traces.traceheader.machineformat;
           trchtmp.tracetype=traces.traceheader.tracetype;
       end
        trctmp=Trace(trchtmp,'extendedheader','no');
        trctmp=trctmp.getTraces();
        if strcmpi(trctmp.traceheader.filefmt,traces.traceheader.filefmt);
            traces.traceheader.nontypecasthdr=[traces.traceheader.nontypecasthdr,trctmp.traceheader.nontypecasthdr];
        else
            
            words=trctmp.traceheader.definitions.values(:,strcmpi(trctmp.traceheader.definitions.keys(),'Name'));
            for k=1:length(words)
                vari=trctmp.traceheader.getheadervalue(words{k},0);
                vari=swapbytes(vari);
                trctmp.traceheader=trctmp.traceheader.setheadervalue(words{k},vari,0);
            end
        end
    end
    
    traces.tracedata.data=[traces.tracedata.data,trctmp.tracedata.data];
end
    
end