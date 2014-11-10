function obj=SEGY_setHeader(obj,word,value,scalevals)
% Texthead=SEGY_setHeader(Texthead,'header',header)
% Binaryhead=SEGY_setHeader(Binaryhead,word,value)
% Tracehead=SEGY_setHeader(Tracehead,word,value)
% Trace=SEGY_setHeader(Trace,word,value)
%
% For a TextHeader Object word must be 'header' and value must be a a 
%   40 x 80 char array.
%
% For a BinaryHeader Object word must be a name from the definitions
%   file.  SEG-Y revision 1 standard words are:
%     'jobid','lino','reno','ntrpr','nart','hdt','dto','hns','nso',
%     'format','fold','tsort','vscode','hsfs','hsfe','hslen','hstyp',
%     'schn','hstas','hstae','htatyp','hcorr','bgrcv','rcvm','mfeet',
%     'polyt','vpol','rev','flen','netfh','unasn1', through 'unasn84'.  
%   Value must be a single value that belonging to the word.
%
% For a TraceHeader Object or a Trace object word must be a name 
%   from the definitions file.  SEG-Y revision 1 standard words are:
%     'tracl','tracr','fldr','tracf','ep','cdp','cdpt','trid','nvs','nhs',
%     'duse','offset','gelev','selev','sdepth','gdel','sdel','swdep',
%     'gwdep','scalel','scalco','sx','sy','gx','gy','counit','wevel',
%     'swevel','sut','gut','sstat','gstat','tstat','laga','lagb','delrt',
%     'muts','mute','ns','dt','gain','igc','igi','corr','sfs','sfe','slen',
%     'styp','stas','stae','tatype','afilf','afils','nofilf','nofils',
%     'lcf','hcf','lcs','hcs','year','day','hour','minute','sec','timbas',
%     'trwf','grnors','grnofr','grnlof','gaps','otrav','cdpx','cdpy',
%     'iline','xline','sp','scalsp','tval','tconstm','tconste','tunit',
%     'devtr','scalt','stypeo','sedm','sede','smmtm','smmte','smmtu',
%     'fbpicks','scalfb';}
%   Value must be a numerical array where the values are described by the 
%   word and the length of value must be the same as the number of traces.
%
% For Trace , TraceHeader, or BinaryHeader objects the values will be
% automatically scaled by the scale factor in the header,  if it is desired
% to not have these values scaled then scalevals should be zero.  An
% example of this is seen below.
%
%    val=SEGY_setHeader(tracehead,'fbpicks',values,0);
%
% To allow the values to be scaled enter either of the following:
%
%    val=SEGY_setHeader(tracehead,'fbpicks',values,1);
%                     or
%    val=SEGY_setHeader(tracehead,'fbpicks',values,);
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

try
   

if isa(obj,'TextHeader');
    if ischar(value) && size(value)==[40,80]
    obj.header=value;
    return
    else
    me=MException('SEGY_getHeader:InsufficentInput',....
        'If obj is a TextHeader header must be a 40 x 80 character array.');
    throw(me) 
    end
end

 if nargin<2
    me=MException('SEGY_getHeader:InsufficentInputs',....
        'word must be one of the names in the definitions file');
    throw(me)
end
if nargin<3
    me=MException('SEGY_getHeader:InsufficentInputs',....
        'Value must be a numerical input if obj is a Traceheader, Trace or Binaryheader');
    throw(me)
end
if nargin<4
    scalevals=1;
end

if isa(obj,'Trace');
    obj.traceheader=obj.traceheader.setheadervalue(word,value,scalevals);
    return
end

if isa(obj,'Header')
    obj=obj.setheadervalue(word,value,scalevals);
end

catch me
    error(me.message);
end












