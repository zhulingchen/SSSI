function hdr=SEGY_getHeader(obj,word,scalevals)
% hdr=SEGY_getHeader(Texthead,'header')
% val=SEGY_getHeader(Binaryhead,word)
% val=SEGY_getHeader(Tracehead,word)
% val=SEGY_getHeader(Trace,word)
%
% If obj is a TextHeader Object, a 40 x 80 char array will be returned;
%
% If obj is a BinaryHeader Object word must be a name from the definitions
%   file.  SEG-Y revision 1 standard words are:
%     'jobid','lino','reno','ntrpr','nart','hdt','dto','hns','nso',
%     'format','fold','tsort','vscode','hsfs','hsfe','hslen','hstyp',
%     'schn','hstas','hstae','htatyp','hcorr','bgrcv','rcvm','mfeet',
%     'polyt','vpol','rev','flen','netfh','unasn1', through 'unasn84'.  
%   It will return the value associated with this word.
%
% If obj is a TraceHeader Object or a Trace object word must be a name 
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
%     'fbpicks','scalfb';
%   It will return the values associated with this word.  One value for
%   every trace.
% 
% For Trace , TraceHeader, or BinaryHeader objects the values will be
% automatically scaled by the scale factor in the header,  if it is desired
% to not have these values scaled then scalevals should be zero.  An
% example of this is seen below.
%
%    val=SEGY_getHeader(tracehead,'fbpicks',0);
%
% To allow the values to be scaled enter either of the following:
%
%    val=SEGY_getHeader(tracehead,'fbpicks',1);
%                     or
%    val=SEGY_getHeader(tracehead,'fbpicks');
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
    hdr=obj.header;
    return
end

if nargin<3
    scalevals=1;
end

if isa(obj,'Trace');
    obj=obj.traceheader;
end
if nargin<2
    me=MException('SEGY_getHeader:InsufficentInputs',....
        'word must be one of the names in the definitions file');
    throw(me)
end
if isa(obj,'Header')
    hdr=obj.getheadervalue(word,scalevals);
    hdr=double(hdr);
end

catch me
    error(me.message);
end












