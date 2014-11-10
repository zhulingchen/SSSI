function obj=SEGY_getTraces(tracehead,varargin)
% traces=SEGY_getTraces(tracehead)
% traces=SEGY_getTraces(tracehead,word,value);
% traces=SEGY_getTraces(tracehead,word1,value1,word2,value2,word2,value3,...);
%
% SEGY_getTraces is a function that allows the user to only read in some of
% the traces.  This is done by comparing value to the traceheader data.
% word must be one of the names defined in the definitions file.  For SEG-Y
% revision 1 standards the word must be one of the following:
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
%
% It is possible to search using multiple words.  An example would be 
%  traces=SEGY_getTraces(tracehead,'gx',1:50,'cdpx',60);  
% this would search for traces that had a x coordinate from 1:50 and has a
% common depth point of 60.
%
% Input:
%  tracehead= can be either a TraceHeader object which can be created by 
%     SEGY_readHeaders or a filename.
%  word= a name defined in the definitions file, SEG-Y Revision 1 Standard
%     words are given above
%  value= a value is either scalar values or numerical arrays
%
% Output:
%  traces= a Trace object.  To get the traceheader values use
%       SEGY_getHeader.  To get the trace data values use SEGY_getData.
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

obj=Trace(tracehead);
obj=obj.getTraces(varargin);
end
