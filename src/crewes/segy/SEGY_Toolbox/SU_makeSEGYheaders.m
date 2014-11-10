function [texthead, binhead]=SU_makeSEGYheaders(filein)
% [texthead,binaryhead]=SU_makeSEGYheaders(filein)
%
% SU_makeSEGYheaders creates a textheader and binary header that can be
%   used with the SEGY_Toolbox.  This function creates a textheader that is
%   a 40 x 80 character array conformining to SEG-Y Revision 1 standards as
%   well as contains information that the data in the file has been created
%   using Elastic Finite Difference Software.  This function creates a 
%   numerical binary header using data in the first traceheader that 
%   conforms to the SEG-Y Revision 1 standards.
%
% Inputs:
%   filein=the file name of the SU file containing at least one trace and
%     traceheader of data.
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

%% CREATE THE OBJECTS
    texthead=TextHeader(filein,'new','1');
    binhead=BinaryHeader(filein,'new','1','machineformat','ieee-le');
    tracehead=TraceHeader(filein,'extendedheaderflag','no','hdroffset','0',...
        'new','1','machineformat','ieee-le');
    fseek(tracehead.fid,0,'bof');
    tracehead.nontypecasthdr=fread(tracehead.fid,tracehead.hdrsize,'*uint8');
    tracehead.filefmt='L';
%% TEXTHEADER
% SEG-Y rev 1
%      12345678901234567890123456789012345678901234567890123456789012345678901234567890
C1  = 'C 1 CLIENT                        COMPANY  CREWES               CREW NO         ';
C2  = 'C 2 LINE            AREA                        MAP ID                          ';
C3  = 'C 3 REEL NO           DAY-START OF REEL     YEAR      OBSERVER                  ';
C4  = 'C 4 INSTRUMENT: MFG            MODEL            SERIAL NO                       ';
C5  = 'C 5 DATA TRACES/RECORD        AUXILIARY TRACES/RECORD         CDP FOLD          ';
C6  = 'C 6 SAMPLE INTERVAL         SAMPLES/TRACE       BITS/IN      BYTES/SAMPLE       ';
C7  = 'C 7 RECORDING FORMAT        FORMAT THIS REEL        MEASUREMENT SYSTEM          ';
C8  = 'C 8 SAMPLE CODE: FLOATING PT     FIXED PT     FIXED PT-GAIN     CORRELATED      ';
C9  = 'C 9 GAIN  TYPE: FIXED     BINARY     FLOATING POINT     OTHER                   ';
C10 = 'C10 FILTERS: ALIAS     HZ  NOTCH     HZ  BAND     -     HZ  SLOPE    -    DB/OCT';
C11 = 'C11 SOURCE: TYPE            NUMBER/POINT        POINT INTERVAL                  ';
C12 = 'C12     PATTERN:                           LENGTH        WIDTH                  ';
C13 = 'C13 SWEEP: START     HZ  END     HZ  LENGTH      MS  CHANNEL NO     TYPE        ';
C14 = 'C14 TAPER: START LENGTH       MS  END LENGTH       MS  TYPE                     ';
C15 = 'C15 SPREAD: OFFSET        MAX DISTANCE        GROUP INTERVAL                    ';
C16 = 'C16 GEOPHONES: PER GROUP     SPACING     FREQUENCY     MFG          MODEL       ';
C17 = 'C17     PATTERN:                           LENGTH        WIDTH                  ';
C18 = 'C18 TRACES SORTED BY: RECORD     CDP     OTHER                                  ';
C19 = 'C19 AMPLITUDE RECOVERY: NONE      SPHERICAL DIV       AGC    OTHER              ';
C20 = 'C20 MAP PROJECTION                      ZONE ID       COORDINATE UNITS          ';
C21 = 'C21 PROCESSING:                                                                 ';
C22 = 'C22 PROCESSING:                                                                 ';
C23 = 'C23                                                                             ';
C24 = 'C24 THE DATA IN THIS FILE HAS BEEN CREATED BY ELASTIC FINITE DIFFERENCING       ';
C25 = 'C25 SOFTWARE                                                                    ';
C26 = 'C26                                                                             ';
C27 = 'C27 CREWES USES BIG ENDIAN IEEE FLOATING POINT FORMATING FOR TRACES IN THIS FILE';
C28 = 'C28 TRACEHEADER BYTES 233-236 IS USED FOR FIRST BREAK PICKS                     ';
C29 = 'C29 TRACEHEADER BYTES 237-240 IS USED FOR FIRST BREAK SCALING FACTOR            ';
C30 = 'C30                                                                             ';
C31 = 'C31                                                                             ';
C32 = 'C32                                                                             ';
C33 = 'C33                                                                             ';
C34 = 'C34                                                                             ';
C35 = 'C35                                                                             ';
C36 = 'C36                                                                             ';
C37 = 'C37                                                                             ';
C38 = 'C38                                                                             ';
C39 = 'C39 SEG Y REV1                                                                  ';
C40 = 'C40 END TEXTUAL HEADER                                                          ';

%Combine lines to create textual header

texthead.header = [C1; C2; C3; C4; C5; C6; C7; C8; C9; C10; C11; C12; C13; C14;...
    C15; C16; C17; C18; C19; C20; C21; C22; C23; C24; C25; C26; C27;...
    C28; C29; C30; C31; C32; C33; C34; C35; C36; C37; C38; C39; C40];
texthead.format='ascii';
%% BINARYHEADER
k=dir(filein);
ns=SEGY_getHeader(tracehead,'ns');
ep=SEGY_getHeader(tracehead,'ep');
dt=SEGY_getHeader(tracehead,'dt');
ntrpr=(k.bytes/(4*double(ns)+double(tracehead.hdrsize)))/double(ep);
binhead=SEGY_setHeader(binhead,'ntrpr',ntrpr);
binhead=SEGY_setHeader(binhead,'hdt',dt);
binhead=SEGY_setHeader(binhead,'hns',ns);
binhead=SEGY_setHeader(binhead,'format',5);
binhead=SEGY_setHeader(binhead,'fold',ep);
binhead=SEGY_setHeader(binhead,'tsort',1);
binhead=SEGY_setHeader(binhead,'mfeet',1);
binhead=SEGY_setHeader(binhead,'rev',1);
binhead=SEGY_setHeader(binhead,'flen',1);
binhead=SEGY_setHeader(binhead,'netfh',0);


end