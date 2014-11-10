function [aryout]=conv45(aryin)
% CONV45: used by KIRK_MIG for 45 degree phase shift
%  CONV45 performs 45 degree phase shift for the
%         traces of the input dataset. The process
%         is fulfilled in time domain.
%
%  USAGE:
%
%  function [aryout]=conv45(aryin)
%
%  By Xinxiang Li.
%  CREWES Project, University of Calgary
%  1996
% 
%  The filter vector is from a relative
%  FORTRAN routine by Dr. John C. Bancroft.
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

itrans = 0;
[nrow,nvol]=size(aryin);
if nrow == 1
	aryin=aryin';
	nrow = nvol;
	nvol = 1;
	itrans = 1;
end
aryout=zeros(size(aryin));

filt = [-0.0010 -0.0030,-0.0066,-0.0085,-0.0060, -0.0083, -0.0107,...
-0.0164,-0.0103,-0.0194,-0.0221,-0.0705,0.0395,-0.2161,-0.3831,...
0.5451,0.4775,-0.1570,0.0130,0.0321,-0.0129]';

for j = 1:nvol
	conv1=conv(aryin(:,j),filt);
	aryout(:,j)=conv1(16:nrow+15);
end
if itrans
	aryout = aryout';
end
