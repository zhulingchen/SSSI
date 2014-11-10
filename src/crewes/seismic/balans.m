function trout=balans(trin,trref,tz)
% BALANS: match the rms power of one trace to another
%
% trout=balans(trin,trref,tz)
% trout=balans(trin,trref)
% trout=balans(trin)
%
% BALANS adjusts the rms power of trin to equal that of trref
%
% trin= input trace or gather to be balanced
% trref= input reference trace
% ******** default= ones(trin) *******
% tz = vector of indicies specifying a time zone over which balancing is to
%     be done.
% trout= balanced output trace
%
% NOTE: If trin is a matrix, then traces are assumed to be in the columns.
% Each column is balanced w.r.t. the reference trace. trref may not be a
% matrix.
%
% by G.F. Margrave, June 1991 and 2009
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


if nargin<2
 trref=ones(size(trin));
end
if nargin<3
   if(length(trref)>=length(trin))
        tz=[1:length(trin)];
   else
        tz=[1:length(trref)];
   end
end 
if(~isvector(trref))
        error('trref must be a vector')
end
if(isvector(trin))
    %single channe case
    trout= trin*norm(trref(tz),2)/norm(trin(tz),2);
else
    trout=zeros(size(trin));
    
    top=norm(trref(tz),2);
    for k=1:size(trin,2)
        trout(:,k)=trin(:,k)*top/norm(trin(tz,k),2);
    end
end
  
