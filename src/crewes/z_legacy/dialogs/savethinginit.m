function savethinginit(transfer,q1,name)
% savethinginit(transfer,q1,name)
%
% SAVETHINGINIT initiates a dialog which allows the user to decide whether to
% save something.
% transfer = a string matlab command to be evaluated with EVAL at the 
%            completion of the dialog. Usually this will re-invoke the 
%            calling program
% q1 = question to appear at the top of the dialog. Usually this will be 
%      something like 'Save the Universe first?'
% name = current name by which the thing to be save is known. User will have
%        a chance to change this
%
% G.F. Margrave February 1994
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
l=max([length(transfer) length(q1) length(name)]);
qs=ones(3,l);
qs(1,1:length(q1))=q1;
qs(2,1:length(name))=name;
qs(3,1:length(transfer))=transfer;
set(gca,'userdata',qs);
savething('init');
