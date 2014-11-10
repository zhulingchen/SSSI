function traj=matraj(amat,ksamps)
% traj=matraj(amat,ksamps)
%
% MATRAJ evaluates a matrix along a columnwise trajectory.
% That is, ksamps contains one entry per column of amat giving
% a desired sample number. Traj is a row vector, one entry
% per column of amat, containing those samples desired.
% If ksamps contains the entry 0, then zero will be returned
% in traj. No checking is done to ensure ksamps points to
% values between 1 and nrows. This must be done externally. 
%
% G.F. Margrave, University of Calgary, 1996
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
[nr,nc]=size(amat);
nsamps=nr*nc;
%compute indicies of first sample in each column
i1=1:nr:nr*nc-1;
%indicies of desired samples
isamp= i1+ksamps-1;
%initialize trajectory
traj=zeros(1,nc);
ind=find(ksamps~=0);
traj(ind)= amat(isamp(ind));
