function snaps=afd_makesnapshots(dx,dtstep,velocity,snap1,snap2,tsnaps,laplacian,boundary,wavelet)
% AFD_MAKESNAPSHOTS ... returns a cell array of snapshots of wave propagation
%
% snaps=afd_makesnapshots(dx,dtstep,velocity,snap1,snap2,tsnaps,laplacian,boundary,wavelet)
%
% AFD_MAKSNAPSHOTS propogates a wavefield forward in time saveing snapshots
% at specified times. This function simply calls afd_snapn to accomplish
% the propagation.
%
% delx = the horizontal AND vertical bin spacing in consistent units
% dtstep = time step in seconds
% velocity = the input velocity matrix in consisnent units
%          = has a size of floor(zmax/delx)+1 by floor(xmax/delx)+1
% snap1 = the wavefield at time=0 - delt (same size as velocity matrix)
% snap2 = the wavefield at time = 0 (same size as velocity matrix.
% NOTE: Commonly snap1 is zeros(size(velocity)) while snap2 is all zeros
%       except at source locations where is it usually 1.0.
% tsnaps =  vector of times at which snapshots will be returned. These must
%       be monotonically increasing.
% laplacian - an option between two approximation to the laplacian operator
%           - 1 is a 5 point approximation (2nd order)
%           - 2 is a nine point approximation (4th order)
% boundary = indicate whether all sides of the matrix are absorbing
%          = '1' indicates all four sides are absorbing
%          = '2' choses three sides to be absorbing, and the top one not to be
%             this enables sources to be put on the surface
% wavelet = wavelet specified as a time series. Must be sampled at delt and
%       must be as long as the desired simulation (i.e. the maximum of
%       tsnaps).
% NOTE: At each time step, the corresponding sample of the wavelet is
% injected into the simulation at the non-zero locations in snap2
%
% snaps = cell array of length(tsnaps) containing one snapshot for each
% time in tsnaps. All smapshots are indentical in size.
%
% by Gary Margrave, January 2014
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

test=diff(tsnaps);
ind=find(test<=0, 1);
if(~isempty(ind))
    error('tsnaps must be monotonically increasing')
end

%check for t=0 in tsnaps, we always include such
ind=find(tsnaps==0, 1);
if(isempty(ind))
    tsnaps=tsnaps(:)';%force row vector
    tsnaps=[0 tsnaps];
end


tin=tsnaps(1);
snapm=snap1;
snapn=snap2;
snaps=cell(size(tsnaps));
snaps{1}=snap2;
tprev=clock;
tused=0;
for k=2:length(tsnaps)
    tout=tsnaps(k);
    [snapn,snapm]=afd_snapn(dx,dtstep,velocity,snapm,snapn,snap2,tin,tout,laplacian,boundary,wavelet);
    snaps{k}=snapn;
    tin=tout;
    tnow=clock;
    tstep=etime(tnow,tprev);
    tused=tused+tstep;
    tneeded=(length(tsnaps)-1)*tstep;
    tremaining=tneeded-tused;
    disp(['Completed snapshot ' int2str(k) ' of ' int2str(length(tsnaps)) ...
        '. time used= ' int2str(tused) 's, time remaining = ' num2str(tremaining) ' s.'])
    tprev=tnow;
end