function seisf=assemble_chunks(seischunk,fchunk,xchunk,ichunk,dx,nf,nx,vcrit)
% ASSEMBLE_CHUNKS ... Reassemble the complete wavefield from is spatially resampled chuncks
%
% seisf=assemble_chunks(seischunk,fchunk,xchunk,ichunk,dx,nf,nx,vcrit)
%
% seischunk ... cell array of the seismic data as partitioned into frequency chunks. 
%       Length(seischunck) is the number of chuncks. (These are usually
%       created by CREATE_CHUNKS). seisf=seischunk{j} is a seimic matrix for the 
%       jth frequency band.
% fchunk ... cell array of the frequency coordinate vectors for each chunk.
% xchunk ... cell array of the x coordinate vectors for each chunk. These
%       will each have a different sample rate. This means that
%       seisf=seischunk{j} is a matrix of size nf-by-nx where
%       nf=length(fchunk{j}) and nx=length(xchunk{j}).
% ichunk ... cell array of the frequency coordinate indicies of each chunk.
%       That is ifreqs=ichunk{k} gives the row numbers of the position of
%       seisf=seischunk{j} in the re-assembled matrix.
% dx ... desired spatial sample rate of the output
% nf ... total number of frequencies desried in the output
% nx ... total number of spatial samples in the output
% vcrit ... the 'critical' velocity that drives the spatial resampling. If
%       you have no better choice, make vcrit=min(v(:)), i.e. the slowest
%       velocity in your model.
% seisf ... re-assembled seismic matirx in the frequency domain. Size is nf-by-nx
% 
%
% G.F. Margrave, CREWES/POTSI 2004
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
nchunks=length(seischunk);
if(nargin<8)
    vcrit=ones(nchunks,1);%just make it a small number so that no effect
end
seisf=zeros(nf,nx);
for kchunk=1:nchunks
    xtmp=xchunk{kchunk};
    kcrit=fchunk{kchunk}/vcrit(kchunk);
    seisf(ichunk{kchunk},:)=kresample(seischunk{kchunk},xtmp,dx,kcrit);
end