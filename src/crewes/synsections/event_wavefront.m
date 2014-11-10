function amat=event_wavefront(amat,t,x,tnot,xnot,v,amp,aper)
% EVENT_WAVEFRONT: inserts a wavefront (circle) event in a matrix.
%
% amat=event_wavefront(amat,t,x,tnot,xnot,v,amp,flag,aper)
%
% EVENT_WAVEFRONT inserts a wavefront circle event in a matrix. This is
% done in the time domain so that the reciprocal nature of diffraction
% hyperbolae and wavefront circles can be easily demonstrated.
%
% amat ... the matrix of size nrows by ncols
% t ... vector of length nrows giving the matrix t coordinates
% x ... vector of length ncols giving the matrix x coordinates
% tnot ... t coordinate of the wavefront nadir (lowest point)
% xnot ... x coordinate of the wavefront nadir (lowest point)
% v ... velocity of wavefront. Value is divided by two to simulate post
%       stack.
% amp ... amplitude of the wavefront
% aper ... truncate wavefront beyond this aperture
% ******** default is no truncation ********
%
% G.F. Margrave, CREWES Project, University of Calgary, 2014
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

if( nargin <8 )
	aper = inf;
end

v=v/2;
r=v*tnot;%radius of the circle
r=min([r aper]);
x1=xnot-r;
x2=xnot+r;
ix=near(x,x1,x2);%these are the columns of interest

%loop over columns
[nsamp,nc]=size(amat);

dt=t(2)-t(1);
tmin=t(1);

for k=ix
    xoff=x(k)-xnot;
    tk = sqrt(r^2-xoff^2)/v;
    a=amp;
    ik=(tk-tmin)/dt+1;
    if( between(1,nsamp,ik) )
        ik1=floor(ik);
        ik2=ceil(ik);
        if(ik1==ik2)
            %exactly on a sample
            amat(ik1,k)=amat(ik1,k)+a;
        else
            %a simple interpolation
            amat(ik1,k)=amat(ik1,k)+a*(ik-ik2)/(ik1-ik2);
            amat(ik2,k)=amat(ik2,k)+a*(ik-ik1)/(ik2-ik1);
        end
    end
end
