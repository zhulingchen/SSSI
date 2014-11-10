function [grout, x, y] = ifftgrid(grin, kx, ky, nxout, nyout, pctaper)
% [grout, x, y] = ifftgrid(grin, kx, ky, nxout, nyout, pctaper)
% [grout, x, y] = ifftgrid(grin, kx, ky, nxout, nyout)
% [grout, x, y] = ifftgrid(grin, kx, ky, matout)
% [grout, x, y] = ifftgrid(grin, kx, ky)
%  
% ifftgrid computes the 2-D inverse fft of a complex valued grid. It is assumed
% that the spectral grid being inverse transformed was originally created
% in fftgrid. Thus there should be no need for padding and tapering prior
% to the inverse transform. (If this is not so, then perform these actions
% prior to ifftgrid.)
%
% grin = input complex grid (matrix)
% kx = vector of coordinates labeling the columns of grin (ordinarily created
%      in fftgrid)
%    !!! length(x) must equal the number of columns of grin !!!
% ky = vector of coordinates labeling the rows of grin
%    !!! length(x) must equal the number of rows of grin !!!
% nxout = desired number of columns of output. This is used to remove the
%    pad applied on the forward transform
%    *********** default nxout = number of columns of grin ********
% nyout = desired number of rows of output. This is used to remove the
%    pad applied on the forward transform
%    *********** default nyout = number of rows of grin ********
% matout = any matrix whose size is used to determine nxout and nyout
% pctaper = percentage oc cosine taper which was applied on the forward
% 	transform and is to be removed now.
%    ********** default = 0 ********** (means no taper removal)
%
% grout = real valued output grid created from the inverse transform. (Should
%	the inverse transform produce and non-zero imaginary part, it is
%	discarded.)
% x = row vector giving the x coordinates of the columns
%      of spec
% y = column vector giving the y coordinates of the rows 
%      of spec 
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
[rowin colin]=size(grin);
if( nargin < 6) pctaper = 0; end
if( nargin < 4 ) nxout = colin; end
if( nargin < 5 ) nyout = rowin; end
if( nargin == 4) [nxout, nyout]=size(nxout); end
% test for valid size of nxout & nyout
if( nxout > colin ) nxout = colin; end
if( nyout > rowin ) nyout = rowin; end
		
% compute the ifft and shift it
		grout = ifft2(fftshift(grin));
		grin=[];
		grout = real(grout(1:nyout,1:nxout));
% remove any taper that was applied
if( pctaper > 0 )
   stab = .0001;
   rowtaper = mwindow(nxout, pctaper);
   ind=find(rowtaper<stab);
   rowtaper(ind)=stab*ones(size(ind));% apply a stab factor
   rowtaper = ones(nyout,1)*rowtaper;
   coltaper = mwindow(nyout,pctaper)';
   ind=find(coltaper<stab);
   coltaper(ind)=stab*ones(size(ind));
   coltaper = coltaper*ones(1,nxout);
   coltaper=coltaper.*rowtaper;
   rowtaper=[];
   grout = grout./coltaper;
   coltaper=[];
end
		
% compute the output coordinate vectors
		xnyq = 1./( kx(2)-kx(1) );
		ynyq = 1./( ky(2)-ky(1) );
		dx = 1./(2.*kx(length(kx)));
		dy = 1./(2.*ky(length(ky)));
		x = xcoord(0.,dx,nxout);
		y = xcoord(0.,dy,nyout)';
