function plotcolourk(Ux,Uz,frClip,xvect,zvect,Dxz,ibar,iblack)
%
% Plot vector data with direction coded around the spectrum.
%function plotcolourk(Ux,Uz,frClip,xvect,zvect,Dxz,ibar,iblack)
%
%Ux  ...... %Matrix of displacements in the x direction
%Uz  ...... %Matrix of displacements in the z direction
%frClip ... %Reduces amplitudes higher than frClip to frClip level,
                %where the maximum amplitude is translated to 1.
%xvect  ... %X-coordinate values of the Ux and Uz arrays
%zvect  ... %Z-coordinate values of the Ux and Uz arrays
%Dxz  ..... %The sample rate in space
%ibar  .... %Control whether to plot colour key, and at top or bottom
%iblack  .. %Control for colour of the plot background
%
% P.M. Manning, Oct 2008
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

if nargin<7
  ibar = 0;     %0 - no colourbar  >0 - bottom,  <0 - top
end
if nargin<8
  iblack = 0;   %<=0 - white background,  >0 - black background
end

if frClip<=0
  frClip = 1;
end
if frClip>1
  frClip = 1;
end
sz = size(Ux);
nx = sz(1);
nz = sz(2);
Arr = zeros(nz,nx,3);
amp = sqrt(Ux.^2+Uz.^2);
maxa = max(max(abs(amp)));
%Allow for clipping at fraction frClip of the maximum
maxnew = maxa*frClip;
if maxnew<=0
  maxnew = 1;
end
ampnew = amp;
Uxnew = Ux;
Uznew = Uz;
if frClip<1
  large = find(abs(amp)>maxnew);
  ampnew(large) = maxnew;
  Uxnew(large) = Ux(large)*maxnew./amp(large);  %Reduce clipped amplitudes in
  Uznew(large) = Uz(large)*maxnew./amp(large);  % proportion for Ux and Uz
end
r3o2 = sqrt(3)/2;

%include colour ring for a key
if ibar~=0
  nuse = min([nx,nz]);
  nbx = 2*round(nuse/20)+1;
  ibx = nx+1-nbx;
  nbz = nbx;
  if ibar<1
    ibz = 1;
    lbz = nbz;
  else
    %ibz = nz+1-nbz;
    ibz = nz+1-nbz;
    lbz = nz;
  end
  
  ixmid = (nbx+1)/2;
  izmid = ixmid;
  maxr2 = izmid^2;
  Uxmat = zeros(nbx,nbz);
  Uzmat = zeros(nbx,nbz);
  for iz=1:nbz
    for ix=1:nbx
        Uzh = iz-izmid;
        Uxh = ix-ixmid;
        rad2 = Uzh^2+Uxh^2;
        if rad2<maxr2
            %rad = sqrt(rad2);xmat,Uzmat,1,1:10,1:10
            Uxmat(ix,iz) = Uxh/ixmid;
            Uzmat(ix,iz) = Uzh/izmid;
        end
    end
  end

  ampmat = sqrt(Uxmat.^2+Uzmat.^2);
  maxmat = max(max(abs(ampmat)));
  Uxnew(ibx:nx,ibz:lbz) = Uxmat*maxnew/maxmat;
  Uznew(ibx:nx,ibz:lbz) = Uzmat*maxnew/maxmat;
  ampnew = sqrt(Uxnew.^2+Uznew.^2);
end

if iblack<=0
  Arr(:,:,1)=1-(ampnew'-Uznew')/(maxnew*2);
  Arr(:,:,2)=1-(ampnew'+(.5*Uznew'-r3o2*Uxnew'))/(maxnew*2);
  Arr(:,:,3)=1-(ampnew'+(.5*Uznew'+r3o2*Uxnew'))/(maxnew*2);
else
  Arr(:,:,1)=(ampnew'-Uznew')/(maxnew*2);
  Arr(:,:,2)=(ampnew'+(.5*Uznew'-r3o2*Uxnew'))/(maxnew*2);
  Arr(:,:,3)=(ampnew'+(.5*Uznew'+r3o2*Uxnew'))/(maxnew*2);
end

less0 = find(Arr<0);
Arr(less0) = 0;
%Arr(:,:,1)=uint8(-(Uz'+amp)/(maxa*2));
%Arr(:,:,2)=uint8(-(amp-(.5*Uz'+r3o2*Ux'))/(maxa*2));
%Arr(:,:,3)=uint8(-(amp-(.5*Uz'-r3o2*Ux'))/(maxa*2));
image(xvect,zvect,Arr)

if ibar~=0      % Draw arrows
    hold on
    aleng = Dxz*ixmid/2;
    ixmid = ixmid+ibx-1;
    izmid = izmid+ibz-1;
    xv = zeros(1,8);
    zv = zeros(1,8);
    uv = zeros(1,8);
    vv = zeros(1,8);
    anginc = pi/4;
    ang = 0;
    xmid = (ixmid-1)*Dxz+xvect(1);
    zmid = (izmid-1)*Dxz+zvect(1);
    for iang=1:8
        xv(iang) = cos(ang)*aleng/2+xmid;
        zv(iang) = sin(ang)*aleng/2+zmid;
        uv(iang) = cos(ang)*aleng;
        vv(iang) = sin(ang)*aleng;
        ang = ang+anginc;
    end
    quiver(xv,zv,uv,vv,0,'k')
    hold off
end
%plotcolourk(ampx(:,:,1),ampz(:,:,1),clip,0:Dxz:(lx-1)*Dxz,0:Dxz:(lz-1)*Dxz,Dxz,1)