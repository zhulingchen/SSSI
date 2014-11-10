function [xwid,zwid]=resosim(x,t,vo,c,f,alpha,flag,xs,ts)
% RESOSIM: rake sesolution estimates for v(z) 
%
% [xwid,zwid]=resosim(x,t,vo,c,f,alpha,flag,xs,ts)
% 
% resosim estimates the resolution properties of a zero offset section.
% It analyzes aperture, record length, and spatial aliasing limits in
% the case velocity linear with depth and makes a display.
%
% Note: spatial aliasing limits are currently not functional
%
% x ... vector of x coordinates for seismic line
% t ... vector of t coordinates for seismic line
% vo ... constant velocity
% c ... accelaertor in velocity function: v(z) = vo + c*z
% f ... freqeuncy for aliasing computations
% alpha ... resolution scaling factor
%    ********* default = 2 ***********
% flag ... 0 -> produce a time display output
%          1 -> produce a depth display output
%    ********* default = 0 ***********
% ts ... vector of times at which resolution is to be simulated
%    ********* default = 10 evenly spaced times *********
% xs ... vector of x coordinates at which resolution is to be simulated
%    ********* default = 10 evenly spaced positions ********
%
% Gary F. Margrave, 1997, The CREWES Project, The University of Calgary
%
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

if(nargin<6) 
	%resolution scale factor
	alpha = 2;
end
if(nargin<7) flag = 0; end

% migration algorithm limit
angmig= 90;

nt=length(t);
nx=length(x);

dx=x(2)-x(1);
xmin=x(1);
xmax=x(nx);

T=t(nt);
dt=t(2)-t(1);

if(nargin<7)
	ntpts=10;
	it=linspace(2,length(t)-1,ntpts);
	ts=(it-1)*dt;
else
	ntpts=length(ts);
	it=ts/dt+1;
end

if(nargin<8)
	nxpts=10;
	ix=linspace(2,length(x)-1,nxpts);
	xs=(ix-1)*dx;
else
	nxpts=length(xs);
	ix=xs/dx+1;
end


if( c~=0)
	z=vo*(exp(c*t/2)-1)/c;
	zmax= z(nt);
else
	z=vo*t/2;
	zmax= z(nt);
end

xwid=zeros(length(it),length(ix));
zwid=xwid;

for k1=1:ntpts
	tnow=t(it(k1));

	%current depth
	if( c~=0)
		znow= vo*(exp(c*tnow/2)-1)/c;
	else
		znow= vo*tnow/2;
	end
	
	
	%velocity
	vnow= vo+c*znow;
	
	%record length
	angrec= threc(T,vo,c,[0 znow]);
	angrec=angrec(2);
	
	%spatial aliasing
	angalias= thalias(dx,f,vo,c,znow);
	
	%vertical resolution
	delz = alpha*vnow/(4*f);
	if(delz>zmax)delz=zmax;end
	
	for k2=1:nxpts
		xnow=x(ix(k2));
		Aleft= abs(xmin-xnow);
		Aright= abs(xmax-xnow);
		
		angleft= thaper(Aleft,vo,c,znow);
		angright=thaper(Aright,vo,c,znow);
	
		%angmin= -min([angrec angalias angleft angmig]);
		%angmax= min([angrec angalias angright angmig]);
		angmin= -min([angrec angleft angmig]);
		angmax= min([angrec angright angmig]);
		
		angband= angmax-angmin;
		angnot=(angmax+angmin)/2;
		
		% Postulate: delx is determined by angband/2
		delx = alpha*vnow/(4*f*sin(pi*angband/(2*180)));
		if(delx>xmax)delx=xmax;end
		
		%draw some little lines
		if (k1+k2==2) figure; end
		sn=sin(pi*angnot/180);
		cs=cos(pi*angnot/180);
		
		if(flag)
			factor=1;
		else
			factor=tnow/znow;
		end
		
		rotmat=[cs -sn; sn cs];
		p1p=[-delx/2 0]';
		p2p=[delx/2 0]';
		p3p=[0 -delz/2]';
		p4p=[0 delz/2]';
		p1=rotmat*p1p;
		p2=rotmat*p2p;
		p3=rotmat*p3p;
		p4=rotmat*p4p;
		
		p5p=[-delx/2 -delz/2]';
		p6p=[-delx/2 delz/2]';
		p7p=[delx/2 delz/2]';
		p8p=[delx/2 -delz/2]';
		p5=rotmat*p5p;
		p6=rotmat*p6p;
		p7=rotmat*p7p;
		p8=rotmat*p8p;
		
		%line([xnow-delx*cs/2, xnow+delx*cs/2],...
		%	factor*[znow-delx*sn/2,znow+delx*sn/2],'color','w');
		%line([xnow-delz*sn/2, xnow+delz*sn/2],...
		%	factor*[znow+delz*cs/2,znow-delz*cs/2],'color','w');
		%patch([xnow-delx*cs/2 xnow-delz*sn/2 xnow+delx*cs/2  xnow+delz*sn/2],...
		%	factor*[znow-delx*sn/2 znow+delz*cs/2 znow+delx*sn/2 znow-delz*cs/2],...
		%	'w');
		%patch(xnow+[p1(1) p3(1) p2(1) p4(1)],factor*(znow+[p1(2) p3(2) p2(2) p4(2)]),'w');
		patch(xnow+[p5(1) p6(1) p7(1) p8(1)],factor*(znow+[p5(2) p6(2) p7(2) p8(2)]),'w');
		zwid(k1,k2)=delz;
		xwid(k1,k2)=delx;
		
	end
end

if(flag)
	set(gca,'xlim',[0 xmax],'ylim',[0 zmax]);
else
	set(gca,'xlim',[0 xmax],'ylim',[0 T]);
end
flipy;
whitefig;
grid;
		
		
		
		
		
