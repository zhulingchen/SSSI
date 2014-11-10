function fact=fmplot(fmat,vaflag,fact,flipy,kolor,hax)

% fact=fmplot(fmat,vaflag,fact,flipy,kolor,hax)
%
% FMPLOT does a quick plot of a flexi-mat in a figure window (made by fmplot). 
% It plots the columns as wiggle traces centered at their x coordinates.
%
%	fmat ... the flexi-mat to be plotted
%	vaflag ... if 0, then traces are plotted at wt (wiggle traces) 
%		if 1, then the traces are plotted wtva
%      if 2 the rev. pol. wtva
%		************* default = 1 *********
%	fact ... scaling factors. Make fact(1) bigger for bigger wiggle
%		traces. fact(2) controls the overall plot scale. If not
%		provided, it is computed as max(abs(matrix)). To scale two
%		fmplots the same, capture the return value from the first and
%		provide is as fact(2) for the second
%		************* default 1.5 ***********
%	flipy ... if 1, then the y axis is reverse so that it goes from top
%		of window to the bottom
%	   ************* default 1 ***********
%  kolor ... color to plot the traces
%     ************* default = [1 0 0] (red) ************
%  hax ... handle of the axis to plot the fleximat in. If -1, then a new figure
%          is created.
%     ************* default = -1 **************
%Note: To scale two fmplots with respect to the maximum absolute value
% on the first plot, capture the return value from the first
% and provide it as the third argument for the second plot:
%	fact=fmplot(fmat1);
%	fmplot(fmat2,1,fact);
%
% G.F. Margrave, June 1994
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

if(nargin < 6)
	hax=-1;
end
if(nargin<5)
	kolor=[1 0 0];
end
if(nargin<4)
	flipy=1;
end
if(nargin<3)
	fact=1.5;
end
if(nargin<2)
	vaflag=1;
end
%unpack
mat=fmget(fmat,'mat');
x=fmget(fmat,'x');
y=fmget(fmat,'y');

if(length(x)>1)
	bnds=(max(x)-min(x))/(length(x)+1);
end

if( length(fact)<2 )
	ilive=find(~isnan(mat));
	trcmin=min(mat(ilive));
	trcmax=max(mat(ilive));
	s=max([abs(trcmax) abs(trcmin)]);
	fact=[fact s];
end

if( hax==-1 )
	figure;
else
	hfig=get(hax,'parent');
	figure(hfig);
	set(hfig,'currentaxes',hax);
end

if(flipy)
	set(gca,'ydir','reverse');
end

for k=1:length(x)
	trc=mat(:,k)/fact(2);
	ilive=find(~isnan(trc));
	m=mean(trc(ilive));
	trc=(trc-m)*bnds*fact(1) + x(k);
	if(~vaflag)
		line(trc(ilive),y(ilive),'color',kolor);
	elseif( vaflag==1)
		%replace NaN's with mean values
		if(length(ilive)< length(trc))
			idead=find(isnan(trc));
			trc(idead)=x(k)*ones(size(idead));
		end
		wtva(trc,y,kolor,x(k),1,1,1);
	elseif(vaflag==2)
		wtva(trc(ilive),y(ilive),kolor,x(k),-1,1,1);
	end
end
