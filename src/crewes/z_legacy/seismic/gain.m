function trout=gain(trin,t,atten,tmin,tmax)
% trout=gain(trin,t,atten,tmin,tmax)
% trout=gain(trin,t,atten)
% gain(trin,t,0,tmin,tmax)
%
% GAIN applies exponential gain to the input trace, trin. 
%  If atten is negative, the exponential decay will result.
% GAIN with no return variable assigned creates an interactive
% analysis plot to allow the estimation of gain parameters. In this
% case, atten is not used.
%
% trin= input trace
% t= input time coordinate vector
% atten= exponentioal gain constant in db/sec
% tmax= time at which applied gain ceases to be exponential
%       and simply becomes constant
%   ********** default = t(length(t)) ***********
% tmin = start time for gain corrections
%   ********** default = 0.0 *********
% trout = output trace with gain applied
%
% by G.F. Margrave, May 1991
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
%test for dbanal curve
if(isstr(trin))
	if(strcmp(trin,'dbcurve') )
		dat=get(gcf,'userdata');
		if( isempty(dat))
			pt=get(gca,'currentpoint');
			pt=pt(1,1:2);
			set(gcf,'userdata',pt);
			title('Click again');
			return;
		end
		pt1x=dat(1); pt1y=dat(2);
		pt=get(gca,'currentpoint');
		pt=pt(1,1:2);
		m=(pt(2)-pt1y)/(pt(1)-pt1x);
		x=[pt1x pt(1)];
		y=[pt1y pt(2)];
		line('xdata',x,'ydata',y,'color','r','linestyle','-.');
		htxt=text(pt(1),pt(2),[num2str(abs(m)) ' db/sec']);
		set(htxt,'color','r');
		set(gcf,'userdata',[]);
		title('Click twice to define db/sec curve');
		xlabel([' Chosen atten = ' num2str(abs(m))])
		return;
	end
end
 
% set defaults
 if nargin<5
   tmax=t(length(t));
 end
 if nargin<4
   tmin=t(1);
 end
 if(tmin>=tmax)
		error(' tmin and tmax specified incorrectly ');
	end
 if(nargout==0) %create gain analysis plot
		trinew=padpow2(trin,0);
		env=abs(hilbm(trinew));
		env=env(1:length(trin));
		env=todb(env);
		figure;
		iz=near(t,tmin,tmax);
		plot(t,env);
		ylabel('Envelope Amplitude (db)');
		p=polyfit(t(iz),env(iz),1);
		fit=p(1)+p(2)*t(iz);
		hfit=line('xdata',t(iz),'ydata',fit,'color','r','linestyle','-.');
		htxt=text(t(iz(length(iz))),fit(length(fit)),...
			[num2str(abs(p(2))) ' db/sec']);
		set(htxt,'color','r');
		set(gcf,'windowbuttondownfcn','gain(''dbcurve'')');
		title('Click twice to define db/sec curve');
		xlabel(['Least sqs fit gives atten = ',num2str(abs(p(2)))]);
		set(gcf,'userdata',[]);
		grid
		return;
	end
% compute gain trace
 alpha=log(10)*atten/20.;
 iz=near(t,tmin,tmax);
 tgain=ones(size(t));
 tgain(iz)=exp(alpha*(t(iz)-tmin));
 tgain(iz(length(iz))+1:length(t))=tgain(iz(length(iz)))*...
		ones(size(iz(length(iz))+1:length(t)));
% apply gain
 trout=trin.*tgain;
