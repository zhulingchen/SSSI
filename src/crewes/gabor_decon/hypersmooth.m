function [tvshyp,smooth,smodev,tflevels]=hypersmooth(tvs,trow,fcol,nbins)
% HYPERSMOOTH: smooth a time variant spectrum along hyperbolic contours.
%
% [tvshyp,smooth,smodev,tf]= hypersmooth(tvs,trow,fcol,nbins)
%
% Hypersmooth defines nbins bins from 0 to tmax*fmax. Each sample in the
% input tvs is summed into one of thes bins depending on the time-frequency
% value of its coordinates. The number of samples summed into each bin is
% counted and the mean value of the tvs in each bin is then determined. The
% result is a 1D curve giving amplitude as a function of t*f and determined
% by nbins points. Finally, a tvs is populated by interpolating values off
% this curve.
%
% tvs ... time variant spectrum to be smoothed
% trow ... time coordinates for the rows of tvs
% fcol ... frequency coordinates for the columns of tvs
% nbins ... number of levels or bins in the smoothing
% tvshyp ... output hyperbolically smoothed tvs
% smooth ... vector of average values, one per bin
% smodev ... standard deviation for each bin
% tflevels ... time-frequency coordinates for smooth and smodev
%
% G.F. Margrave, CREWES, 2002
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

%nbins must not be greater than nf
nt=length(trow);
nf=length(fcol);
t1=trow(ceil(.05*nt));%pick a value 5 percent along the t axis
f1=fcol(ceil(.05*nf));%pick a similar f value
if(nbins>nf); nbins=nf; end
if(nargin<4); nbins=nf; end

%maximum tf value
tfmax=trow(end)*fcol(end);
deltf=tfmax/nbins;

%make a matrix of tf values
tf=trow(:)*(fcol(:)');

% tic
% tfnow=0;
% smooth=zeros(1,nbins+1);
% tflevels=zeros(1,nbins+1);
% tflevels(end)=tfmax;
% for k=1:nbins
%    tf1=tfnow+deltf;
%    tflevels(k)=(tfnow+tf1)/2;
%    ind=find((tf<=tf1)&(tf>tfnow));
%    smooth(k)=mean(tvs(ind));
%    tflevels(k)=tfnow;
%    tfnow=tfnow+deltf;
% end
% toc

%tfnow=0;
smooth=zeros(1,nbins);
fold=smooth;
small=1000*eps;
%tflevels=linspace(0,tfmax,nbins+1);
%choose tflevels logarithmically
dtflog=(log(tfmax)-log(t1*f1))/(nbins-1);
tfloglevels=log(t1*f1):dtflog:log(tfmax);
tflog=log(tf+100*small);

%this loop computes the average value in each bin
for k=1:nf
    for kk=1:nt
        tf1=tflog(kk,k);
        if(tf1<tfloglevels(1))
            tf1=.9*tfloglevels(1);%dump very small or negative values into bin 1
        end
        klevel=floor((tf1+small-tfloglevels(1))/dtflog)+1;
        if(klevel>nbins); klevel=nbins; end
        if(klevel<1); klevel=1; end
        smooth(klevel)=smooth(klevel)+tvs(kk,k);
        fold(klevel)=fold(klevel)+1;
   end
end
ind=find(fold==0);
if(~isempty(ind))
    fold(ind)=1;
end
smooth=smooth./fold;
%here we compute the standard deviation in each bin
if(nargout>2)
	smodev=zeros(1,nbins);
	%fold=smodev;
	for k=1:nf
        for kk=1:nt
            tf1=tflog(kk,k);
            %klevel=floor(tf1/deltf)+1;
            if(tf1<tfloglevels(1))
                tf1=.9*tfloglevels(1);%dump very small or negative values into bin 1
            end
            klevel=floor((tf1+small-tfloglevels(1))/dtflog)+2;
            if(klevel>nbins); klevel=nbins; end
            if(klevel<1); klevel=1; end
            smodev(klevel)=smodev(klevel)+(tvs(kk,k)-smooth(klevel))^2;
            %fold(klevel)=fold(klevel)+1;
       end
	end
	smodev=sqrt(smodev./fold);
end

%populate the output matrix
tvshyp=zeros(size(tvs));
tflogbincenters=tfloglevels-.5*dtflog;%the mean values determined previously are assigned to bincenters
for k=1:nt
    for kk=1:nf
        tf1=tflog(k,kk);
        if(tf1<tflogbincenters(1))
            tvshyp(k,kk)=smooth(1);
        elseif(tf1>tflogbincenters(end))
            tvshyp(k,kk)=smooth(end);
        else
            kbin1=floor((tf1-tflogbincenters(1))/dtflog)+1;
            kbin2=kbin1+1;
            %linear interpolation
            tvshyp(k,kk)=...
                smooth(kbin1)*(tf1-tflogbincenters(kbin2))/(tflogbincenters(kbin1)-tflogbincenters(kbin2)) + ...
                smooth(kbin2)*(tf1-tflogbincenters(kbin1))/(tflogbincenters(kbin2)-tflogbincenters(kbin1));
        end
    end
end
 tflevels=exp(tfloglevels);
