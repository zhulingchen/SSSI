function sout=predlowf(s,f,flower,fupper,n)
% sout=predlowf(s,f,flower,fupper,n)
%
% PREDLOWF uses Burg prediction filtering to predict the low portion of
% a frequency spectrum.
% Algorithm:
%	- Copy the input spectrum into the output spectrum and zero
%	all values below flower
%	- Loop over zero'd frequencies starting with the highest
%	- Design a 1 lag Burg prediction filter to predict the
%	first zero'd frequency. Use only frequencies between the
%	first non-zero'd one and fupper in the design of the filter.
%	- predict the highest zero'd spectral component.
%	- use the predictions of previous iterations in the current
%	prediction
%
% s ... input complex spectrum
% f ... frequency coordinate vector for s
% flower ... lowest frequency to keep
% fupper ... higest frequency to use in the prediction filter design
% n ... number of points in the prediction operator
% sout ... returned spectrum with low-end replaced with predicted
%	samples.
%
% G.F. Margrave May 1995
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
% check for row vectors and transpose if needed
aa=size(s);
bb=size(f);
if aa(1)==1
	s=s';
end
if bb(1)==1
	f=f';
end
ind=find(f<flower);
ind=flipud(ind); %start with the higher f's
iup=find(f<=fupper);
iup=max(iup);
sout=s;
%zero the frequencies to be predicted
sout(ind)=zeros(size(ind));
%handle real and imaginary part separately
tempr=real(sout);
tempi=imag(sout);
for k=1:length(ind)
	%design filter
	pfilt=(burgpr(tempr(ind(k):iup).',n)).';
	%convert from prediction error to prediction filter
	pfilt(1)=pfilt(1)+1;
	
	%predict tempr(ind(k))
	p= sum(tempr(ind(k):ind(k)+n-1).*pfilt);
	tempr(ind(k))=p;
end
for k=1:length(ind)
	%design filter
	pfilt=(burgpr(tempi(ind(k):iup).',n)).';
	%convert from prediction error to prediction filter
	pfilt(1)=pfilt(1)+1;
	
	%predict tempi(ind(k))
	p= sum(tempi(ind(k):ind(k)+n-1).*pfilt);
	tempi(ind(k))=p;
end
sout=tempr+i*tempi;
