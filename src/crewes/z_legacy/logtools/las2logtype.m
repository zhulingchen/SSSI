function itype=las2logtype(mnem)

% itype=las2logtype(mnem)
%
% Convert a 4 letter las mnemonic identifying a log to a numeric logtype.
% The master list of possible numeric logtypes is to be found in 
% the source code of this function.
%
% As of Dec 7, 1995, these were:
% -1 ... unknown or un-specified
% 0  ... p-wave sonic
% 1  ... bulk density
% 2  ... formation denisty
% 3  ... apparent density
% 4  ... gamma ray
% 5  ... spontaneous potential
% 6  ... caliper
% 7  ... s-wave sonic
% 8  ... neutron porosity
% 9  ... apparent porosity
% 10 ... porosity density (LS)
% 11 ... porosity effective
% 12 ... porosity total
% 13 ... focussed resistivity
% 14 ... medium induction
% 15 ... deep induction
% 16 ... SFL resistivity
% 17 ... mel caliper
% 18 ... micronormal
% 19 ... microinverse
% 20 ... porosity density (SS)
% 21 ... Poissons ratio
% 22 ... VP/VS ratio
% 23 ... Youngs Modulus
% 24 ... Lames Lamda Constant
% 25 ... Lames Rigidity
% 26 ... Bulk Modulus
% 27 ... P-wave velocity
% 28 ... S-wave velocity
% 29 ... S-wave from array sonic 
% 30 ... P-wave from array sonic
% 31 ... Gamma ray density
% 32 ... Gamma ray porosity
% 33 ... Photoelectric cross section
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
% The term'SOFTWARE' refers to the Matlab source code, translations to
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

    if (~ischar(mnem))
        error('mnem must be a character array (string)');
    end

    % strip trailing spaces and convert to uppercase
	mnem = upper(deblank(mnem));
	
	%strip any trailing digits from the mnemonic
	p = length(mnem);
	c = mnem(p);
	while (l2lt_isdigit(c) & (p>1)) 
        p = p - 1;
        c = mnem(p);
	end
    % deblank handles the case of 'DT 5'
	mnem = deblank(mnem(1:p));
	
	switch mnem
	case {'AU','AC','DT','CO','SON','PSON','DTP','ITT','DT4P','DTCR','DTCO'...
            'DTC'}
        itype=0;
	case {'RHOB','DENE','DEN'}
        itype=1;
	case 'RHGF';
        itype=2;
	case 'RHGA'
        itype=3;
	case {'GRC','GR'}
        itype=4;
	case 'SP'
        itype=5;
	case 'CALI'
        itype=6;
	case {'DTSW','SSON','SDT','DTS','DT2','DTSD','DTSM'}
        itype=7;
	case {'NPHI','PHIN'}
        itype=8;
	case 'PHIA'
        itype=9;
	case {'PHID','DPHI'}
        itype=10;
	case {'EPHI','PHIE'}
        itype=11;
	case 'PHIT'    % as in a phit of rage
        itype=12;
	case {'SFLU','SFL'}
        itype=13;
	case 'ILM'
        itype=14;
	case 'ILD'
        itype=15;
	case 'SFLR'
        itype=16;
	case 'UNVI'
        itype=17;
	case 'MNOR'
        itype=18;
	case 'MINV'
        itype=19;
	case 'DPSS';
        itype=20;     
	case 'POIS';
        itype=21;     
	case 'VPVS';
        itype=22;     
	case 'YNGM';
        itype=23;     
	case 'LMDA';
        itype=24;     
	case 'RIDG';
        itype=25;     
	case 'BMOD';
        itype=26;     
	case 'VP';
        itype=27;     
	case 'VS';
        itype=28;
	case 'ASSW';
        itype=29;     
	case 'ASPW';
        itype=30;     
	case 'GRD';
        itype=31;     
	case 'GRP';
        itype=32;     
	case 'PEF';
        itype=33;     
	otherwise
        itype=-1;
	end
    

% Returns true if d is a digit '0'..'9', false otherwise
function x = l2lt_isdigit(d)
    x = ((d >= '0') & (d <= '9'));
