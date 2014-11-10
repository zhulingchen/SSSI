function trout=aec(trin,sampint,op_length,trip)
% AEC: automatic envelope correction, a better AGC.
%
% Syntax
%  trout=aec(trin);
%  trout=aec(trin,sampint,op_length,trip)
%
% Description 
%   AEC performs an automatic amplitude adjustment.
%
% Method
%   1) Compute Hilbert envelope of the input trace TRIN
%   2) Convolve envelope with triangular smoother of half-length
%      OP_LENGTH
%   3) Divide input trace by smoothed envelope
%   4) Balance the output to a maximum of 1.0
%
% Inputs
%   trin      = input trace or gather of traces.
%   t or sampint   = sample interval for trin 
%               For backwards compatibility, if sampint is supplied as a time coordinate 
%               vector, the difference of the first two elements is used as
%               the sample interval.
%               Default is 0.001 seconds.  (1 millisecond)
%   op_length = half-length of triangular smoother in seconds
%               default is 1/8th of the trace length
%               ***** must be less than half the trace length *****
%   trip      = front end time before which the smoothed envelope is
%               set to a constant
%               default is op_length/10
% Outputs
%   trout     = output trace
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
   
% set defaults
if nargin < 2 || isempty(sampint)
    sampint = 0.001;   % sample interval in seconds
else
    % Backwards compatibility with time coordinate vector
    if length(sampint) > 1
        sampint = sampint(2) - sampint(1);
    end
end
if nargin < 3 || isempty(op_length)
    op_length = sampint*length(trin)/8;
end
if nargin < 4 || isempty(trip);
    trip=op_length/10.;
end

% Number of traces
ntr = size(trin,2);


% Turn scalar op_length and trip into vectors if necessary
if length(op_length) ~= ntr    
    op_length=ones(1,ntr) * op_length;
end
if length(trip) ~= ntr
    trip = ones(1,ntr) * trip;
end


% Handle 2-D or 1-D invocation differently
if numel(trin) ~= length(trin)
    % 2-D invocation
    if any((2*op_length/sampint) >= (length(trin)-2))
        error('Operator length too long, must be less than 1/2 trace length.');
    end

    for k=1:ntr
        trout(:,k) = aec_vector(trin(:,k)',sampint,op_length(k),trip(k));
    end
    
else
    % 1-D invocation
    pivotted=(ntr==1);
    if pivotted
        trin=trin.';
    end
    trout= aec_vector(trin,sampint,op_length,trip);
    if pivotted
        trout = trout';
    end
end


function trout_ = aec_vector(trin_,sampint_,op_length_,trip_)

% double the operator length
op2=op_length_*2;
% form new trace padded to a power of 2
trinew=padpow2(trin_,0);
% compute the envelope
env=abs(hilbm(trinew));
env=env(1:length(trin_));
% compute the smoothed envelope
nop=round(op2/sampint_)+1;
envsm=conv(env,triang(nop));
% grab the central length(trin) samples
envsm=envsm(round(nop/2):length(trin_)+round(nop/2)-1);
% stabilize the envelope at the ends
ntrip=round(trip_/sampint_)+1;

%for some reason, envsm(ntrip+1:length(envsm) is sometimes a column vector
%and sometimes a row vector. KWH, 2014
envsm_part = envsm(ntrip+1:length(envsm));
[m n] = size(envsm(ntrip+1:length(envsm)));
if m > 1
    envsm_part = envsm_part';
end

envsm=[envsm(ntrip)*ones(1,ntrip) envsm_part];
envsm=[envsm(1:length(envsm)-ntrip) envsm(length(envsm)-ntrip)*ones(1,ntrip)];
% correct the trace
if(sum(abs(envsm))==0)
    trout_=zeros(size(trin_));
else
    trout_=trin_./envsm;
    % balance the output to have a maximum of 1
    trout_=trout_/max(abs(trout_));
end


% balance the output to have the same mean power as input
%trout=balans(trout,trin);
