function ieee = ibm2ieee (ibm, fmt)
%
% function ieee=ibm2ieee (ibm, fmt)
% Where:
%   ibm  = single precision ibm float(s) as uint32
%   ieee = single precision ieee float(s) as fmt
%    fmt = output format: uint32 (suitable for fwrite), single or 
%          double precision
%
% Example (Tested on little-endian 64-bit Windows PC, Matlab R2007b:
%   >> fid = fopen('ibmtest.sgy','r','ieee-be');
%   >> fseek(fid,3200+400+240);
%   >> ibm = fread(fid,4,'uint32=>uint32');
%   >> ibm2ieee(ibm)
%
% References
%   http://bytes.com/topic/c/answers/...
%        221981-c-code-converting-ibm-370-floating-point-ieee-754-a
%        => C code contributed by lawrence.jones@ugs.com
%
% Modified and tested in Matlab 2009b on Windows Vista 64-bit PC
%   Kevin Hall, CREWES, University of Calgary
%   February 2010
%

% Check input type
if(isa(ibm,'uint32'))
    val = ibm;
else
    error('Input must be of type ''uint32''');
end

% Set output format if not specified
if(nargin == 1)
    fmt = 'double';
end

% Get signbit, exponent and fraction
signbit  = bitshift(val,-31);     % get sign bit
fraction = bitshift(val,1);       % shift sign bit out
exponent = bitshift(fraction,-25);% store exponent
fraction = bitshift(fraction,7);  % shift exponent out

% Adjust exponents from base 16 offset 64 to base 2 offset 127
% (exponent -64)*4 +127 -1 == exponent*4 -130 == bitshift(exponent,2) -130
%  - Needs to be int32 to allow for negative exponents in the (re)normalize
%  step
exponent = int32(bitshift(exponent,2) -130);
  
for i=1:numel(fraction) 
    % For every number in the array; required because of Matlab array
    % behaviour in the (re)normalize while loop
    if (fraction(i) == 0)
        % Found zero
        exponent(i)=0;
    else
        % (Re)normalize
        while (fraction(i) < uint32(hex2dec('80000000')))
           exponent(i) = exponent(i) -1;
           fraction(i) = bitshift(fraction(i),1);
        end %end while    

        % Sanity checks
        if (exponent(i) <= 0) % Undeflow        
            if (exponent(i) < -24) % Complete underflow
                warning('IBM:underflow:completeunderflow',...
                    'Complete underflow - returning zero')
                fraction(i) = 0;
            else % Partial underflow
                warning('IBM:underflow:partialunderflow',...
                    'Partial underflow - returning denormalized number')
                shift = -1 *exponent(i);
                fraction(i) = bitshift(fraction(i),-1*shift);
            end
            exponent(i) = 0;
        elseif (exponent(i) >= 255) % Overflow
            warning('IBM:overflow',...
                'Overflow - returning Inf')
            fraction(i) = 0;
            exponent(i) = 255;
        else
            % Regular number: remove implied IEEE digit (bit 24)
            fraction(i) = bitshift(fraction(i),1);
        end %end sanity checks
    end  
end %end for

% Assemble the answer
signbit  = bitshift(signbit,31);
exponent = bitshift(uint32(exponent),23);
fraction = bitshift(fraction,-9);
ieee = bitor(signbit,bitor(exponent,fraction));

% Format for output
switch(fmt)
    case('uint32')
        %return uint32
    case('single')
        ieee = typecast(ieee,'single');
    case('double')
        ieee = double(typecast(ieee,'single'));
    otherwise
        %return uint32
end





