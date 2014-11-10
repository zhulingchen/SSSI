function [out, temp] = slr(in,sw,lw)
% This function calculate the STA/LTA ratio for a single trace

% check x dimension
if size(in,2) ~= 1
    in = in';
end

% deal with zero mean problem
%in = in - mean(in) + 1;
in = detrend(in,'constant')+eps;

% short and long half windows
short = floor(sw/2);
long = floor(lw/2);

% zero pad
x = padarray(in,[long, 0], 'both');

% initialization
out = zeros(1,length(in));
temp = out;

% calculate the ratio sample by sample
for i = 1 : length(in)
    a1 = i;
    a2 = i + 2 * long ;
    b1 = i + long - short;
    b2 = i + long + short;
    LTA = mean(x(a1:a2));
    STA = mean(x(b1:b2));
    out(i) = STA/(LTA+eps);
    temp(i) = LTA;
end
end
