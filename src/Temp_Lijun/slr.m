function out = slr(in,sw,lw)
% This function calculate the STA/LTA ratio for a single trace

% check x dimension
if size(in,1) ~= 1
    in = in';
end

% zero pad
x = padarray(in,[0, ceil(lw/2)], 'both');
for 
end
