function [f0, f1] = ffilters(h0, h1);
% FFILTERS	Fan filters from diamond shape filters
%
%	[f0, f1] = ffilters(h0, h1);
        
f0 = cell(1, 4);
f1 = cell(1, 4);

% For the first half channels
f0{1} = modulate2(h0, 'r');
f1{1} = modulate2(h1, 'r');

f0{2} = modulate2(h0, 'c');
f1{2} = modulate2(h1, 'c');

% For the second half channels, 
% use the transposed filters of the first half channels
f0{3} = f0{1}';
f1{3} = f1{1}';

f0{4} = f0{2}';
f1{4} = f1{2}';
