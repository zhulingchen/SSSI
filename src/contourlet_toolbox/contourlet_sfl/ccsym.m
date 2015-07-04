%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Yue M. Lu and Minh N. Do
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%	ccsym.m
%	
%	First created: 08-14-05
%	Last modified: 07-13-09
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = ccsym(x, k, type)

% Exploit the complex conjugate symmetry in the fourier transform of
% real-valued signals
%
% INPUT
%
%   x: the input signal
%   
%   type: 'c' compaction 'e' expansion
%
%   k: along which dimension
%
% OUTPUT
%
%   y: the compacted (or expanded) signal in the frequency domain

% Dimensions of the problem
N = ndims(x);
szX = size(x);

if type == 'c'
    % initialize the subscript array
    sub_array = repmat({':'}, [N, 1]);  
    sub_array{k} = 1 : szX(k) / 2 + 1;
    y = x(sub_array{:});
else
    % subscript mapping for complex conjugate symmetric signal recovery
    szX(k) = (szX(k)-1) * 2;
    sub_conj = cell(N, 1);

    for m = 1 : N
        sub_conj{m} = [1 szX(m):-1:2];
    end
    
    sub_conj{k} = [szX(k)/2 : -1 : 2];
    % recover the full signal
    y = cat(k, x, conj(x(sub_conj{:})));

end
