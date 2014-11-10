function value=iseven(vector)
% value=iseven(vector);
% ISEVEN returns an array the same size as vector of logicals indicating 
% true(1) when a the number is an even number and false(0) when the number
% is an odd number.

value=~isodd(vector);
end