function value=isodd(vector)
% value=isodd(vector);
% ISODD returns an array the same size as vector of logicals indicating 
% true(1) when a the number is an odd number and false(0) when the number
% is an even number.

vecdiv2=double(vector)./2;
vecround=floor(vecdiv2);
value=logical(vecdiv2-vecround);
end