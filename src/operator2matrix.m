function [PhiSyn, PhiAna] = operator2matrix(synOp, anaOp, atomLen)
% OPERATOR2MATRIX converts synthesis and analysis operators into matrices


coefLen = length(anaOp(zeros(atomLen, 1)));

PhiSyn = zeros(atomLen, coefLen);
for ic = 1:coefLen
    c = zeros(coefLen, 1);
    c(ic) = 1;
    a = synOp(c);
    PhiSyn(:, ic) = a(:);
end

PhiAna = zeros(coefLen, atomLen);
for ia = 1:atomLen
    a = zeros(atomLen, 1);
    a(ia) = 1;
    c = anaOp(a);
    PhiAna(:, ia) = c(:);
end