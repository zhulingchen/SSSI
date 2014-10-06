function y = learnedOp(x, baseSynOp, baseAnaOp, A, mode)
% LEARNEDOP learned synthesis / analysis operator combined with B and A
% that can be used as a function handle for spgl1 toolbox


if (mode == 1) % inverse transform or synthesis
    y = baseSynOp(A*x);
elseif (mode == 2) % transform or analysis
    y = A'*baseAnaOp(x);
else
    error('Wrong mode!');
end

end