function y = learnedOp(x, mask, baseSynOp, baseAnaOp, A, mode)
% LEARNEDOP learned synthesis / analysis operator combined with B and A
% that can be used as a function handle for spgl1 toolbox


if (mode == 1) % inverse transform or synthesis
    if (isa(baseSynOp, 'function_handle'))
        y = mask * (baseSynOp(A * x));
    else
        y = mask * (baseSynOp * (A * x));
    end
elseif (mode == 2) % transform or analysis
    if (isa(baseAnaOp, 'function_handle'))
        y = A' * baseAnaOp(pinv(mask) * x);
    else
        y = A' * (baseAnaOp * (pinv(mask) * x));
    end
else
    error('Wrong mode!');
end

end