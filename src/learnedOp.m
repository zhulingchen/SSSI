function y = learnedOp(x, mask, baseSynOp, baseAnaOp, A, mode)
% LEARNEDOP learned synthesis / analysis operator combined with B and A
% that can be used as a function handle for spgl1 toolbox

if (mode == 1) % inverse transform or synthesis
    if (isa(baseSynOp, 'function_handle'))
        y = baseSynOp(A * x);
    else
        y = baseSynOp * (A * x);
    end
    if (isempty(mask))
        mask = ones(size(y));
    end
    y = mask .* y;
elseif (mode == 2) % transform or analysis
    if (isempty(mask))
        mask = ones(size(x));
    end
    x = mask .* x;
    if (isa(baseAnaOp, 'function_handle'))
        y = A' * baseAnaOp(x);
    else
        y = A' * (baseAnaOp * x);
    end
else
    error('Wrong mode!');
end

end