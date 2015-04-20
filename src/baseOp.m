function y = baseOp(x, baseSynOp, baseAnaOp, mode)
% BASEOP synthesis / analysis operator with B that can be used as a
% function handle for spgl1 toolbox

if (mode == 1) % inverse transform or synthesis
    if (isa(baseSynOp, 'function_handle'))
        y = baseSynOp(x);
    else
        y = baseSynOp * x;
    end
elseif (mode == 2) % transform or analysis
    if (isa(baseAnaOp, 'function_handle'))
        y = baseAnaOp(x);
    else
        y = baseAnaOp * x;
    end
else
    error('Wrong mode!');
end

end