function y = baseOp(x, baseSynOp, baseAnaOp, mode)
% BASEOP synthesis / analysis operator with B that can be used as a
% function handle for spgl1 toolbox

if (mode == 1) % inverse transform or synthesis
    y = baseSynOp(x);
elseif (mode == 2) % transform or analysis
    y = baseAnaOp(x);
else
    error('Wrong mode!');
end

end