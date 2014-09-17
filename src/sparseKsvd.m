function [A, X, err] = sparseKsvd(Y, baseSynOp, baseAnaOp, A0, trainIter, blkSize, blkNum, atomSpThres, sigSpThres)
% SPARSEKSVD Sparse K-SVD dictionary training
%  SPARSEKSVD runs the sparse K-SVD dictionary training algorithm on the
%  specified set of training signals Y based on a known structured
%  transform (with synthesis operator baseSynOp and analysis operator
%  baseAnaOp), starting from an initial dictionary A0, returning the sparse
%  dictionary representation matrix A, and the corresponding signal
%  representation matrix X and representation error err.
%
%      min  |Y-B*A*X|_F^2      s.t.  |X_i|_1 <= T
%      A,X

dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

atomLen = blkSize * blkSize;
coefLen = length(baseAnaOp(zeros(atomLen, 1)));

%% create block training data
% size: blkSize * blkSize
% amount: blkNum
Y_bak = Y;
idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, blkNum, 'eqdist');
Y = sampgrid(Y, blkSize, idx{:});

A = A0;
X = zeros(coefLen, blkSize);

%% main loop for dictionary learning
for iter = 1:trainIter
    % lasso for each block
    % X_i = argmin_x ||Y_i - B*A*x||_2^2 s.t. ||x||_1 < Thres
    for iblk = 1:blkNum
        opts = spgSetParms('verbosity', 1, 'optTol', 1e-12);
        X(:, iblk) = spg_lasso(@(x, mode) combOp(x, baseSynOp, baseAnaOp, A, mode), Y(:, iblk), sigSpThres, opts);
    end
    
    % dictionary learning and updating
    for icoef = 1:coefLen
        
    end
end


end


%% combined synthesis / analysis operator with B and A that can be used as a function handle for spgl1 toolbox
function y = combOp(x, baseSynOp, baseAnaOp, A, mode)

if (mode == 1) % inverse transform or synthesis
    y = baseSynOp(A*x);
elseif (mode == 2) % transform or analysis
    y = A'*baseAnaOp(x);
else
    error('Wrong mode!');
end

end