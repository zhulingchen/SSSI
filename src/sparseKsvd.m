function [A, X, err] = sparseKsvd(Y, baseSynOp, baseAnaOp, A0, blkSize, blkNum, atomSpThres, sigSpThres)
% SPARSEKSVD Sparse K-SVD dictionary training
%  SPARSEKSVD runs the sparse K-SVD dictionary training algorithm on the
%  specified set of training signals Y based on a known structured
%  transform (with synthesis operator baseSynOp and analysis operator
%  baseAnaOp), starting from an initial training dictionary A0, returning
%  the sparse dictionary representation matrix A, and the corresponding
%  signal representation matrix X and representation error err.
%
%      min  |Y-B*A*X|_F^2      s.t.  |X_i|_1 <= T
%      A,X

%% create block training data
% size: blkSize * blkSize
% amount: blkNum
dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, blkNum, 'eqdist');
blkData = sampgrid(Y, blkSize, idx{:});