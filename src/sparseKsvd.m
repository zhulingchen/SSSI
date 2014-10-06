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
X = zeros(coefLen, blkNum);

%% main loop for dictionary learning
for iter = 1:trainIter
    fprintf('Training Iteration %d\n', iter);
    
    % lasso for each block
    % X_i = argmin_x ||Y_i - B*A*x||_2^2 s.t. ||x||_1 <= sigSpThres
    for iblk = 1:blkNum
        opts = spgSetParms('verbosity', 1, 'optTol', 1e-20);
        X(:, iblk) = spg_lasso(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, A, mode), Y(:, iblk), sigSpThres, opts);
    end
    
    % dictionary learning and updating
    for iatom = 1:coefLen
        fprintf('Learning Atom %d\n', iatom);
        
        A(:, iatom) = zeros(coefLen, 1);
        I = (X(iatom, :) ~= 0); % I indicates the indices of the signals in Y whose representations use A(:, iatom)
        g = X(iatom, I).';
        g = g / norm(g, 2);
        YI = Y(:, I);
        XI = X(:, I);
        E = zeros(atomLen, nnz(I));
        for ii = 1:nnz(I)
            E(:, ii) = YI(:, ii) - baseSynOp(A * XI(:, ii));
        end
        z = E * g;
        
        z2 = Y(:, I) * g - learnedOp(X(:, I) * g, baseSynOp, baseAnaOp, A, 1);
        % a = argmin_a || z - B*a ||_2^2 s.t. ||a||_1 <= atomSpThres
        opts = spgSetParms('verbosity', 1, 'optTol', 1e-20);
        a = spg_lasso(@(x, mode) baseOp(x, baseSynOp, baseAnaOp, mode), z, atomSpThres, opts);
        % normalize vector a
        a = a / norm(baseSynOp(a), 2);
        A(:, iatom) = a;
        X(iatom, I) = (E' * baseSynOp(a)).';
        
        % X(iatom, I) = (Y(:, I)' * baseSynOp(a) - (baseSynOp(A * X(:, I)))' * baseSynOp(a)).';
    end
end

% calculate residue error
err = 0;
for iblk = 1:blkNum
    err = err + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1));
end

end