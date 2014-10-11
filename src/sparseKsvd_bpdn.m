function [A, X, err] = sparseKsvd_bpdn(Y, baseSynOp, baseAnaOp, A0, trainIter, blkSize, blkNum, atomSpThres, sigSpThres)
% SPARSEKSVD Sparse K-SVD dictionary training
%  SPARSEKSVD runs the sparse K-SVD dictionary training algorithm on the
%  specified set of training signals Y based on a known structured
%  transform (with synthesis operator baseSynOp and analysis operator
%  baseAnaOp), starting from an initial dictionary A0, returning the sparse
%  dictionary representation matrix A, and the corresponding signal
%  representation matrix X and representation error err.
%
%      min  |X_i|_1 <= T      s.t.  |Y_i-B*A*X_i|_2 <= tau    for all i
%       X

dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

atomLen = blkSize * blkSize;
coefLen = length(baseAnaOp(zeros(atomLen, 1)));
[PhiSyn, ~] = operator2matrix(baseSynOp, baseAnaOp, atomLen);

%% create block training data
% size: blkSize * blkSize
% amount: blkNum
Y_bak = Y;
idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, blkNum, 'eqdist');
Y = sampgrid(Y, blkSize, idx{:});
blkNum = size(Y, 2);

A = A0;
X = zeros(coefLen, blkNum);

spgOptTol_sig = 1e-6;
spgOptTol_atom = 1e-6;
% spgOptTol = 1e-3;

hFigTrainedDict = figure;
errBpdn = zeros(trainIter, 1);
errLasso = zeros(trainIter, 1);
errBefore = zeros(coefLen, trainIter);
errAfter = zeros(coefLen, trainIter);

%% main loop for dictionary learning
for iter = 1:trainIter
    fprintf('Learning Iteration %d\n', iter);
    
    % lasso for each block
    % X_i = argmin_x ||Y_i - B*A*x||_2^2 s.t. ||x||_1 <= sigSpThres
    for iblk = 1:blkNum
        opts = spgSetParms('verbosity', 1, 'optTol', spgOptTol_sig);
        X(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, A, mode), Y(:, iblk), sigSpThres, opts);
    end
    
    % calculate residue error
    errBpdn(iter) = 0;
    errLasso(iter) = 0;
    for iblk = 1:blkNum
        errBpdn(iter) = errBpdn(iter) + norm(X(:, iblk), 1);
        errLasso(iter) = errLasso(iter) + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2);
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
            E(:, ii) = YI(:, ii) - learnedOp(XI(:, ii), baseSynOp, baseAnaOp, A, 1);
        end
        z = E * g;
        % z = Y(:, I) * g - learnedOp(X(:, I) * g, baseSynOp, baseAnaOp, A, 1);
        
        % a = argmin_a || z - B*a ||_2^2 s.t. ||a||_1 <= atomSpThres
        opts = spgSetParms('verbosity', 0, 'optTol', spgOptTol_atom);
        a = spg_lasso(@(x, mode) baseOp(x, baseSynOp, baseAnaOp, mode), z, atomSpThres, opts);
        % baseCell = {baseSynOp, baseAnaOp};
        % a2 = OMP(baseCell, z, atomSpThres);
        % normalize vector a
        a = a / norm(baseSynOp(a), 2);
        
        % calculate residue error
        errBefore(iatom, iter) = 0;
        for iblk = 1:blkNum
            errBefore(iatom, iter) = errBefore(iatom, iter) + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2);
        end
        fprintf('err before update = %.6e\n', errBefore(iatom, iter));
        
        A(:, iatom) = a;
        X(iatom, I) = (E' * baseSynOp(a)).';
        
        % calculate residue error
        errAfter(iatom, iter) = 0;
        for iblk = 1:blkNum
            errAfter(iatom, iter) = errAfter(iatom, iter) + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2);
        end
        fprintf('err after update = %.6e\n', errAfter(iatom, iter));
        
        % X(iatom, I) = (Y(:, I)' * baseSynOp(a) - X(:, I)' * learnedOp(baseSynOp(a), baseSynOp, baseAnaOp, A, 2)).';
    end
    
    % show trained dictionary
    dictimg = showdict(PhiSyn * A, [1 1]*sqrt(size(PhiSyn * A, 1)), round(sqrt(size(PhiSyn * A, 2))), round(sqrt(size(PhiSyn * A, 2))), 'whitelines', 'highcontrast');
    figure(hFigTrainedDict); imshow(imresize(dictimg,2,'nearest')); title('Trained Dictionary');
    
end

% calculate residue error
err = 0;
for iblk = 1:blkNum
    err = err + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2);
end

end