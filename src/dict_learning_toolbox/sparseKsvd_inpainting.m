function [A, X, errLasso, errBpdn] = sparseKsvd_inpainting(Y, mask, baseSynOp, baseAnaOp, A0, trainIter, blkSize, nBlocks, atomSpThres, sigSpThres, option)
% SPARSEKSVD Sparse K-SVD dictionary learning
%  SPARSEKSVD runs the sparse K-SVD dictionary training algorithm on the
%  specified set of training signals Y based on a known structured
%  transform (with synthesis operator baseSynOp and analysis operator
%  baseAnaOp), starting from an initial dictionary A0, returning the sparse
%  dictionary representation matrix A, and the corresponding signal
%  representation matrix X and representation error err.
%
%      BPDN:    min  |X_i|_1             s.t. |Y_i-B*A*X_i|_2 <= tau    for all i
%               A,X
% or
%      Lasso:   min  |Y_i-B*A*X_i|_2     s.t. |X_i|_1 <= tau            for all i
%               A,X
%
%
% This matlab source file is free for use in academic research.
% All rights reserved.
%
% Written by Lingchen Zhu (zhulingchen@gmail.com)
% Center for Signal and Information Processing, Center for Energy & Geo Processing
% Georgia Institute of Technology

SPGOPTTOL_SIG = 1e-6;
SPGOPTTOL_ATOM = 1e-6;
N_DICTRUN_ITER = 10;
MUTCOH_THRES = 0.99;
USE_THRES = 4; % the atom must be used by this number of blocks so as to be kept
VAL_THRES = 1e-8;

dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

atomLen = blkSize * blkSize;

if (xor(isa(baseSynOp, 'function_handle'), isa(baseAnaOp, 'function_handle')))
    error('Synthesis Operator and Analysis Operator must be both function handles or both not.');
end

if (isa(baseSynOp, 'function_handle') && isa(baseAnaOp, 'function_handle'))
    coefLen = length(baseAnaOp(zeros(atomLen, 1)));
    [PhiSyn, PhiAna] = operator2matrix(baseSynOp, baseAnaOp, atomLen); % transform base transform algorithm operator into dictionary (otherwise might be too slow)
else
    coefLen = length(baseAnaOp * zeros(atomLen, 1));
    PhiSyn = baseSynOp;
    PhiAna = baseAnaOp;
end


%% create block training data
% size: blkSize * blkSize
% amount: blkNum
% Y_bak = Y;
idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, nBlocks, 'eqdist');
Y = sampgrid(Y, blkSize, idx{:});
mask = sampgrid(mask, blkSize, idx{:});
% normalization
% Y = Y - repmat(mean(Y), atomLen, 1);
nBlocks = size(Y, 2);

A = A0;
X = zeros(coefLen, nBlocks);

hFigTrainedDict = figure;
errLasso = zeros(trainIter, 1);
errBpdn = zeros(trainIter, 1);

%% main loop for dictionary learning
for iter = 1:trainIter
    fprintf('Learning Iteration %d...\t', iter);
    
    %% solve BPDN problem for each block
    % X_i = argmin_x ||Y_i - B*A*x||_2 s.t. ||x||_1 <= sigSpThres
    for iblk = 1:nBlocks
        if (option.verbosity)
            fprintf('Updating coefficients of block %d\n', iblk);
        end
        
        opts = spgSetParms('verbosity', option.verbosity, 'optTol', SPGOPTTOL_SIG);
        y = mask(:, iblk) .* Y(:, iblk);
        sigma = sqrt(nnz(mask(:, iblk)) / (blkSize * blkSize)) * sigSpThres;
        switch lower(option.method)
            case 'lasso'
                x = spg_lasso(@(x, mode) learnedOp(x, mask(:, iblk), PhiSyn, PhiAna, A, mode), y, sigma, opts);
            case 'bpdn'
                x = spg_bpdn(@(x, mode) learnedOp(x, mask(:, iblk), PhiSyn, PhiAna, A, mode), y, sigma, opts);
            otherwise
                error('Invalid optimization option! Should be either ''lasso'' or ''bpdn''');
        end
        % X(:, iblk) = OMP({@(x) PhiSyn*A*x, @(x) A'*PhiAna*x}, Y(:, iblk), sigSpThres);
        X(:, iblk) = x;
    end
    
    %% dictionary learning and updating
    unusedSig = 1:nBlocks;  % track the signals that were used to replace "dead" atoms.
    replacedAtom = zeros(1, coefLen);  % mark each atom replaced by optimize_atom
    for iatom = 1:coefLen
        if (option.verbosity)
            fprintf('Updating Atom %d\n', iatom);
        end
        
        A(:, iatom) = zeros(coefLen, 1);
        
        I = (X(iatom, :) ~= 0); % I indicates the indices of the signals in Y whose representations use B*A(:, iatom)
        % the case when no signal in Y is using B*A(:, iatom) in its representation
        if (nnz(I) <= 1)
            % err = zeros(length(unusedSig), 1);
            % for iblk = 1:length(unusedSig)
            %     err(iblk) = norm(Y(:, unusedSig(iblk)) - learnedOp(X(:, unusedSig(iblk)), PhiSyn, PhiAna, A, 1), 2)^2;
            % end
            err = sum((Y(:, unusedSig) - PhiSyn * A * X(:, unusedSig)).^2, 1);
            [~, idxErr] = max(err);
            opts = spgSetParms('verbosity', option.verbosity, 'optTol', SPGOPTTOL_ATOM);
            
            y = mask(:, unusedSig(idxErr)) .* Y(:, unusedSig(idxErr));
            a = spg_lasso(@(x, mode) baseOp(x, mask(:, iblk), PhiSyn, PhiAna, mode), y, atomSpThres, opts);
            if (norm(PhiSyn * a, 2) == 0)
                warning('a is null!');
                a(randperm(coefLen, 1)) = atomSpThres;
            end
            a = a / norm(PhiSyn * a, 2);
            A(:, iatom) = a;
            unusedSig = unusedSig([1:idxErr-1, idxErr+1:end]);
            replacedAtom(iatom) = 1;
            continue;
        end
        
        g = X(iatom, I).';
        a = A(:, iatom);
        if (isnan(g))
            error ('g is NaN!');
        end
        
        % YI = Y(:, I);
        % XI = X(:, I);
        % E = zeros(atomLen, nnz(I));
        % for ii = 1:nnz(I)
        %     E(:, ii) = YI(:, ii) - learnedOp(XI(:, ii), PhiSyn, PhiAna, A, 1);
        % end
        % z = E * g;
        E = Y(:, I) - PhiSyn * A * X(:, I);
        for iter_dictrun = 1:N_DICTRUN_ITER
            % normalize g such that ||g||_2 = 1
            g = g / norm(g, 2);
            
            z = (mask(:, I) .* E + (ones(size(mask(:, I))) - mask(:, I)) .* (PhiSyn * a * g')) * g;
            
            % a = argmin_a || z - B*a ||_2 s.t. ||a||_1 <= atomSpThres
            opts = spgSetParms('verbosity', option.verbosity, 'optTol', SPGOPTTOL_ATOM);
            a = spg_lasso(@(x, mode) baseOp(x, [], PhiSyn, PhiAna, mode), z, atomSpThres, opts);
            % a = OMP({@(x) (PhiSyn*x), @(x) (PhiAna*x)}, z, atomSpThres);
            % normalize a such that ||PhiSyn * a||_2 = 1
            if (norm(PhiSyn * a, 2) == 0)
                warning('a is null!');
                a(randperm(coefLen, 1)) = atomSpThres;
            end
            a = a / norm(PhiSyn * a, 2);
            
            g = (mask(:, I) .* E + (ones(size(mask(:, I))) - mask(:, I)) .* (PhiSyn * a * g'))' * (PhiSyn * a);
        end
        A(:, iatom) = a;
        X(iatom, I) = g.';
    end
    
    %% dictionary clearing
    % err = zeros(1, blkNum);
    % for iblk = 1:blkNum
    %     err(iblk) = norm(Y(:, iblk) - learnedOp(X(:, iblk), PhiSyn, PhiAna, A, 1), 2)^2;
    % end
    err = sum((Y - PhiSyn * A * X).^2, 1);
    
    numClearedAtom = 0;
    useCount = sum(abs(X)>VAL_THRES, 2);
    
    for iatom = 1:coefLen
        % compute mutual coherence
        mutCoh = learnedOp(mask(:, iatom) .* (PhiSyn * A(:, iatom)), mask(:, iatom), PhiSyn, PhiAna, A, 2);
        mutCoh(iatom) = 0; % excluding self coherence (=1)
        
        % replace atoms if they do not meet requirements
        if ( (max(abs(mutCoh))>MUTCOH_THRES || useCount(iatom) < USE_THRES) && ~replacedAtom(iatom) )
            [~, idxErr] = max(err(unusedSig));
            y = mask(:, unusedSig(idxErr)) .* Y(:, unusedSig(idxErr));
            opts = spgSetParms('verbosity', option.verbosity, 'optTol', SPGOPTTOL_ATOM);
            a = spg_lasso(@(x, mode) baseOp(x, mask(:, iblk), PhiSyn, PhiAna, mode), y, atomSpThres, opts);
            if (norm(PhiSyn * a, 2) == 0)
                warning('a is null!');
                a(randperm(coefLen, 1)) = atomSpThres;
            end
            a = a / norm(PhiSyn * a, 2);
            A(:, iatom) = a;
            unusedSig = unusedSig([1:idxErr-1, idxErr+1:end]);
            numClearedAtom = numClearedAtom + 1;
        end
    end
    
    %% calculate residue error
    errLasso(iter) = norm(Y - PhiSyn * A * X, 'fro');
    errBpdn(iter) = sum(abs(X(:)));
    fprintf('Lasso error |Y - B*A*X|_2 = %f, BPDN error |X|_1 = %f\n', errLasso(iter), errBpdn(iter));
    
    %% show trained dictionary
    dictImg = showdict(PhiSyn * A, [1 1]*sqrt(size(PhiSyn * A, 1)), round(sqrt(size(PhiSyn * A, 2))), round(sqrt(size(PhiSyn * A, 2))), 'whitelines', 'highcontrast');
    figure(hFigTrainedDict); imshow(imresize(dictImg, 2, 'nearest')); title(sprintf('Trained Dictionary (Iteration %d)', iter));
    
end

% % calculate residue error
% % err = 0;
% % for iblk = 1:blkNum
% %     err = err + norm(Y(:, iblk) - learnedOp(X(:, iblk), PhiSyn, PhiAna, A, 1), 2);
% % end
% err = norm(Y - PhiSyn * A * X, 'fro');

end