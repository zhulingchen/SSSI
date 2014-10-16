function [A, X, err] = sparseKsvd(Y, baseSynOp, baseAnaOp, A0, trainIter, blkSize, blkNum, atomSpThres, sigSpThres, option)
% SPARSEKSVD Sparse K-SVD dictionary training
%  SPARSEKSVD runs the sparse K-SVD dictionary training algorithm on the
%  specified set of training signals Y based on a known structured
%  transform (with synthesis operator baseSynOp and analysis operator
%  baseAnaOp), starting from an initial dictionary A0, returning the sparse
%  dictionary representation matrix A, and the corresponding signal
%  representation matrix X and representation error err.
%
%      BPDN:    min  |X_i|_1             s.t.  |Y_i-B*A*X_i|_2 <= tau    for all i
%               A,X
% or
%      Lasso:   min  |Y_i-B*A*X_i|_2     s.t. |X_i|_1 <= tau             for all i
%               A,X

SPGOPTTOL_SIG = 1e-6;
SPGOPTTOL_ATOM = 1e-6;
MUTCOH_THRES = 0.99;
USE_THRES = 3; % the atom must be used by this number of blocks so as to be kept
VAL_THRES = 1e-8;

dim = ndims(Y);
if ( dim < 2 || dim > 3 )
    error('Only 2-D and 3-D signals are supported!');
end

atomLen = blkSize * blkSize;
coefLen = length(baseAnaOp(zeros(atomLen, 1)));
[PhiSyn, PhiAna] = operator2matrix(baseSynOp, baseAnaOp, atomLen); % transform base transform algorithm operator into dictionary (otherwise might be too slow)

%% create block training data
% size: blkSize * blkSize
% amount: blkNum
% Y_bak = Y;
idx = cell(dim, 1);
[idx{:}] = reggrid(size(Y)-blkSize+1, blkNum, 'eqdist');
Y = sampgrid(Y, blkSize, idx{:});
blkNum = size(Y, 2);

A = A0;
X = zeros(coefLen, blkNum);

hFigTrainedDict = figure;
errBpdn = zeros(trainIter, 1);
errLasso = zeros(trainIter, 1);

%% main loop for dictionary learning
for iter = 1:trainIter
    fprintf('Learning Iteration %d\n', iter);
    
    %% solve BPDN problem for each block
    % X_i = argmin_x ||Y_i - B*A*x||_2^2 s.t. ||x||_1 <= sigSpThres
    for iblk = 1:blkNum
        opts = spgSetParms('verbosity', 1, 'optTol', SPGOPTTOL_SIG);
        switch lower(option)
            case 'lasso'
                X(:, iblk) = spg_lasso(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, A, mode), Y(:, iblk), sigSpThres, opts);
            case 'bpdn'
                X(:, iblk) = spg_bpdn(@(x, mode) learnedOp(x, baseSynOp, baseAnaOp, A, mode), Y(:, iblk), sigSpThres, opts);
            otherwise
                error('Invalid optimization option! Should be either ''lasso'' or ''bpdn''');
        end
        % X(:, iblk) = OMP({@(x) baseSynOp(A*x), @(x) A'*baseAnaOp(x)}, Y(:, iblk), sigSpThres);
    end
    
    % calculate residue error
    errBpdn(iter) = sum(abs(X(:)));
    errLasso(iter) = norm(Y - PhiSyn * A * X, 'fro');
    
    %% dictionary learning and updating
	unusedSig = 1:blkNum;  % track the signals that were used to replace "dead" atoms.
    replacedAtom = zeros(1, coefLen);  % mark each atom replaced by optimize_atom
    for iatom = 1:coefLen
        fprintf('Learning Atom %d\n', iatom);
        
        A(:, iatom) = zeros(coefLen, 1);
        
        I = (X(iatom, :) ~= 0); % I indicates the indices of the signals in Y whose representations use B*A(:, iatom)
        % the case when no signal in Y is using B*A(:, iatom) in its representation
        if (nnz(I) == 0)
            % err = zeros(length(unusedSig), 1);
            % for iblk = 1:length(unusedSig)
            %     err(iblk) = norm(Y(:, unusedSig(iblk)) - learnedOp(X(:, unusedSig(iblk)), baseSynOp, baseAnaOp, A, 1), 2)^2;
            % end
            err = sum((Y(:, unusedSig) - PhiSyn * A * X(:, unusedSig)).^2, 1);
            [~, idxErr] = max(err);
            opts = spgSetParms('verbosity', 1, 'optTol', SPGOPTTOL_ATOM);
            a = spg_lasso(@(x, mode) baseOp(x, baseSynOp, baseAnaOp, mode), Y(:, unusedSig(idxErr)), atomSpThres, opts);
            a = a / norm(baseSynOp(a), 2);
            if (isnan(a))
                error ('a is NaN!');
            end
            A(:, iatom) = a;
            unusedSig = unusedSig([1:idxErr-1, idxErr+1:end]);
            replacedAtom(iatom) = 1;
            continue;
        end
        g = X(iatom, I).';
        g = g / norm(g, 2);
        if (isnan(g))
            error ('g is NaN!');
        end
        % YI = Y(:, I);
        % XI = X(:, I);
        % E = zeros(atomLen, nnz(I));
        % for ii = 1:nnz(I)
        %     E(:, ii) = YI(:, ii) - learnedOp(XI(:, ii), baseSynOp, baseAnaOp, A, 1);
        % end
        % z = E * g;
        z = Y(:, I) * g - PhiSyn * A * X(:, I) * g;
        
        % a = argmin_a || z - B*a ||_2^2 s.t. ||a||_1 <= atomSpThres
        opts = spgSetParms('verbosity', 1, 'optTol', SPGOPTTOL_ATOM);
        a = spg_lasso(@(x, mode) baseOp(x, baseSynOp, baseAnaOp, mode), z, atomSpThres, opts);
        % a = OMP({baseSynOp, baseAnaOp}, z, atomSpThres);
        % normalize vector a
        a = a / norm(baseSynOp(a), 2);
        if (isnan(a))
            error ('a is NaN!');
        end
        
        A(:, iatom) = a;
        % X(iatom, I) = (E' * baseSynOp(a)).';
        
        X(iatom, I) = (Y(:, I)' * PhiSyn * a - (PhiSyn * A * X(:, I))' * PhiSyn * a).';
    end
    
    %% dictionary clearing
    % err = zeros(1, blkNum);
    % for iblk = 1:blkNum
    %     err(iblk) = norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2)^2;
    % end
    err = sum((Y - PhiSyn * A * X).^2, 1);
    
    numClearedAtom = 0;
    useCount = sum(abs(X)>VAL_THRES, 2);
    
    for iatom = 1:coefLen
        % compute mutual coherence
        mutCoh = learnedOp(baseSynOp(A(:, iatom)), baseSynOp, baseAnaOp, A, 2);
        mutCoh(iatom) = 0; % excluding self coherence (=1)
        
        % replace atoms if they do not meet requirements
        if ( (max(abs(mutCoh))>MUTCOH_THRES || useCount(iatom) < USE_THRES) && ~replacedAtom(iatom) )
            [~, idxErr] = max(err(unusedSig));
            opts = spgSetParms('verbosity', 1, 'optTol', SPGOPTTOL_ATOM);
            a = spg_lasso(@(x, mode) baseOp(x, baseSynOp, baseAnaOp, mode), Y(:, unusedSig(idxErr)), atomSpThres, opts);
            a = a / norm(baseSynOp(a), 2);
            if (isnan(a))
                error ('a is NaN!');
            end
            A(:, iatom) = a;
            unusedSig = unusedSig([1:idxErr-1, idxErr+1:end]);
            numClearedAtom = numClearedAtom + 1;
        end
    end
    
    %% show trained dictionary
    dictimg = showdict(PhiSyn * A, [1 1]*sqrt(size(PhiSyn * A, 1)), round(sqrt(size(PhiSyn * A, 2))), round(sqrt(size(PhiSyn * A, 2))), 'whitelines', 'highcontrast');
    figure(hFigTrainedDict); imshow(imresize(dictimg,2,'nearest')); title('Trained Dictionary');
    
end

% calculate residue error
% err = 0;
% for iblk = 1:blkNum
%     err = err + norm(Y(:, iblk) - learnedOp(X(:, iblk), baseSynOp, baseAnaOp, A, 1), 2)^2;
% end
err = norm(Y - PhiSyn * A * X, 'fro');

end