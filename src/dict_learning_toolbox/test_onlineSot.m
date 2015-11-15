close all;
clear;
clc;

addpath(genpath('../'));

% load online data
load dm_list.mat
nList = length(dm_list);
% image size
nz = 120;
nx = 140;
% block size
blkSize = [10, 20];
nBlockRows = floor(nz / blkSize(1));
nBlockCols = floor(nx / blkSize(2));
spgOpts = spgSetParms('optTol', 1e-16);
% gradient classification
nClass = 5;
dDir = 180 / (2 * nClass - 1); % consider [-180:-180+dAngle/2] and [180-dAngle/2:180]
dirRange = -90 + dDir/2 : dDir : 90 - dDir/2;
dirRange = [-90, dirRange, 90];

% initial dictionary (DCT matrix) for each gradient classification
D = cell(1, nClass);
A = cell(1, nClass);
for idxClass = 1:nClass
    D{idxClass} = zeros(prod(blkSize), prod(blkSize));
    iatom = 1;
    for ic = 1:blkSize(2)
        for ir = 1:blkSize(1)
            invC = zeros(blkSize);
            invC(ir, ic) = 1;
            buffer = idct(idct(invC).').';
            D{idxClass}(:, iatom) = buffer(:);
            iatom = iatom+1;
        end
    end
    A{idxClass} = zeros(prod(blkSize), prod(blkSize));
end
nBlocks = (nz - blkSize(1) + 1) * (nx - blkSize(2) + 1);

hFigDictImg = figure(1);
hDm = figure(2);
for ii = 1:nList
    fprintf('Input image no. %d\n', ii);
    Y = dm_list{ii};
    [Y_train_patches, Y_train_patches_dirClass] = getPatches(Y, blkSize, nBlocks, dirRange);
    
    [Y_train_patches2, Y_train_patches_dirClass2] = getPatches(Y, blkSize, size(Y) ./ blkSize, dirRange);
    
    % train dictionaries by classes
    for idxClass = 1:nClass
        idx = find(Y_train_patches_dirClass == idxClass);
        [D{idxClass}, X, A{idxClass}] = onlineSotTrain(Y_train_patches(:, idx), D{idxClass}, ii, A{idxClass}, nBlocks, ...
            struct('verbosity', 0, ...
            'lambda', 0.15, ...
            'iterations', 200, ...
            'tol', 1e-3));
    end
    
%     % test Curvelet
%     tau = 1;
%     Y_curvelet_coeff = fdct_wrapping(Y, 1, 2, 4, 16);
%     [Y_curvelet_coeff_vec, sCurvelet] = curvelet2vec(Y_curvelet_coeff);
%     P = createSampler(nz * nx, nz * nx / 2, 'rdemod');
%     fdctFunc = @(x, mode) wrapper_fdct_wrapping(x, sCurvelet, 1, 4, 16, nz, nx, mode);
%     Phi = zeros(nz * nx, length(Y_curvelet_coeff_vec));
%     icol = 1;
%     for icoeff = 1:length(Y_curvelet_coeff_vec)
%         tmp_coeff = zeros(length(Y_curvelet_coeff_vec), 1);
%         tmp_coeff(icoeff) = 1;
%         tmp_img = fdctFunc(tmp_coeff, 1);
%         Phi(:, icol) = tmp_img(:);
%         icol = icol + 1;
%     end
%     Y_curvelet_coeff_vec_rec_l1 = spg_lasso(Phi, Y(:), tau, spgOpts);
    
    
%     % test L1-norm minimization
%     Y_rec_l0 = zeros(nz, nx);
%     Y_rec_l1 = zeros(nz, nx);
%     for idxBlockCol = 1:nBlockCols
%         idxCol = (idxBlockCol-1)*blkSize(2)+1:idxBlockCol*blkSize(2);
%         % process block by block in the current column
%         for idxBlockRow = 1:nBlockRows
%             idxRow = (idxBlockRow-1)*blkSize(1)+1:idxBlockRow*blkSize(1);
%             block = Y(idxRow, idxCol);
%             x_l0 = D' * block(:);   % for reference
%             P = createSampler(prod(blkSize), prod(blkSize) / 2, 'rdemod');
%             tau = norm(x_l0, 1) * 1.05;
%             x_l1 = spg_lasso(P * D, P * block(:), tau, spgOpts);
%             recBlock_l0 = D * x_l0;
%             recBlock_l0 = reshape(recBlock_l0, blkSize);
%             Y_rec_l0(idxRow, idxCol) = Y_rec_l0(idxRow, idxCol) + recBlock_l0;
%             recBlock_l1 = D * x_l1;
%             recBlock_l1 = reshape(recBlock_l1, blkSize);
%             Y_rec_l1(idxRow, idxCol) = Y_rec_l1(idxRow, idxCol) + recBlock_l1;
%         end
%     end
    
    % test SOT functions
    [Y_coeff, Y_dirClass] = forwardSot(Y, D, blkSize, dirRange);
    nRetainCoeffs = round(prod(blkSize) / 20);
    for idx_r = 1:size(Y_coeff, 1)
        for idx_c = 1:size(Y_coeff, 2)
            c = zeros(size(Y_coeff{idx_r, idx_c}));
            [~, idx] = sort(abs(Y_coeff{idx_r, idx_c}), 'descend');
            c(idx(1:nRetainCoeffs)) = Y_coeff{idx_r, idx_c}(idx(1:nRetainCoeffs));
            Y_coeff{idx_r, idx_c} = c;
        end
    end
    Y_rec = inverseSot(Y_coeff, D, blkSize, Y_dirClass);
    mse = 1/numel(Y) * norm(Y - Y_rec, 'fro')^2;
    
    % test SOT wrapper functions
    Y_coeff_vector = wrapper_sot(Y(:), D, blkSize, nz, nx, [], dirRange, 2); % SOT
    [~, Y_dirClass] = forwardSot(Y, D, blkSize, dirRange);
    Y_rec_vector = wrapper_sot(Y_coeff_vector, D, blkSize, nz, nx, Y_dirClass, dirRange, 1); % inverse SOT
    
%     Y_rec = zeros(size(Y));
%     totalBlockNum = (nz - blkSize + 1) * (nx - blkSize + 1);
%     processedBlocks = 0;
%     % take out blocks in columns
%     for ibatch = 1:blkSize:nx - blkSize + 1
%         fprintf('Batch %d... ', ibatch);
%         
%         % the current batch of blocks
%         blocks = im2colstep(Y(:, ibatch:ibatch+blkSize-1), blkSize * [1, 1], blkSize * [1, 1]);
%         
%         % remove mean values
%         [blocks, dc] = remove_dc(blocks, 'columns');
%         
%         cleanBlocks = zeros(size(blocks));
%         nRetainCoeffs = blkSize * blkSize;
%         
%         % process block by block in the current column
%         for iblk = 1:size(blocks, 2)
%             blockCoeff = D' * blocks(:, iblk);
%             recCoeff = zeros(size(blockCoeff));
%             [~, idx] = sort(abs(blockCoeff), 'descend');
%             recCoeff(idx(1:nRetainCoeffs)) = blockCoeff(idx(1:nRetainCoeffs));
%             cleanBlocks(:, iblk) = D * recCoeff;
%         end
%         
%         % add mean values
%         cleanBlocks = add_dc(cleanBlocks, dc, 'columns');
%         
%         cleanBatch = col2imstep(cleanBlocks, [nz, blkSize], blkSize * [1, 1], blkSize * [1, 1]);
%         Y_rec(:,ibatch:ibatch+blkSize-1) = Y_rec(:,ibatch:ibatch+blkSize-1) + cleanBatch;
%         
%         processedBlocks = processedBlocks + (nz - blkSize + 1);
%         fprintf('Processed %d blocks\n', processedBlocks);
%     end
%     
%     % average the overlapping area
%     cnt = countcover([nz, nx], blkSize * [1, 1], blkSize * [1, 1]);
%     Y_rec = Y_rec ./ cnt;
    
    figure(hFigDictImg);
    nSubFigRows = floor(sqrt(nClass));
    nSubFigCols = ceil(nClass / nSubFigRows);
    for idxClass = 1:nClass
        dictImg = showdict(D{idxClass}, blkSize, blkSize(1), blkSize(2), 'whitelines', 'highcontrast');
        subplot(nSubFigRows, nSubFigCols, idxClass); imshow(dictImg); title(sprintf('Dictionary Learning\n(Class %d, Iteration %d)', idxClass, ii));
    end
    figure(hDm);
    ax(1) = subplot(1, 2, 1); imagesc(Y);
    ax(2) = subplot(1, 2, 2); imagesc(Y_rec); title(sprintf('MSE = %d', mse));
    linkaxes(ax, 'xy');
end