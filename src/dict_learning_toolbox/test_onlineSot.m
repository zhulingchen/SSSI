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
blkSize = 20;

D = eye(blkSize * blkSize, blkSize * blkSize);
A = zeros(blkSize * blkSize, blkSize * blkSize);
nBlocks = 2000;

hFigDictImg = figure(1);
hDm = figure(2);
for ii = 1:nList
    fprintf('Input image no. %d\n', ii);
    Y = dm_list{ii};
    [D, X, A] = onlineSotTrain(Y, D, ii, A, blkSize, nBlocks, struct('verbosity', 0, 'lambda', 0.2, 'iterations', 100, 'tol', 1e-3));
    
    % test wrapper functions
    Y_coeff_vector = wrapper_sot(Y(:), D, nz, nx, 2);
    Y_rec_vector = wrapper_sot(Y_coeff_vector, D, nz, nx, 1);
    
    % test SOT functions
    Y_coeff = forwardSot(Y, D);
    nRetainCoeffs = blkSize * blkSize;
    for idx_r = 1:size(Y_coeff, 1)
        for idx_c = 1:size(Y_coeff, 2)
            c = zeros(size(Y_coeff{idx_r, idx_c}));
            [~, idx] = sort(abs(Y_coeff{idx_r, idx_c}), 'descend');
            c(idx(1:nRetainCoeffs)) = Y_coeff{idx_r, idx_c}(idx(1:nRetainCoeffs));
            Y_coeff{idx_r, idx_c} = c;
        end
    end
    Y_rec = inverseSot(Y_coeff, D);
    
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
    
    dictImg = showdict(D, [1 1]*sqrt(size(D, 1)), round(sqrt(size(D, 2))), round(sqrt(size(D, 2))), 'whitelines', 'highcontrast');
    figure(hFigDictImg); imshow(imresize(dictImg, 2, 'nearest')); title(sprintf('Dictionary Learning Iteration %d', ii));
    figure(hDm);
    ax(1) = subplot(1, 2, 1); imagesc(Y);
    ax(2) = subplot(1, 2, 2); imagesc(Y_rec);
    linkaxes(ax, 'xy');
end