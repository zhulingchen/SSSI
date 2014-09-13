function coeffs = decdemo( im, option )
% DECDEMO   Demonstrates contourlet decomposition and reconstruction.
%
%   DECDEMO shows how to use the contourlet toolbox to decompose
%   and reconstruct an image.  It provides a sample script that uses
%   basic functions such as pdfbdec, pdfbrec, showpdfb.
%
%   It can be modified for applications such as image analysis,
%   image retrieval and image processing.
%
%   While displaying images, the program will pause and wait for your response.
%   When you are ready, you can just press Enter key to continue.
%
%       decdemo( [im, option] )
%
% Input:
%	image:  a double or integer matrix for the input image.
%           The default is the zoneplate image.
%   option: option for the demos. The default value is 'auto'
%       'auto' ------  automtatical demo, no input
%       'user' ------  semi-automatic demo, simple interactive inputs
%       'expert' ----  mannual, complete interactive inputs.
%                      (Not implmented in this version)
%
% Output:
%	coeffs: a cell vector for the contourlet decomposition coefficients.
%
% See also:     PDFBDEC, PDFBREC, SHOWPDFB

% History:
%   10/20/2003  Creation.
%   10/22/2003  Change the user interface for better image display.

close all;

disp('Welcome to the contourlet decomposition demo! :)');
disp('Type help decdemo for help' ) ;
disp('You can also view decdemo.m for details.') ;
disp(' ');

% Input image
if ~exist('im', 'var')
    % Zoneplate image: good for illustrating multiscale and directional
    % decomposition
    im = imread ('zoneplate.png') ;
    % load('im.mat');
end

% Show the input image
[nx, ny] = size(im);
disp( 'Displaying the input image...');
figure; imagesc(im, [0, 255]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
% input( 'Press Enter key to continue...' ) ;
disp( ' ' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image decomposition by contourlets using the
% pyramidal directional filter bank (PDFB).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameteters:
nlevels = [2, 3] ;      % Decomposition level
pfilter = 'pkva' ;      % Pyramidal filter
dfilter = 'pkva' ;      % Directional filter

% Contourlet decomposition
coeffs = pdfbdec( double(im), pfilter, dfilter, nlevels );

% Display the coefficients
disp('Displaying the contourlet coefficients...') ;
figure; imcoeff = showpdfb( coeffs ) ;
title('Contourlet coefficients');
% input('Press Enter key to continue...' ) ;
disp(' ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the dictionary of contourlet transform
% Lingchen Zhu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v, s] = pdfb2vec(coeffs);
nVec = length(v);
Phi = zeros(nx * ny, nVec);
for iv = 1:nVec
    display(iv);
    vtmp = zeros(nVec, 1);
    vtmp(iv) = 1;
    imtmp = pdfbrec( vec2pdfb(vtmp, s), pfilter, dfilter ) ;
    Phi(:, iv) = imtmp(:);
end
im2 = Phi * v;
im2 = reshape(im2, nx, ny);
delta = double(im) - im2;
display(max(abs(delta(:))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyramidal directional filter bank (PDFB) reconstruction.
% This is the inverse of pdfbdec, i.e.
% imrec = pdfbrec(coeffs, pfilter, dfilter);
% would reconstruct imrec = im
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reconstruct image
imrec = pdfbrec( coeffs, pfilter, dfilter ) ;

disp('Displaying the reconstructed image...') ;
disp('It should be a perfect reconstruction' ) ;
disp(' ') ;

% Show the reconstruction image and the original image
figure;
subplot(1,2,1), imagesc( im, [0, 255] );
title('Original image' ) ;
axis image off;
subplot(1,2,2), imagesc( imrec, [0, 255] );
title('Reconstructed image' ) ;
axis image off;

mse = sum( sum( (imrec - double(im)).^2 ) );
mse = mse / prod(size(im));

disp( sprintf('The mean square error is: %f', mse ) );
disp(' ');