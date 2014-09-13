function nlademo( im, option )
% NLADEMO   Demo for contourlet nonlinear approximation. 
%   NLADEMO shows how to use the contourlet toolbox to do nonlinear 
%   approximation. It provides a sample script that uses basic functions 
%   such as pdfbdec, pdfbrec, showpdfb, pdfb_tr, pdfb2vec and vec2pdfb.
%
%   It can be modified for applications such as denoising, compression, 
%   and computer vision.
%
%   While displaying images, the program will pause and wait for your response.
%   When you are ready, you can just press Enter key to continue.
%
%       nlademo( [im, option] )
%
% Input:
%	im:     a double or integer matrix for the input image.
%           The default input is the 'peppers' image.  
%   option: option for the demos. The default value is 'auto'
%       'auto' ------  automtatical demo, no input
%       'user' ------  semi-automatic demo, simple interactive inputs
%       'expert' ----  mannual, complete interactive inputs. 
%                      (It is same as 'user' in this version)
%   
% See also:     PDFBDEC, PDFBREC, SHOWPDFB, PDFB2VEC, VEC2PDFB

% History:
%   10/21/2003  Creation.
%   10/22/2003  Change the user interface for better image display. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to contourlet nonlinear approximation demo! :)');
disp('Type help nlademo for help' ) ;
disp('You can also view nlademo.m for details.') ;
disp(' ');

% Input image
if ~exist('im', 'var')
    im = imread('barbara.png');
    im = double(im) / 256;
end

disp( 'Displaying the input image...');
clf;
imagesc(im, [0, 1]);
title( 'Input image' ) ;
axis image off;
colormap(gray);
input( 'Press Enter key to continue...' ) ;
disp(' ');

% Running option
if ~exist('option', 'var')
    option = 'auto' ;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image decomposition by contourlets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameteters:
nlevels = [0, 0, 0, 4, 5];   % Decomposition level
pfilter = '9-7' ;            % Pyramidal filter
dfilter = 'pkva';            % Directional filter

% Contourlet transform
coeffs = pdfbdec( im, pfilter, dfilter, nlevels );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nonlinear approximation.
% It will keep the most significant coefficients and use these 
% coefficients to reconstruct the image.
% It will show the reconstructed image and calculate the distortion. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display information
nPixels = prod( size(im) );             % number of pixels
nCoeffs = length(pdfb2vec(coeffs));     % number of PDFB coefficients

disp( sprintf('Number of image pixels is %d', nPixels) ) ;
disp( sprintf('Number of coefficients is %d', nCoeffs) ) ;
if strcmp( option, 'auto' )
    nSignif = round(nPixels * 2.5 / 100) ;  % 2.5% of coefficients
    disp( sprintf( 'It will keep %d significant coefficients', nSignif ) ) ;
else
    % Get the input and check the input
    nSignif = -1 ;
    while nSignif < 0 | nSignif > nCoeffs
        nSignif = input( ...
            sprintf('Input the number of retained coefficient (1 to %d): ', ...
                nCoeffs) );
    end;
end;
disp(' ');


% Truncate to only the nSignif most significant coefficients
nla_coeffs = pdfb_tr(coeffs, 0, 0, nSignif);

disp('Displaying the position of the retained coefficients...') ;
showpdfb( nla_coeffs, nSignif ) ;
title('Retained coefficients');
input('Press Enter key to continue...' ) ;
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pyramidal directional filter bank (PDFB) reconstruction.
% from the retained coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reconstructed image
imrec = pdfbrec( nla_coeffs, pfilter, dfilter ) ;

% Display the original image as well as the reconstructed image
subplot(1,2,1), imagesc ( im, [0, 1] ) ; 
title( sprintf('Original image (%d X %d)', size(im))) ;
axis image off;
subplot(1,2,2), imagesc( imrec, [0, 1] );
title(sprintf('Reconstructed image\n(using %d coefs; SNR = %.2f dB)', ...
    nSignif, SNR(im, imrec)));
axis image off;

disp('Comparing the original image with the NLA image by contourlets...') ;
input('Press Enter key to continue...' ) ;
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare with NLA using wavelets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wavelet transfrom using PDFB with "zero" number of DFB level
wcoefs = pdfbdec(im, pfilter, dfilter, zeros(length(nlevels), 1));

% Keep the same number of most significant coefficients as PDFB
nla_wcoefs = pdfb_tr(wcoefs, 0, 0, nSignif);
im_wrec = pdfbrec(nla_wcoefs, pfilter, dfilter);

% Only show a portion of images (size 256 x 256);
ind1 = 201:456;
ind2 = 101:356;

disp('Comparing NLA by wavelets and by contourlets...') ;

subplot(1,2,1), imagesc ( im_wrec(ind1, ind2), [0, 1] ) ; 
title( sprintf('NLA using wavelets\n(M = %d coefs; SNR = %.2f dB)', ...
    nSignif, SNR(im, im_wrec))) ;
axis image off;

subplot(1,2,2), imagesc ( imrec(ind1, ind2), [0, 1] ) ; 
title( sprintf('NLA using contourlets\n(M = %d coefs; SNR = %.2f dB)', ...
    nSignif, SNR(im, imrec))) ;
axis image off;