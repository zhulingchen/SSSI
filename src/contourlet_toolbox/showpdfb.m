function displayIm = showpdfb(y, scaleMode, displayMode, ...
                              lowratio, highratio, coefMode, subbandgap)
% SHOWPDFB   Show PDFB coefficients. 
%
%       showpdfb(y, [scaleMode, displayMode, ...
%                    lowratio, highratio, coefMode, subbandgap])
%
% Input:
%	y:	a cell vector of length n+1, one for each layer of 
%		subband images from DFB, y{1} is the lowpass image
%
%   scaleMode: 
%       scale mode (a string or number):
%           If it is a number, it denotes the number of most significant 
%           coefficients to be displayed.  Its default value is 'auto2'.
%           'auto1' ---   All the layers use one scale. It reflects the real 
%                         values of the coefficients.
%                         However, the visibility will be very poor.
%           'auto2' ---   Lowpass uses the first scale. All the highpass use 
%                         the second scale.
%           'auto3' ---   Lowpass uses the first scale. 
%                         All the wavelet highpass use the second scale.
%                         All the contourlet highpass use the third scale.
%   displayMode: 
%       display mode (a string): 
%           'matlab' --- display in Matlab environment. 
%                        It uses the background color for the marginal
%                        image.
%           'others' --- display in other environment or for print.
%                        It used the white color the marginal image. 
%                        It is the default value.
%   lowratio:
%       display ratio for the lowpass filter (default value is 2).
%       It ranges from 1.2 to 4.0.
%   highratio:
%       display ratio for the highpass filter (default value is 6).
%       It ranges from 1.5 to 10.
%  coefMode: 
%       coefficients mode (a string): 
%           'real' ----  Highpass filters use the real coefficients. 
%           'abs' ------ Highpass filters use the absolute coefficients. 
%                        It is the default value
%   subbandgap:	
%       gap (in pixels) between subbands. It ranges from 1 to 4.
%
% Output:
%	displayIm: matrix for the display image.
%   
% See also:     PDFBDEC, DFBIMAGE, COMPUTESCALE

% History:
%   09/17/2003  Creation.
%   09/18/2003  Add two display mode, denoted by 'displayMode'.
%               Add two coefficients mode, denoted by 'coeffMode'.
%   10/03/2003  Add the option for the lowpass wavelet decomposition.
%   10/04/2003  Add a function computescale in computescales.m. 
%               This function will call it.
%               Add two scal modes, denoted by 'scaleMode'. 
%               It can also display the most significant coefficients.
%   10/05/2003  Add 'axis image' to control resizing.
%               Use the two-fold searching method to find the 
%               background color index.

% Scale mode
if ~exist('scaleMode', 'var')
    scaleMode = 'auto2' ;
elseif isnumeric( scaleMode ) % Denote the number of significant coefficients to be displayed
    if scaleMode < 2
        display( 'Warning! The numbe of significant coefficients must be positive!' ) ;
        scaleMode = 50 ;
    end
elseif ~strcmp(scaleMode,'auto1') & ~strcmp(scaleMode, 'auto2') & ~strcmp(scaleMode, 'auto3')
    display ('Warning! There are only two scaleMode mode: auto1, auto2, auto3! Its defualt value is "auto2"!');
    scaleMode = 'auto2' ;
end

% Display ratio for the lowpass band
if ~exist('lowratio', 'var')
    lowratio = 2 ;
elseif highratio < 1
    display ('Warning! lowratio must be larger than 1!Its defualt value is 2!');
end

% Display ratio for the hiphpass band
if ~exist('highratio', 'var')
    highratio = 6 ;
elseif highratio < 1
    display ('Warning! highratio must be larger than 1! Its defualt value is 6!');
end

% Gap between subbands
if ~exist('subbandgap', 'var')
    subbandgap = 1;  
elseif subbandgap < 1
    display ('Warning! subbandgap must be no less than 1! Its defualt value is 1!');
    subbandgap = 1;
end

% Display mode
if ~exist('displayMode', 'var')
    displayMode = 'others' ;
elseif ~strcmp(displayMode,'others') & ~strcmp(displayMode, 'matlab')
    display ('Warning! There are only two display mode: matlab, others! Its defualt value is "others"!');
    displayMode = 'others' ;
end

% Coefficient mode
if ~exist('coefMode', 'var')
    coefMode = 'abs' ;
elseif ~strcmp(coefMode,'real') & ~strcmp(coefMode, 'abs')
    display ('Warning! There are only two coefficients mode: real, abs! Its defualt value is "abs"!');
    coefMode = 'abs' ;
end

% Parameters for display
layergap = 1 ; % Gap between layers

% Input structure analysis. 
nLayers = length(y); %number of PDFB layers
% Compute the number of wavelets layers.
% We assume that the wavelets layers are first several consecutive layers.
% The number of the subbands of each layer is 3.
fWaveletsLayer = 1;
nWaveletsLayers = 0 ; %Number of wavelets layers.
nInxContourletLayer = 0 ; %The index of the first contourlet layer.
i = 2 ;
while fWaveletsLayer > 0 & i <= nLayers
    if length( y{i} ) == 3
        nWaveletsLayers = nWaveletsLayers + 1 ;
    else
        fWaveletsLayer = 0 ;
    end;
    i = i + 1 ;
end;
nInxContourletLayer = 2 + nWaveletsLayers ;
    
% Initialization 
% Since we will merge the wavelets layers together, 
% we shall decrease the number of display layers.
nDisplayLayers = nLayers - nWaveletsLayers ;
cellLayers = cell (1, nDisplayLayers); % Cell for multiple display layers
vScalesTemp = zeros (1, 2) ;  % Temporary scale vector.
vScales = zeros (nLayers, 2); % Scale vectors for each layer 
nAdjustHighpass = 2 ; % Adjustment ratio for the highpass layers.


if ~isnumeric( scaleMode )
    switch( scaleMode )% Compute the scales for each layer
    case 'auto1'
        vScalesTemp = computescale( y, lowratio, 1, nLayers, coefMode ) ;
        for i = 1 : nLayers
            vScales( i, :) = vScalesTemp ;
        end
    case 'auto2'
        vScales( 1, :) = computescale( y, lowratio, 1, 1, coefMode ) ;
        vScalesTemp = computescale( y, highratio, 2, nLayers, coefMode ) ;
        % Make a slight adjustment. Compared to the lowpass, the highpass shall be insignificant.
        % To make the display more realistic, use a little trick to make the upper bound a little bigger. 
        vScalesTemp (2) = nAdjustHighpass * ( vScalesTemp (2) - vScalesTemp (1) ) + vScalesTemp (1) ;
        for i = 2 : nLayers
            vScales( i, :) = vScalesTemp ;
        end
    case 'auto3'
        vScales( 1, :) = computescale( y, lowratio, 1, 1, coefMode ) ;
        vScalesTemp = computescale( y, highratio, 2, 1+nWaveletsLayers, coefMode ) ;
        % Make a slight adjustment. Compared to the lowpass, the highpass shall be insignificant.
        % To make the display more realistic, use a little trick to make the upper bound a little bigger. 
        vScalesTemp (2) = nAdjustHighpass * ( vScalesTemp (2) - vScalesTemp (1) ) + vScalesTemp (1) ;
        for i = 2 : nWaveletsLayers + 1
            vScales( i, :) = vScalesTemp ;
        end
        vScalesTemp = computescale( y, highratio, nInxContourletLayer, nLayers, coefMode ) ;
        % Make a slight adjustment. Compared to the lowpass, the highpass shall be insignificant.
        % To make the display more realistic, use a little trick to make the upper bound a little bigger. 
        vScalesTemp (2) = nAdjustHighpass * ( vScalesTemp (2) - vScalesTemp (1) ) + vScalesTemp (1) ;
        for i = nInxContourletLayer : nLayers
            vScales( i, :) = vScalesTemp ;
        end
    otherwise % Default value: 'auto2'.
        vScales( 1, :) = computescale( y, lowratio, 1, 1, coefMode ) ;
        vScalesTemp = computescale( y, highratio, 2, nLayers, coefMode ) ;
        for i = 2 : nLayers
            vScales( i, :) = vScalesTemp ;
        end
    end
    % Verify that they are reasonable
    for i = 1: nLayers
        if vScales (i, 2) < vScales (i, 1) + 1.0e-9
            display ('Error! The scale vectors are wrong! Exit!' ) ;
            exit ;
        end;
    end;
    %display ( vScales ) ;
else % Compute the threshold for the display of coefficients
    % Convert the output into the vector format
    [vCoeff, s] = pdfb2vec(y);
    
    % Sort the coefficient in the order of energy.
    vSort = sort( abs( vCoeff ));
    clear vCoeff;
    vSort = fliplr(vSort);
    
    % Find the threshold value based on number of keeping coeffs
    dThresh = vSort( scaleMode );
    clear vSort;
end

% Prepare for the display
colormap(gray);
cmap = get(gcf,'Colormap');
cColorInx = size(cmap,1);

%  Find background color index:
if strcmp( displayMode, 'matlab' )
    % Get the background color (gray value)
    dBgColor = get( gcf, 'Color' ) ;
    
    % Search the color index by 2-fold searching method.
    % This method is only useful for the gray color!
    nSmall = 1 ;
    nBig = cColorInx ;
    while nBig > nSmall + 1
        nBgColor = floor ((nSmall + nBig) / 2) ;
        if dBgColor(1) < cmap (nBgColor, 1)
            nBig = nBgColor ;
        else
            nSmall = nBgColor ;
        end
    end;
    if abs( dBgColor(1) - cmap (nBig, 1) ) > abs ( dBgColor(1) - cmap( nSmall, 1) )
        nBgColor = nSmall ;
    else
        nBgColor = nBig ;
    end
end

% Merge all layers to corresponding display layers.
% Prepare the cellLayers, including the boundary.
% Need to polish with real boudary later!
% Now we add the boundary, but erase the images!!
% First handle the lowpass filter
% White line around subbands

% 1. One wavelets layers. 
gridI = cColorInx - 1 ;
cell4Wavelets = cell(1, 4) ; %Store 4 wavelets subbands.
if isnumeric ( scaleMode ) %Keep the significant efficients
    waveletsIm = cColorInx * double(abs(y{1}) >= dThresh) ;
else
    dRatio = (cColorInx-1) / (vScales(1,2)-vScales(1,1));
    if strcmp( coefMode, 'real' )
        waveletsIm = double( 1 + (y{1}-vScales(1,1))*dRatio ) ;
    else
        waveletsIm = double( 1 + (abs(y{1})-vScales(1,1))*dRatio ) ;
    end
end
% Merge other wavelets layers
if nWaveletsLayers > 0
    for i=2 : nWaveletsLayers + 1
        cell4Wavelets{1} = waveletsIm ; 
        % Compute with the scale ratio.
        if ~isnumeric( scaleMode )
            dRatio = (cColorInx-1) / (vScales(i,2)-vScales(i,1));
        end
        m = length(y{i});
        if m ~= 3
            display('Error! Incorect number of wavelets subbands! Exit!');
            exit;  
        end
        for k = 1:m
            if isnumeric ( scaleMode ) %Keep the significant efficients
                cell4Wavelets{k+1} = cColorInx * double(abs(y{i}{k}) >= dThresh) ;
            else
                if strcmp( coefMode, 'real' )
                    cell4Wavelets{k+1} = double( 1 + (y{i}{k}-vScales(i,1))*dRatio );
                else
                    cell4Wavelets{k+1} = double( 1 + (abs(y{i}{k})-vScales(i,1) )* dRatio );
                end
            end
        end
        waveletsIm = dfbimage(cell4Wavelets, subbandgap, gridI);
    end
end
cellLayers{1} = waveletsIm ;
nHeight = size( cellLayers{1}, 1 );


% 2. All the contourlet layers
for i = nInxContourletLayer : nLayers
    % Compute with the scale ratio.
    if ~isnumeric( scaleMode )
        dRatio = (cColorInx-1) / (vScales( i, 2)-vScales( i, 1));
    end
    m = length(y{i});
    z = cell(1, m);
    for k = 1:m
        if isnumeric ( scaleMode ) %Keep the significant efficients
            z{k} = cColorInx * double(abs(y{i}{k}) >= dThresh) ;
        else
            if strcmp( coefMode, 'real' )
                z{k} = double( 1 + (y{i}{k}-vScales(i,1)) * dRatio );
            else
                z{k} = double( 1 + (abs(y{i}{k})-vScales(i,1) )* dRatio );
            end
        end
    end
    cellLayers{i-nWaveletsLayers} = dfbimage(z, subbandgap, gridI);
    nHeight = nHeight + size(cellLayers{i-nWaveletsLayers}, 1);
end
% Compute the width of the dispaly image.
nWidth = size(cellLayers{nDisplayLayers}, 2);

% Merge all layers and add gaps between layers
nHeight = nHeight + layergap * (nDisplayLayers - 1) ;
% Set the background for the output image
if strcmp( displayMode, 'matlab' )
    displayIm = nBgColor * ones( nHeight, nWidth);
else
    displayIm = (cColorInx-1) * ones( nHeight, nWidth);
end
nPos = 0; %output image pointer
for i = 1 : nDisplayLayers
    [h, w] = size( cellLayers{i} );
    displayIm( nPos+1: nPos+h, 1:w) = cellLayers{i};
    if i < nDisplayLayers
        % Move the position pointer and add gaps between layers
        nPos = nPos + h + layergap ;
    end        
end

hh = image( displayIm );
% title('decompostion image');
axis image off;