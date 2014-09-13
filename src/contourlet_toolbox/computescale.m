function vScales = computescale ( cellDFB, dRatio, nStart, nEnd, coefMode )
% COMPUTESCALE   Comupute display scale for PDFB coefficients
%
%       computescale(cellDFB, [dRatio, nStart, nEnd, coefMode])
%
% Input:
%	cellDFB:	a cell vector, one for each layer of 
%		subband images from DFB.
%   dRatio:
%       display ratio. It ranges from 1.2 to 10.
%
%   nStart:
%       starting index of the cell vector cellDFB for the computation. 
%       Its default value is 1.
%   nEnd:
%       ending index of the cell vector cellDFB for the computation. 
%       Its default value is the length of cellDFB.
%   coefMode: 
%       coefficients mode (a string): 
%           'real' ----  Highpass filters use the real coefficients. 
%           'abs' ------ Highpass filters use the absolute coefficients. 
%                        It's the default value
% Output:
%	vScales ---- 1 X 2 vectors for two scales.
%
% History: 
%   10/03/2003  Creation.
%
% See also:     SHOWPDFB

if ~iscell(cellDFB)
    display ('Error! The first input must be a cell vector! Exit!');
    exit;
end

% Display ratio
if ~exist('dRatio', 'var')
    dRatio = 2 ;
elseif dRatio < 1
    display ('Warning! the display ratio must be larger than 1!Its defualt value is 2!');
end

% Starting index for the cell vector cellDFB
if ~exist('nStart', 'var')
    nStart = 1 ;
elseif nStart < 1 | nStart > length(cellDFB)
    display ('Warning! The starting index from 1 to length(cellDFB)! Its defualt value is 1!');
    nStart = 1 ;
end

% Starting index for the cell vector cellDFB
if ~exist('nEnd', 'var')
    nEnd = length(cellDFB) ;
elseif nEnd < 1 | nEnd > length(cellDFB)
    display ('Warning! The ending index from 1 to length(cellDFB)! Its defualt value is length(cellDFB)!');
    nEnd = length( cellDFB ) ;
end

% Coefficient mode
if ~exist('coefMode', 'var')
    coefMode = 'abs' ;
elseif ~strcmp(coefMode,'real') & ~strcmp(coefMode, 'abs')
    display ('Warning! There are only two coefficients mode: real, abs! Its defualt value is "abs"!');
    coefMode = 'abs' ;
end

% Initialization
dSum = 0 ;
dMean = 0 ;
dAbsSum = 0 ;
nCount = 0 ;
vScales = zeros (1, 2);

if strcmp( coefMode, 'real') %Use the real coefficients
    % Compute the mean of all coefficients
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) % Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m
                dSum = dSum + sum( sum( cellDFB{i}{j} ));
                nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            dSum = dSum + sum( sum( cellDFB{i} ));
            nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	if nCount < 2 | dAbsSum < 1e-10
        display('Error! No data in this unit! Exit!');
        exit;
	else
        dMean = dSum / nCount ;
    end
    
    % Compute the STD.
    dSum = 0 ;
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) %Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m
                dSum = dSum + sum( sum( (cellDFB{i}{j}-dMean).^2 ));
                nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            dSum = dSum + sum( sum( (cellDFB{i}-dMean).^2 ));
            nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	if nCount < 2 | dSum < 1e-10
        display('Error! No data in this unit! Exit!');
        exit;
	else
        dStd = sqrt( dSum / (nCount-1) );
	end
	vScales( 1, 1 ) = dMean - dRatio * dStd ;
    vScales( 2, 2 ) = dMean + dRatio * dStd ;
    
else %Use the absolute coefficients
    % Compute the mean of absolute values
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) % Check whether it is a cell vector
            m = length(cellDFB{i});
            for j=1 : m
                dAbsSum = dAbsSum + sum( sum( abs(cellDFB{i}{j}) ));
                nCount = nCount + prod( size( cellDFB{i}{j} ) );
            end
        else
            dAbsSum = dAbsSum + sum( sum( abs(cellDFB{i}) ));
            nCount = nCount + prod( size( cellDFB{i} ) ); 
        end
    end
	if nCount < 2 | dAbsSum < 1e-10
        display('Error! No data in this unit! Exit!');
        exit;
	else
        dAbsMean = dAbsSum / nCount ;
    end
    
    % Compute the std of absolute values
    dSum = 0 ;
    for i= nStart : nEnd
        if iscell( cellDFB{i} ) %Check whether it is a cell vector
            m = length( cellDFB{i} );
            for j = 1 : m
                dSum = dSum + sum( sum( (abs(cellDFB{i}{j})-dAbsMean).^2 ));
            end
        else
            dSum = dSum + sum( sum( (abs(cellDFB{i})-dAbsMean).^2 ));  
        end
    end
    dStd = sqrt( dSum / (nCount-1) ); 
    
    % Compute the scale values
    if dAbsMean - dRatio * dStd > 0
        vScales( 1 ) = dAbsMean - dRatio * dStd ;
    else % The absolute value will be nonnegative.
        vScales( 1 ) = 0 ;
    end
    vScales( 2 ) = dAbsMean + dRatio * dStd ;	
end
