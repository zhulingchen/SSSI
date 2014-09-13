function [SX,SY,FX,FY,NX,NY] = fdct_wrapping_param(C)

% fdct_wrapping_param - Gives the phase space location of each curvelet
%
% Input
%     C     Curvelet coefficients
%
% Output
%     SX    cell array, SX{s}{w}(i,j), gives the row index of the center of the curvelet
%           index by (i,j)  in scale s and wedge w 
%     SY    cell array, SY{s}{w}(i,j), gives the column index of the center of the curvelet
%           index by (i,j)  in scale s and wedge w 
%     FX    cell array, FX{s}{w} gives the row index of the center frequency of the curvelet
%           in scale s and wedge w
%     FY    cell array, FY{s}{w} gives the column index of the center frequency of the curvelet
%           in scale s and wedge w
%     NX    cell array, NX{s}{w} gives the number of rows of the curvelet matrix
%           in scale s and wedge w
%     FY    cell array, FY{s}{w} gives the number of columns of the curvelet matrix
%           in scale s and wedge w
%
  
  [m,n] = size(C{end}{1});
  nbscales = floor(log2(min(m,n)))-3;
  nbangles_coarse = 16;
  allcurvelets = 0;
  
  [SX,SY,FX,FY,NX,NY] = fdct_wrapping_param_mex(m, n, nbscales, nbangles_coarse, allcurvelets);
  
  %rescale
  for s=1:nbscales
    for w=1:length(FX{s})
      nx = NX{s}{w};      ny = NY{s}{w};
      cx = ceil((nx+1)/2);      cy = ceil((ny+1)/2);
      sx = SX{s}{w};      sy = SY{s}{w};
      [IX, IY] = ndgrid(1:nx, 1:ny);
      SX{s}{w} = 1 + m*(sx*(IX-cx)+0.5);
      SY{s}{w} = 1 + n*(sy*(IY-cy)+0.5);
    end
  end
  
  
