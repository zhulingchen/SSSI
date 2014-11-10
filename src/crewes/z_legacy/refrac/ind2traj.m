function indtraj=ind2traj(irow,icol,nrows)
%
% indtraj=ind2traj(irow,icol,nrows)
%
% IND2TRAJ converts a set of row and column indicies to a set
% equivalent 1-D indicies which evaluate the matrix along the
% trajectory specified by the row/column indicies. Given irow
% as a vector and icol as a vector, direct evaluation of a 
% A as: B=A(irow,icol), will result in B as a matrix which
% has dimensions length(irow) by length(icol). However, there
% is an alternate interpretation of irow and icol as a set of
% index pairs which specify a trajectory through A. In this
% case, we would like to write something like: B=A(indtraj)
% and have B result as a vector which is the matrix A sliced
% along the trajectory. This routine does jsut that.
%
%	irow ... vector of row indicies of the trajectory
%	icol ... vector of column indicies of the trajectory
%	nrows ... number of rows of the matrix to be sliced
%	indtraj ... index of the trajectory
%
% example: A=magic(3);
%	indtraj=ind2traj(1:3,1:3,3);
%	A(indtraj) % returns the main diagonal of A
%
% G.F. Margrave, University of Calgary, 1996
%
if(length(icol)~=length(irow))
	error('irow and icol must be the same length');
end
indtraj = (icol(:)-1)*nrows + irow(:);
