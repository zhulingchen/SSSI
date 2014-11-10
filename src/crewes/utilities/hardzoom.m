function [m2,y2,x2]=hardzoom(m1,y1,x1,lims);
% HARDZOOM: subset a matrix according to axis limit specs
%
% [m2,y2,x2]=hardzoom(m1,y1,x1,lims);
%
% Subset a matrix and its coordinate vectors according
% to an axis specification.
% m1,y1,x1 ... input matrix and its row and column coordinates
% m2,y2,x2 ... output matrix and its row and column coordinates
% lims ... vector of length four giving [xmin xmax ymin ymax]
%	(specified as in a call to axis)
%
%
indx=between(lims(1),lims(2),x1);
indy=between(lims(3),lims(4),y1);
y2=y1(indy);
x2=x1(indx);
m2=m1(indy,indx);

