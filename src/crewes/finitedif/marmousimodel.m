function [vel,x,z]=marmousimodel(dx)
% MARMOUSIMODEL ... return the 2D p-wave Marmousi model
%
% [vel,x,z]=marmousimodel(dx)
%
% dx ... must be either 5, 10 or 25. This specifies the grid size in both x
%       and z for the output model. 
% vel ... matrix of p-wave velocities
% x ... x (column) coordinate for vel
% z ... z (row) coordinate for vel
% 
% You can plot the result with
% figure;
% imagesc(x,z,vel);colorbar
%
if(dx~=5 && dx~= 10 && dx~=25)
    error('only dx of 5, 10, and 25 are available')
end
vel=[];
x=[];
z=[];
if(dx==25)
    s=which('raymarmousi_demo');
    if(isempty(s))
        error('25 m Marmousi model not found, you need to load and install the CREWES toolbox')
    end
    ind = strfind(s,'raymarmousi_demo');
    sm=[s(1:ind-1) 'marmousi_mod'];
    disp(['Marmousi model loaded from ' sm])
    load(sm)
    ind=find(x>9200);
    if(~isempty(ind))
        x(ind)=[];
        vel(:,ind)=[];
    end
elseif(dx==10)
    s=which('afd_snap');
    if(isempty(s))
        error('10m Marmousi model not found, you need to load and install the CREWES toolbox')
    end
    ind = strfind(s,'afd_snap');
    sm=[s(1:ind-1) 'marmousi_dz10'];
    disp(['Marmousi model loaded from ' sm])
    load(sm)
elseif(dx==5)
    s=which('afd_snap');
    if(isempty(s))
        error('5m Marmousi model not found, you need to load and install the CREWES toolbox')
    end
    ind = strfind(s,'afd_snap');
    sm=[s(1:ind-1) 'marmousi_dz5'];
    disp(['Marmousi model loaded from ' sm])
    load(sm)
end
%make sure first x is zero
x=x-x(1);