function [X,Y] = GetMeshCoords(meshfile,Nx,Ny)
% 
% Mapping proceeds with a linear mapping of the z data from [-1,1]
% interval to the  [zmin,zmax] range
%
% Synopsis:  [zm] = Mapping(z,zmin,zmax);
%
% Inputs:    z = data
%            zmin = minimum bound of the range of zm
%            zmax = maximum bound of the range of zm
% Output:    zm = mapped data
%
load(meshfile); X6 = xy_001(:,1); Y6 = xy_001(:,2); clear xy_001;
X = zeros(2*Nx-1,2*Ny-1); Y = zeros(2*Nx-1,2*Ny-1);
X = reshape(X6,2*Ny-1,2*Nx-1)'; Y = reshape(Y6,2*Ny-1,2*Nx-1)'; clear X6 Y6;