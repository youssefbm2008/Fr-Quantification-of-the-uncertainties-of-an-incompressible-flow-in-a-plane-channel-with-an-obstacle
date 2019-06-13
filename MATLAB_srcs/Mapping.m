function [zm] = Mapping(z,zmin,zmax);
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
zm = (zmax-zmin).*z./2 + (zmax+zmin)./2;