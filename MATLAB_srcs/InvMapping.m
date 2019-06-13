function [zm] = InvMapping(z,zmin,zmax)
% 
% InvMapping proceeds with a linear mapping of the z data from [zmin,zmax]
% interval to the [-1,1] range
%
% Synopsis:  [zm] = InvMapping(z,zmin,zmax);
%
% Inputs:    z = data
%            zmin = minimum bound of the range of z
%            zmax = maximum bound of the range of z
% Output:    zm = mapped data
%
zm = (2.*z - (zmax+zmin))./(zmax-zmin);