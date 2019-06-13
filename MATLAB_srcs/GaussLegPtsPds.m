function [x,w] = GaussLegPtsPds(n)
% GaussLegPtsPds  1D Nodes in [-1,1] and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [x,w] = GaussLegPtsPds(n)
%
% Input:     n = order of quadrature rule
%
% Output:    x = vector of quadrature points
%            w = vector of quadrature weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V 