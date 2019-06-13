function [Up] = PCSampling(Um,Pol,Ne,T,M);
% 
% PCSampling computes Ne individual time T (or space X) random realizations
% given the first M random gPC coefficients Um(t/x) and P polynomials
% Psi(xi) of the process; and
% 
% Synopsis:  [Up] = PCSampling(Um,Pol,Ne,T,M);
%
% Inputs:    Um = gPC coefficients; size(Um) = [T M]
%            Pol = gPC polynomials evaluated at size(Pol) = [Ne M]
%            Ne = number of random dimensions
%            T = 
%            M = total number of gPC coefficients
% Output:    Up = realisations du processus stochastique; size(Up) = [Ne T]
%

[r,c]  = size(Pol);
[rr,cc] = size(Um);

if (r~=Ne || c~=M)
    error('The gPC polynomial must have dimension [N,M]');
end

if (rr~=T || c~=M)
    error('The gPC random coefficients must have dimension [T,M]');
end


Um  = repmat(Um(1:T,1:M),[1 1 Ne]);
Pol = repmat(Pol(1:Ne,1:M),[1 1 T]);
Um  = permute(Um,[3 1 2]);
Pol = permute(Pol,[1 3 2]);

Up = sum(Um.*Pol,3); clear Um Pol;


