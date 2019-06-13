function [Index] = PCmultindex(N,P,k);
%
% PCmultindex  returns the multi-index sequence of the polynomial order
% contribution along each random dimension i=1..N of the k-th polynomial in
% the basis
%
% Synopsis:  [Index] = PCmultindex(N,P,k);
%
% Inputs:    N = order of quadrature rule
%            P = Legendre polynomial order
%            k = Polynomial number in the basis
% Output:    Index = multi-index sequence of polynomial orders
%

M = PCnumbterms(N,P);
k = k + 1;

if k == 1
    Index = zeros(N,1)';
    return
end

if N == 1
    Index = k-1;
    return
end

pv = (0:P);

Psize = zeros(1,length(pv));
for i = 1:length(pv)
    Psize(i) = PCnumbterms(N,pv(i));
end

POrder = 0;
while (k-Psize(POrder+1) > 0)
    POrder = POrder + 1;
end

PPos = k-Psize(POrder);
PTotal = POrder;

PosCurr = PPos;
for j = 1:N-1
    pv = (0:PTotal);
    dn = N - j -1;
    
    Psize = zeros(1,length(pv));
    for i = 1:length(pv)
        Psize(i) = PCnumbterms(dn,pv(i));
    end
   
    Total = 0;
    for i = 1:length(Psize)
        Total = Total + Psize(i);
        if PosCurr <= Total
            Total = Total - Psize(i);
            PosCurr = PosCurr - Total;
            break
        end
    end
    pindex = (PTotal:-1:0);
    Index(j) = pindex(i);
    PTotal = PTotal - Index(j);
end

Index(j+1) = POrder - sum(Index);

