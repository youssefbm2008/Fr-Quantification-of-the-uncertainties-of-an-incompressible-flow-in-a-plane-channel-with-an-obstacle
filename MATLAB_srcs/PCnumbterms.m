
function [M]= PCnumbterms(P,N)
% 
% PCnumbterms returns the number of terms M in the polynomial chaos
% developpement of dimensions N and total order at max P
%
% Synopsis:  [M] = PCnumbterms(N,P);
%
% Inputs:    N = number of random variables
%            P = maximum polynomial order
% Output:    M = number of terms in PC expansion
%
% Remark:   Robust implementation of factorial(N+P)/factorial(N)/factorial(P)

  M = 0;
  for s=1:P,
      tmp = 1.0;
      for r=0:s-1,
        tmp = tmp.*(N+r)/(r+1.0);
      end
      M = M + floor(tmp+0.1);
  end
  M = M+1;
