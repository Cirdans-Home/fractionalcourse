function T = kacmatrix(n,alpha)
%KACMATRIX Gives back the Kac-Murdock-Szego matrices
e = ones(n,1);
T = spdiags(([-alpha,1+alpha^2,-alpha]./(1-alpha^2)).*e,-1:1,n,n);
end

