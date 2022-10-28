function [ omega ] = computeoptimalomega( T )
%COMPUTEOPTIMALOMEGA compute the optimal omega for an omega circulant
%preconditioner given a Toeplitz matrix.
%   T := input toeplitz matrix

n = size(T,1);
tc = T(2:end,1)';
tr = T(1,end:-1:2);
j = 1:n-1;
nmj = n-j;
omega = -angle(sum(nmj.*j.*conj(tc).*tr));

end

