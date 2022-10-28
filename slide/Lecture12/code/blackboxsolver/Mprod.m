function [y] = Mprod(A,B,L,h,x)
%MPROD This function computes the fast product with the matrix:
%           M = A \otims Im -h*B \otimes L
% by the vector x, that is y = Mx, where A and B are m x m matrices 
% and In and Ln are n x n matrices

[m,~] = size(A);
[n,~] = size(L);

X = reshape(x,n,m);
Y = X*A' - h*(L*X*B');
y = reshape(Y,n*m,1);

end