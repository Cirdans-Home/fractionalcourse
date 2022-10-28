function [y] = Mprod(A,B,Jm,h,x)
%MPROD This function computes the fast product with the matrix:
%           M = A \otims Im -h*B \otimes Jm
% by the vector x, that is y = Mx, where A and B are (s+1)x(s+1) matrices 
% and Im and Jm are m x m matrices

[sp1,~] = size(A);
[m,~] = size(Jm);

X = reshape(x,m,sp1);
Y = X*A' - h*(Jm*X*B');
y = reshape(Y,m*sp1,1);

end