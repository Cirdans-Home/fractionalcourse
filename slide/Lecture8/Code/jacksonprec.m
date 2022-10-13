%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : jacksonprec.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Construction of Jackson-kernel based circulant preconditioner
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ c ] = jacksonprec( t,r )
%JACKSONPREC Circulant preconditioner based on the Jackson Kernel

n = length(t);            % size of the preconditioner
t1 = conj(t(n:-1:2));     % is the last column of T_n(f)
if r == 2 || r == 3 || r == 4
    coef = convol(n,r).';
    c = [t(1)*coef(1) (coef(2:n).*t(2:n)+coef(n:-1:2).*t1).'];
    c = fft(c)';
else
    error('r needs to be 2, 3 or 4');
end

c = real(c);

end

function coef = convol(n,r)
% convol computes the coefficients for the generalized Jackson kernel.
% K_{m,2r}, see (4.2).
% n: the size of the Toeplitz matrix T;
% r: 2, 3, or 4;
% coef: the coefficients {b_kˆ(m,r):|k|<=r(m-1)} of
% the generalized Jackson kernel K_{m,2r}.

m = floor(n/r);
a = 1:-1/m:1/m;
r0 = 1;
coef = [a(m:-1:2) a];
while r0 < r
    M = (2*r0+3)*m;
    b1 = zeros(M,1);
    c = zeros(M,1);
    c(1:m) = a;
    c(M:-1:M-m+2) = a(2:m);
    b1(m:m+2*r0*(m-1)) = coef;
    tp = ifft(fft(b1).*fft(c));
    coef = real(tp(1:2*(r0+1)*(m-1)+1));
    r0 = r0+1;
end
M = r*(m-1)+1;
coef = [coef(M:-1:1)' zeros(1,n-M)]';
coef = coef';
end
