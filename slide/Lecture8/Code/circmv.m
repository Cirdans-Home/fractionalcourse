%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : circmv.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Circulant Matrix-Vector product O(nlog(n))
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ y ] = circmv( C,x )
%CIRCMV Circulant matrix-vector product
%   If C is a circulant matrix computes y = C x by the complete algorithm
%   Otherwise, if C is a vector its assume that it is the vector of the
%   eigenvalues of C and behaves consequently.
%   TEST:
%   n = 1000;
%   e = 1:n;
%   C = gallery('circul',e);
%   x = rand(n,1);
%   yc = C*x;
%   y = circmv(C,x);
%   norm(yc-y)
if ismatrix(C) && size(C,2) == size(C,1)
    c = C(:,1); 
    C = fft(c(:));
end
y = ifft(C.*fft(x));
end

