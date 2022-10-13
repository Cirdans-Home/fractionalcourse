%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : cprec.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Circulant Preconditioner Apply O(n log(n))
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function y = cprec(c,x)
%%CPREC solves the preconditioning (circulant) system Cy=x.
% where x is the right-hand side vector
%       c is the vector of the eigenvalues of the circulant matrix C (see
%       function gencircprec)
%       y the solution of the system
y = ifft(fft(x)./c);
if norm(imag(y)) < 1.0e-14
y = real(y);
end
