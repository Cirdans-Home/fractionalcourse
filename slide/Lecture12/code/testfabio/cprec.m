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