function G = glmatrix(N,alpha)
%%GLMATRIX produces the GL discretization of the Riemann-Liouville
%derivative
g = gl(N,alpha);
c = zeros(N,1);
r = zeros(1,N);
r(1:2) = g(2:-1:1);
c(1:N) = g(2:end);
G = toeplitz(c,r);
end