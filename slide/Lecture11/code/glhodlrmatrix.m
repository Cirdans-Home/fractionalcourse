function G = glhodlrmatrix(N,alpha,tol)
%%GLMATRIX produces the GL discretization of the Riemann-Liouville
%derivative in HODLR format
g = gl(N,alpha);
c = zeros(N,1);
r = zeros(1,N);
r(1:2) = g(2:-1:1);
c(1:N) = g(2:end);

hodlroption( 'threshold', tol);
G = hodlr('toeplitz',c,r);
end