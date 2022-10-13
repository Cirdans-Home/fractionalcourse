function T = invkacmatrix(n,alpha)
%INVKACMATRIX Gives back the 1/Kac-Murdock-Szego matrices
f = @(th) (1 - alpha^2)./((1-alpha*exp(1i*th)).*(1-alpha*exp(-1i*th)));
c = zeros(n,1); r = zeros(1,n);
for k=1:n
    r(k) = integral(@(th) f(th).*exp(1i*th*(k-1)),0,2*pi)/(2*pi);
    c(k) = integral(@(th) f(th).*exp(-1i*th*(k-1)),0,2*pi)/(2*pi);
end
T = real(toeplitz(r,c));
end

