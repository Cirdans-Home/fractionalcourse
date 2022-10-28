function [ y ] = hankenagy( Lambdainv,x,n,beta,omega )
%HANKENAGY Hanke and Nagy inverse Toeplitz preconditioner 
%   Apply the Hanke and Nagy inverse Toeplitz preconditioner to the vector
%   x, i.e., y = Px, n and beta are the true size of the system and its
%   band amplitude.

r = [x;zeros(beta-1,1)];
OMEGA = omega.^((0:n+beta-2)'/(n+beta-1));
y = conj(OMEGA).*ifft(Lambdainv.*fft(OMEGA.*r));
y = real(y(1:n));


end

