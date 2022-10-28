function [ Lambdainv ] = gennagy( Ki,beta,omega,n )
%GENNAGY Summary of this function goes here
%   Detailed explanation goes here

Lambdainv = zeros(n+beta-1,1);
OMEGA = omega.^((0:n+beta-2)'/(n+beta-1));
c = [Ki(1:n,1);omega*Ki(1,beta:-1:2)'];
Lambda = real(fft(full(OMEGA.*c)));
for j=1:(n+beta-1)        % Compute the inverse element removing the
    if Lambda(j) > 0      % zero entries
        Lambdainv(j) = 1/Lambda(j);
    end
end

end

