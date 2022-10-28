function [ lambdamat ] = pbetagen( A,beta,omega )
%PBETAGEN Genera l'inversa approssimata Pbeta a partire dalla matrice di
%tutto il problema per la banda desiderata.
%   Detailed explanation goes here

n = size(A,1);
e = speye(n);
A = spdiags(spdiags(A,-beta:beta),-beta:beta,n,n);
%A = A'*A;
Lambdainv = zeros(n+beta-1,1);
lambdamat = zeros(n+beta-1,n);
OMEGA = omega.^((0:n+beta-2)'/(n+beta-1));
for i=1:n
    r = e(:,i)'*A;          % e_i^T A
    Ki = toeplitz(r);       % e_i^T A = e_i^T Ki
    c = [Ki(1:n,1);omega*Ki(1,beta:-1:2)'];
    Lambda = real(fft(full(OMEGA.*c)));
    for j=1:(n+beta-1)        % Compute the inverse element removing the
        if Lambda(j) > 0      % zero entries
            Lambdainv(j) = 1/Lambda(j);
        end
    end
    lambdamat(:,i) = Lambdainv;
end




end

