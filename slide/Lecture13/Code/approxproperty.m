%% Approximation property

clear; clc; close all;

alpha = 0.5;
n = 2^10;

x = linspace(0,1,n).';
h = x(2)-x(1);
e  = ones(n-2,1);
An = spdiags([-e,2*e,-e]/h^2,-1:1,n-2,n-2);
f  = sin(pi*x);
alpha = 0.5;
b = ones(n-2,1);
ytrue = mpower(full(An),-alpha)*b;

l = 300;
[Q,H] = arnoldi(An,b,l);
for j=1:l
    y = Q(:,1:j)*(mpower(H(1:j,1:j),-alpha)*(Q(:,1:j)'*b));
    err(j) = norm(y-ytrue)./norm(ytrue);
end

%% Convergence
figure(1)
semilogy(1:l,err);
xlabel('l')
ylabel('Relative Error')


