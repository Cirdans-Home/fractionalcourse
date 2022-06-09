%% Field of Value Example

clear; clc; close all;

load('ibm32.mat');
A = Problem.A;
b = ones(size(A,1),1);

fovals(full(A),500);

%% Build solution

alpha = 0.5;
beta = 1;
fA = ml_matrix(full(A),alpha,beta);

%% Rational Krylov
xi = inf(1,size(A,1)-1);
[V, K, H] = rat_krylov(A, b, xi);

for m=1:size(V,2)-1
   ym = V(:,1:m)*ml_matrix(full( V(:,1:m)'*A*V(:,1:m) ),alpha,beta)*V(:,1:m)'*b;
   err(m) = norm(fA*b - ym,2);
end

%% Rational Krylov
xi = 0.25*ones(1,size(A,1)-1);
[V, K, H] = rat_krylov(A, b, xi);

for m=1:size(V,2)-1
   ym = V(:,1:m)*ml_matrix(full( V(:,1:m)'*A*V(:,1:m) ),alpha,beta)*V(:,1:m)'*b;
   err_sai(m) = norm(fA*b - ym,2);
end

%% Error plot
figure(1)
semilogy(1:m,err,'o-',...
    1:m,err_sai,'x-','LineWidth',2)
