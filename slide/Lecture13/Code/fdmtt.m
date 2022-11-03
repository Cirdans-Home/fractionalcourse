%% Finite Difference Matrix-Transfer Techniques

clear; clc; close all;

n = 100;
x = linspace(0,1,n).';
h = x(2)-x(1);
e  = ones(n-2,1);
An = spdiags([-e,2*e,-e]/h^2,-1:1,n-2,n-2);
f  = sin(pi*x);
u = An\f(2:n-1);
u = [0;u;0];
[U,L,V] = eig(full(An));
ualpha = @(alpha) [0;V'\(L.^alpha\(U\f(2:n-1)));0];


plot(x,f,'--', ...
    x,u,...
    x,ualpha(0.1),...
    x,ualpha(0.5),...
    x,ualpha(0.9),...
    "LineWidth",2)
legend('Right-Hand Side','Laplacian','\alpha = 0.1','\alpha = 0.5','\alpha = 0.9')

