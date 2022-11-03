%% POLYNOMIAL KRYLOV

clear; clc; close all;

addpath('funmkrylov/')
fid = fopen('results.txt','w');

for k = 10
    %% Discretization
    n = 2^k;
    x = linspace(0,1,n).';
    h = x(2)-x(1);
    e  = ones(n-2,1);
    An = spdiags([-e,2*e,-e]/h^2,-1:1,n-2,n-2);
    f  = sin(pi*x);
    alpha = 0.5;
%     uexact = mpower(full(An),-alpha)*f(2:n-1);

    %% Krylov solution
    tol = 1e-9;
    maxit = n;
    [y,res,time] = MatFunLanczos(An,f(2:n-1),-2*alpha,tol,maxit);
    figure(1)
    plot(x,[0;y;0],x,f,'--')
    pause()
end
fclose(fid);