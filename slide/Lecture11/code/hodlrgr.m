%% HODLR GRUNWALD-LETNIKOV

clear; clc; close all;

N = 6000;
alpha = 1.5;
tol = 1e-13;
G = glhodlrmatrix(N,alpha,tol);

figure(1)
spy(G)
axis square
set(gcf,'Color','white')

BN = @(N,tol) 2 + 2*ceil( (2/pi^2).*log(4*N/pi).*log(32/tol) );
fprintf("qsrank(%d,%e)=%d\n",N,tol,BN(N,tol));
