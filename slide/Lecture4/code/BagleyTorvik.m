%% Bagley-Torvik

clear; clc; close all;

m = 10; A = 6; K = 3;
mu = 2; rho = 2;
alpha = [2 3/2] ;
lambda = [m 2*A*sqrt(mu*rho)] ;
f_fun = @(t,X) -K*X;
J_fun = @(t,X) -K;
t0 = 0 ; T = 100 ;
X0 = [0 , 2 ];
h = 1e-2;
[t, X] = mt_fde_pi1_ex(alpha,lambda,f_fun,t0,T,X0,h);

figure(1)
plot(t,X,'-','LineWidth',2)