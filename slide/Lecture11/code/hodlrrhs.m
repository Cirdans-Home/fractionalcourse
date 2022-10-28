%% FUNCTION WITH HODLR RHS

clear; clc; close all;

N = 2000;
x = linspace(0,1,N);
y = linspace(0,1,N);
[X,Y] = meshgrid(x,y);
tau = 1;

C = log(tau + abs(X-Y));
figure(1)
mesh(X,Y,C)
xlabel('x')
ylabel('y')
zlabel('f(x,y)')

hC = hodlr(C);
figure(2)
spy(hC)
axis square
set(gcf,"Color","none")