%% TEST WITH BDF

clear; clc; close all;

L = [-2, 1; 0, -3];
y0 = [1;1];
n = length(L);
m = 100;
T = linspace(0,10,m); h = T(2)-T(1);

%% ODE Solver
tic;
[tt,yy] = ode15s(@(t,y) L*y,T,y0);
timeode = toc;

%% BUILDING THE ALL-AT-ONCE SOLUTION
r = zeros(m-1,1);
c = zeros(m-1,1);
r(1:7) = [147/60,-6,15/2,-20/3,15/4,-6/5,1/6];
c(1) = 147/60;
A = toeplitz(r,c);
% Fix boundary conditions
A(1,1) = 1;
A(2,1) = -2; A(2,2) = 3/2;
A(3,1) = 3/2; A(3,2) = -3; A(3,3) = 11/6;
A(4,1) = -4/3; A(4,2) = 3; A(4,3) = -4; A(4,4) = 25/12;
A(5,1) = 5/4; A(5,2) = -10/3; A(5,3) = 5; A(5,4) = -5; A(5,5) = 137/60;

In = speye(n,n);
Im = speye(m-1,m-1);

%% SOLVING
M = kron(A,In) - h*kron(Im,L);
b = zeros((m-1)*n,1); 
b(1:2) = y0;
b(3:4) = -1/2*y0;
b(5:6) = 1/3*y0;
b(7:8) = -1/4*y0;
b(9:10) = 1/5*y0;
b(11:12) = -1/6*y0;
tic;
x = M\b;
timeall=toc;
%% DID I SOLVE ANYTHING?
X = reshape(x,n,m-1);
X = [y0,X];

figure(1)
plot(T,X(1,:),'r-',T,X(2,:),'b-',...
    T,yy(:,1),'ro',...
    T,yy(:,2),'bo');