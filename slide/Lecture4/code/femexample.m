%% FEM Example

clear; clc; %close all;

n = 32;
R = 'L';
G = numgrid(R,n);
spy(G)
title('A Finite Difference Grid')

L = delsq(G);
spy(L)
title('The 5-Point Laplacian')

N = sum(G(:)>0) % Number of interior points

rhs = ones(N,1);
if (R == 'N') % For nested dissection, turn off minimum degree ordering.
    spparms('autommd',0)
    u = L\rhs;
    spparms('autommd',1)
else
    u = L\rhs; % This is used for R=='L' as in this example
end

U1 = G;
U1(G>0) = full(u(G(G>0)));
figure(1)
clabel(contour(U1));
prism
axis square ij

figure(2)
mesh(U1)
axis([0 n 0 n 0 max(max(U1))])
axis square ij


h = 1/(n-1);
f = @(t,y) -0.05*L*y/h^2;

%% Fractional Time
y0 = u;
Jfdefun = @(t,y) -0.05*L/h^2;
%[t, y] = flmm2(alpha,fdefun,Jfdefun,t0,tfinal,y0,h,param,method,tol,itmax)
[T,Y1] = flmm2(0.3,f,Jfdefun,0,1,y0,1e-2,[],2,1e-6,100);

for i=1:length(T)
    u = Y1(:,i);
    U1 = G;
    U1(G>0) = full(u(G(G>0)));
    figure(1)
    mesh(U1)
    title(sprintf('T = %1.2e',T(i)));
    axis([0 n 0 n 0 max(max(U1))])
    axis square ij
    pause(1/25)
end

%% Ordinary Time
[T,Y2] = ode15s(f,T,y0);

for i=1:length(T)
    u = Y2(i,:);
    U1 = G;
    U1(G>0) = full(u(G(G>0)));
    figure(2)
    mesh(U1)
    title(sprintf('T = %1.2e',T(i)));
    axis([0 n 0 n 0 max(max(U1))])
    axis square ij
    pause(1/25)
end

%% Comparison

for i=1:10:length(T)
    u = Y2(i,:);
    U1 = G;
    U1(G>0) = full(u(G(G>0)));
    figure(3)
    subplot(1,3,1)
    mesh(U1)
    title(sprintf('Ordinary T = %1.2e',T(i)));
    axis([0 n 0 n 0 35])
    axis square ij
    u = Y1(:,i);
    U2 = G;
    U2(G>0) = full(u(G(G>0)));
    subplot(1,3,2)
    mesh(U2)
    title(sprintf('Fractional T = %1.2e',T(i)));
    axis([0 n 0 n 0 35])
    axis square ij
    subplot(1,3,3)
    mesh(abs(U1-U2))
    axis([0 n 0 n 0 10])
    title('Absolute Error');
    axis square ij
    set(gcf,'Color','None');
    export_fig(sprintf('diff_%03d.png',i));
    pause(1/25)
end