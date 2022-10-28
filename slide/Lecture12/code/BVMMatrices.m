%% BVM Matrices

clear; clc; close all;

m = 200;      % Number of time steps
k = 8;       % Number of GBDF steps

% Discretization of the time-dependent part:
[Alpha,Beta] = mab(k,m);
A = Alpha(:,2:m+1);
B = Beta(:,2:m+1);

[f,lineSegs] = fov(full(A));
ev = eig(full(A));

figure(1)
hold on
plot(f,'k-','LineWidth',2);
plot(lineSegs,'k-','LineWidth',2);
plot(real(ev),imag(ev),'rx');
hold off
hline(0,'k--');
vline(0,'k--');
axis square tight