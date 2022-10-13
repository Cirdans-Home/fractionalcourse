%% CLUSTERING

%% THE NON SYMMETRIC CASE

clear; clc; close all;
dofov = false;

alpha = 1.5;
N = 2000;
w0 = @(x) 5*x.*(1-x);
hN = 1/(N-1); x = 0:hN:1;
dt = hN; t = 0:dt:1;
dplus = @(x,t) gamma(3-alpha).*x.^alpha;
dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
%% Discretize
G = glmatrix(N,alpha);
Gr = G; Grt = G.';
Dplus = diag(dplus(x,0));
Dminus = diag(dminus(x,0));
I = eye(N,N);
% Left-hand side
nu = hN^alpha/dt;
A = nu*I - (Dplus*Gr + Dminus*Grt);

w = w0(x).';
[ev,evt] = sunprec(N,alpha);
c = nu + mean(dplus(x,0))*ev + mean(dminus(x,0))*evt;
C = gallery('circul',ifft(c));

ev = eig(A);
evp = eig(A,C);

figure(1)
plot(real(ev),imag(ev),'x',real(evp),imag(evp),'o');
legend({'Original Eigenvalues','Preconditioned'},'FontSize',14)
title(sprintf('\\alpha = %1.1f N = %d',alpha,N))
set(gcf,'Color','none')
export_fig(sprintf("ev%d_%d.pdf",alpha*10,N))

