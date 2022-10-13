%% Inverse of a Toeplitz matrix
clear;
N = 10;
e = ones(N,1);
T = spdiags([-e,2*e,-e],-1:1,N,N);

heatmap(T)
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
Ax.ColorbarVisible = 'off';
set(gcf,'Color','none')
export_fig 'originalmat.png'

heatmap(inv(T))
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
Ax.ColorbarVisible = 'off';
set(gcf,'Color','none')
export_fig 'invmat.png'

%% Symmetric distribution

N = 100; alpha = 1.3;
hN = 1/(N-1); dt = hN;
G = glmatrix(N,alpha);
I = eye(N,N);
nu = dt;
A = nu^(alpha-1)*I-0.5*(G+G');
ev = eig(A);
f = @(t) -exp(-1i*t).*(1-exp(1i*t)).^alpha;
p = @(t) nu^(alpha-1) + 0.5*(f(t) + conj(f(t)));
t = linspace(-pi,pi,N);
figure(2)
plot(t,ev,'o',t,sort(p(t),'ascend'),'-')
axis tight

%% CG failing

clear; clc; close all;

for N = [100,200,500,1000,5000,10000]

    alpha = 1.2;
    hN = 1/(N-1); dt = hN;
    G = glmatrix(N,alpha);
    I = eye(N,N);
    nu = dt;
    A = nu^(alpha-1)*I-0.5*(G+G');
    b = nu^(alpha-1)*ones(N,1);
    [x,flag,relres,iter,resvec] = pcg(A,b,1e-6,N);

    fprintf('%d & %d\n',N,iter);
end