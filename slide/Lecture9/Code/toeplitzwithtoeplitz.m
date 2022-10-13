%% Toeplitz with Toeplitz

clear; clc; close all;

alpha = 0.8;
tol = 1e-6;
for k = 5:11
    N = 2^k;
    b = rand(N,1);
    T = kacmatrix(N,alpha);
    Tinv = invkacmatrix(N,alpha);
    [x,~,~,itnpc] = pcg(T,b,tol);
    [x,~,~,itpc] = pcg(T,b,tol,N,@(x) Tinv*x);
    fprintf('%d & %d & %d \\\\\n',N,itnpc,itpc);
end

%% FIGURE
N = 256;
alpha = 0.5;
T = kacmatrix(N,alpha);
Tinv = invkacmatrix(N,alpha);
I = eye(size(T));
figure(1)
spy(T)
set(gcf,'Color','none')
export_fig('T.pdf')

figure(2)
spy(Tinv)
set(gcf,'Color','none')
export_fig('Tinv.pdf')

figure(3)
spy(I)
set(gcf,'Color','none')
export_fig('I.pdf')

figure(4)
sv = svd(Tinv*T - I);
loglog(1:N,sv,'o--','LineWidth',2)
set(gcf,'Color','none')
export_fig('svd.pdf')