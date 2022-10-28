%% ILLUSTRAZIONE DEL LEMMA SULLO SHORT-MEMORY PRINCIPLE

clear;
clc;

addpath('../generatori_FDE');

N = 1000;
diagNdiv = 5;
x0 = 0; y0 = 0;
x_max = 2; yN = 0;
l = @(x) 0.6*ones(size(x));
r = @(x) 0.5*ones(size(x));
b = @(x) zeros(size(x));
c = @(x) zeros(size(x));
fracexp = 1.2;
T = 2;
t0 = 0;
[L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);

diagN = ceil(N/diagNdiv);
FPR = spdiags( spdiags(L, -diagN:diagN), -diagN:diagN, N-1, N-1);

PM = FPR\L;
eigPM = eig(PM);
svdPM = svd(PM);

figure(1)
%subplot(1,2,1); 
semilogx(real(eigPM),imag(eigPM),'k+')
axis square
%subplot(1,2,2); loglog(svdPM,'ko-');