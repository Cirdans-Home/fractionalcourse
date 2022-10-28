%% TEST PER PROVARE A VEDERE DOVE STRANG NEL CASO FRAZIONARIO DIVENTA SINGOLARE

clear all
clc

x0 = 0;
xN = 1;
N = 501;
y0 = 0;
yN = 0;
l = @(x) 0.5*ones(size(x));
r = @(x) 0.6*ones(size(x));
b = r;
c = r;
alpha = 1.8;

[ A,F,h] = FADE1Dgen( x0,xN,N,y0,yN,l,r,b,c,alpha);
A = -A;

t = A(:,1);
c = gencircprec(A,'strang');
PA = cprec(c,A);
eigA = eig(A);
eigPA = eig(PA);
indneg = find(eigPA < 0);
indpos = find(eigPA > 0);
figure(1)
subplot(1,2,1); 
plot(real(eigA),imag(eigA),'kx'); 
title('Original Spectrum','Interpreter','LaTeX')
subplot(1,2,2); 
plot(real(eigPA(indneg)),imag(eigPA(indneg)),'rx',real(eigPA(indpos)),imag(eigPA(indpos)),'kx'); 
title('Strang Preconditioned Spectrum','Interpreter','LaTeX')
legend('Negative eigenvalues','Positive eigenvalues');

diagN = ceil(N/3);
P = spdiags( spdiags(A, -diagN:diagN), -diagN:diagN, N-1, N-1);
PPA = P\A;
eigPPA = eig(PPA);
figure(2)
subplot(1,2,1); 
semilogx(real(eigA),imag(eigA),'kx'); 
title('Original Spectrum','Interpreter','LaTeX')
subplot(1,2,2); 
semilogx(real(eigPPA),imag(eigPPA),'kx'); 
title('Band Toeplitz Preconditioned Spectrum','Interpreter','LaTeX')
