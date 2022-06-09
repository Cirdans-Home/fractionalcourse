%% Laplace Transform Inversion

clear; clc; close all;

alpha1 = 0.6;
alpha2 = 0.8;
c = 2.0;

a = @(t) alpha2 + (alpha1-alpha2).*exp(-c*t);
A = @(s) (alpha2*c + alpha1*s)./(s.*(c+s));

tlin = linspace(0,4,50);
t = [logspace(-6,-2,20),tlin(2:end)];
figure(1)
plot(t,a(t),'b-','LineWidth',2);

figure(2)
plot(t,t.^(t.*A(t)-1),'-',t,t.^(-t.*A(t)),'-','LineWidth',2)
legend('\Phi_\alpha(s)','\Psi_\alpha(s)')

%% Solution

y0 = 1;
lambda = 0.5;
Psi = @(s) s.^(-s.*A(s));
F = @(s) y0./(s.*(1 + lambda*Psi(s)));

L = -log(eps);
N = ceil(4*L/3/pi);
h = 2*pi/L + L/2/pi/N^2;
p = L^3/4/pi^2/N^2;
u = (0:N)*h;
f = zeros(size(t));
for n = 1:length(t)
    mu = p/t(n);
    z = mu*(u*1i + 1).^2;
    z1 = 2*mu*(1i-u);
    G = exp(z.*t(n)).*F(z).*z1;
    f(n) = ( imag(G(1))/2 + sum(imag(G(2:N+1))))*h/pi;
end


figure(3)
plot(t,f,'k-','LineWidth',2)

%% Comparisons

f_fun = @(t,y) -lambda*y; 
J_fun = @(t,y) -lambda;
t0 = 0; T = 4; h = 1e-2;
alpha = alpha1;
[t1, y1] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h);
alpha = alpha2;
[t2, y2] = fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h);

figure(3)
plot(t,f,'k-',...
    t1,y1,'r--',...
    t2,y2,'b--','LineWidth',2)
legend('Variable Order','\alpha = 0.6', '\alpha = 0.8');