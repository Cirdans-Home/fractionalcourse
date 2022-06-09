%% ABSOLUTE STABILITY

clear; clc; close all;

alpha = 0.5;
y0 = 1;
t0 = 0;
T = 1;
lambda = -5;

f = @(t,y) lambda*y;
Jf = @(t,y) lambda;
ye = @(t) ml(lambda*t.^alpha,alpha)*y0;

h = 1e-3;
[T, Y] = fde_pi1_ex(alpha,f,t0,T,y0,h);

figure(1)
plot(T,Y,'o',T,ye(T),'-','LineWidth',2);
xlabel('t')
ylabel('y(t)');
legend({'Numerical Approximation','True Solution'},...
    'Location','northeast','FontSize',13);

% hval = 2.^(-(1:14));
% for i=1:length(hval)
%     h       = hval(i);
%     t0      = 0;
%     T       = 1;
%     [T, Y] = fde_pi1_ex(alpha,f_fun,t0,T,y0,h);
%     err(i) = norm(Y - ye(T),'inf');
%     %err(i) = abs(Y(end)-ye(T(end)));
% end