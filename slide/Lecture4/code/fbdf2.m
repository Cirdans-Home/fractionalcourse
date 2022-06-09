%% Convergence Order for FDE with Explicit Trapezoidal Rule

clear; clc; close all;

t0      = 0;
T       = 2;
alpha   = 0.5;

f_fun = @(t,y) -2*y;
J_fun = @(t,y) -2;

y0 = 1;
ye = @(t) ml(-2*t.^alpha,alpha)*y0;


h = 1e-2;
[T, Y] = flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],3,eps,100);

figure(1)
plot(T,Y,'o',T,ye(T),'-','LineWidth',2);
xlabel('t')
ylabel('y(t)');
legend({'Numerical Approximation','True Solution'},...
    'Location','southeast','FontSize',13);

%% Fractional Rectangular Rule
hval = 2.^(-(1:14));
for i=1:length(hval)
    h       = hval(i);
    t0      = 0;
    T       = 2;
    [T, Y] = flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],3,eps,100);
    err(i) = norm(Y - ye(T),'inf');
    %err(i) = abs(Y(end)-ye(T(end)));
end

%% Check convergence order
for i=1:length(hval)
    if i == 1
        fprintf('%1.2f & $2^{%d}$ & %1.2e & \\\\\n',alpha,round(log2(hval(i))),err(i));
    else
        q = log2(err(i-1)/err(i));
        fprintf('     & $2^{%d}$ & %1.2e & %1.2f \\\\\n',round(log2(hval(i))),err(i),q);
    end
end