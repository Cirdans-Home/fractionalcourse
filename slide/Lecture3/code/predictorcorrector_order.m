%% Convergence Order for FDE with Explicit Trapezoidal Rule

clear; clc; close all;

t0      = 0;
T       = 1;
alpha   = 0.25;

f_fun = @(t,y) 40320/gamma(9-alpha)*t.^(8-alpha) ...
    -3*gamma(5+alpha/2)/gamma(5-alpha/2)*t.^(4-alpha/2) ...
    + 9/4*gamma(alpha+1) + (3/2*t.^(alpha/2) - t^4)^3 - (y).^(3/2);
J_fun = @(t,y) (3/2)*t.^(3/2-1);

ye = @(t) t.^8 - 3*t.^(4+alpha/2) + 9/4*t.^alpha;
y0 = 0;

h = 1e-2;
[T, Y] = fde_pi1_ex(alpha,f_fun,t0,T,y0,h);

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
    T       = 1;
    [T, Y] = fde_pi1_ex(alpha,f_fun,t0,T,y0,h);
    err(i) = norm(Y - ye(T),'inf');
    %err(i) = abs(Y(end)-ye(T(end)));
end

%% Check convergence order
for i=1:length(hval)
    if i == 1
        fprintf('%1.2f & %1.2e & %1.2e & \\\\\n',alpha,hval(i),err(i));
    else
        q = log2(err(i-1)/err(i));
        fprintf('     & %1.2e & %1.2e & %1.2f \\\\\n',hval(i),err(i),q);
    end
end

%% Predictor-Corrector Method
hval = 2.^(-(1:14));
for i=1:length(hval)
    h       = hval(i);
    t0      = 0;
    T       = 1;
    [T, Y] = fde_pi12_pc(alpha,f_fun,t0,T,y0,h,[],1);
    err(i) = norm(Y - ye(T),'inf');
    %err(i) = abs(Y(end)-ye(T(end)));
end

%% Check convergence order
for i=1:length(hval)
    if i == 1
        fprintf('%1.2f & %1.2e & %1.2e & \\\\\n',alpha,hval(i),err(i));
    else
        q = log2(err(i-1)/err(i));
        fprintf('     & %1.2e & %1.2e & %1.2f \\\\\n',hval(i),err(i),q);
    end
end
