%% FRACTIONAL MODEL FOR OPTION PRICING
% FLMS model
clear;
clc;

% Parametri del modello
r = 1;
sigma = 0.35;
K = 80;
alpha = 1.9;
x0 = 0;
xN = 140;
T = 2;
L = xN - x0;
ic = @(x) -x.*(x-xN);    % Initial condition

% Parametri di discretizzazione
N = 100;
h = L/N;
x = x0:h:xN;

% Initial data for CALL option
for j = 2:N-1
    y0(j,1) = max([x(j)-K,0]);
end
% Condizione al bordo destra per CALL OPTION
b = @(t) max([y0-K*exp(-r*(t));0]);    

% Costruzione degli operatori
e = ones(N-1,1);
c1 = r + 0.5*(sigma^alpha)*sec(0.5*alpha*pi);
Tt = spdiags([-c1*e./(2*h),c1*e./(2*h)],[-1,1],N-1,N-1);
vr_left_T = zeros(N-1,1); vr_left_T(N-1) = c1/(2*h);
c2 = -0.5*(sigma^alpha)*sec(0.5*alpha*pi);
omega = cumprod([1, 1 - ((alpha+1) ./ (1:N))]); %Grunwald-Letnikov coeff.
Tleft = c2*toeplitz(omega(2:N),[omega(2),omega(1),zeros(1,N-3)])./(h^alpha);
vr_left = zeros(N-1,1); vr_left(N-1) = -c2*omega(1)/h^(alpha);
D = -r*speye(N-1);

J = Tt + Tleft + D; figure(1); eigJ=eig(J); plot(real(eigJ),imag(eigJ),'kx');
g = @(t) (vr_left + vr_left_T).*[ones(N-2,1);b(t)];


F = @(t,y) J*y + g(t);

[TT,Y] = ode23s(@(t,y) F(t,y),[0,T],y0);

figure(2)
pause()
for j=1:length(TT)
   plot(x,[0,Y(j,:),b(TT(j))],'bo-'); 
   title(sprintf('tau = %f, t = %f',T-TT(j),TT(j)));
   pause(1/25);
end

figure(3)
for j=1:length(TT)
   gg(j) = b(TT(j));
end
plot(TT,gg,'rx-')

