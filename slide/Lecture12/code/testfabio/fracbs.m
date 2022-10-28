%% ESEMPIO FMLS MODEL
clear;
clc;
%% Parametri del problema
alpha = 1.9;
sigma = 0.35;
r = 1;
K = 80;
x0 = 0;
xN = 140;
T = 1/2;
%% Parametri di discretizzazione
N = 200;
A = (r + 0.5*(sigma^alpha)*sec(alpha*pi*0.5));
b = @(x) -0.5*(sigma^alpha)*sec(alpha*pi*0.5)*ones(size(x));
c = @(x) zeros(size(x));
f = @(x) ones(size(x));
h = @(x) zeros(size(x));
d = -r;
%% Discretizzazione
dx = (xN-x0)/(N-1);
x = (x0:dx:xN)';
% Termine di trasporto 
e = ones(N,1);
Tr = spdiags([-A*e./(2*dx),A*e./(2*dx)],[-1,1],N,N); % Transport Term
% Derivata frazionaria sinistra e destra
B = spdiags(b(x),0,N,N);
C = spdiags(c(x),0,N,N);
omega = cumprod([1, 1 - ((alpha+1) ./ (1:(N)))]);
Tleft = -toeplitz(omega(2:N+1),[omega(2),omega(1),zeros(1,N-2)])./(dx^alpha);
fvec = f(x');
Tright = Tleft';
hvec = h(x');
% Termine di reazione
R = d*speye(N);
% Intervallo temporale
m = 100;
tt = linspace(0,T,m);
dt = tt(2)-tt(1);
% Parti della matrice per CN
P = zeros(N,N);
Q = zeros(N,N);
L = zeros(N,N);
P(2:N-1,2:N-1) = (dt/2)*B(2:N-1,2:N-1)*Tleft(2:N-1,2:N-1).*fvec(2:N-1);
Q(2:N-1,2:N-1) = (dt/2)*C(2:N-1,2:N-1)*Tright(2:N-1,2:N-1).*hvec(2:N-1);
L(2:N-1,2:N-1) = (dt/2)*Tr(2:N-1,2:N-1)+eye(N-2)+dt*R(2:N-1,2:N-1);
I = eye(N);
LEFT = I -P-Q;
RIGHT = L+P+Q;
V(:,1) = x';
z = zeros(N,1);

figure(1)
plot(x,V(:,1),'o-');
pause()
for j=2:m
    z(end) = max([x(end)-K*exp(-r*(T-tt(j))),0])
    V(:,j) = LEFT\(RIGHT*V(:,j-1)+z);
    plot(x,V(:,j),'o-');
    pause(1/25);
end

