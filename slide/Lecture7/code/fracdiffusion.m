%% Fractional Diffusion Equation

N = 100;
alpha = 1.5;

g = gl(N,alpha);
bound = @(N,alpha) abs(-alpha/gamma(1-alpha).*N.^(-alpha-1));

v1 = abs(gl(N,1.2));
v2 = abs(gl(N,1.5));
v3 = abs(gl(N,1.8));
figure(1)
loglog(0:10:N,v1(1:10:end),'o',80:N,bound(80:N,1.2),'k--',...
    0:10:N,v2(1:10:end),'x',80:N,bound(80:N,1.5),'k--',......
    0:10:N,v3(1:10:end),'v',80:N,bound(80:N,1.8),'k--','LineWidth',2)
legend('\alpha=1.1','bound','\alpha=1.5','bound','\alpha=1.8','bound')
axis tight

%% Problem data
theta = 0.9;
alpha = 1.8;
w0 = @(x) 5*x.*(1-x);
%% Discretization data
hN = 1/(N-1); x = 0:hN:1;
dt = hN; t = 0:dt:1;
%% Discretize
G = glmatrix(N,alpha); Gt = glmatrix(N,alpha).';
I = eye(N,N);
% apply B.C.
G(1,:) = -I(1,:); G(N,:) = -I(N,:);
Gt(1,:) = -I(1,:); Gt(N,:) = -I(N,:);
% Left-hand side
A = I - dt/hN^alpha*(theta*G + (1-theta)*Gt);
% Right-hand side
w = w0(x).'; 

W = zeros(N,N);
W(:,1) = w;

figure(2)
plot(x,w,'-','LineWidth',2);
axis([0 1 0 2])
pause();
for i=2:N
    w = A\w;
    W(:,i) = w;
    plot(x,w,'-','LineWidth',2);
    axis([0 1 0 2]);
    pause(dt);
end

figure(3)
[X,T] = meshgrid(x,t);
mesh(X,T,W)
ylabel('x')
xlabel('t')
zlabel('w(x,t)')
%set(gcf,'Color','none')
%export_fig(sprintf('sol_alpha%d_skew.eps',alpha*10));

%% Distribution
N = 100;
alpha = 1.5;
G = glmatrix(N,alpha);
s = @(t) -exp(-1i*t).*(1 + exp(1i*(t+pi))).^alpha; 
sv = svd(G);
th = linspace(0,2*pi,N);

figure(4)
plot(th,sv,'o',th,sort(abs(s(th)),'descend'),'-','LineWidth',2);
axis tight

%% Smallest SV

alpha = 1.5;
sv = [];
for N = [100,500,1000,5000]
    G = glmatrix(N,alpha);
    sv = [sv,svds(G,1,'smallest')];
end
