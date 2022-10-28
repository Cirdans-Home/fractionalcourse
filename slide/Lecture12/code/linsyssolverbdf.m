%% MATRIX EQUATION SOLVER

clear; clc; close all;

m = 10;     % Number of time steps
k = 2;          % Number of GBDF steps
n = 100; % Size of the space problem
alpha = 1.5; % Order of the fractional derivative

%% Discretization of the time-dependent part:
[Alpha,Beta] = mab(k,m);
A = Alpha(:,2:m+1);
B = Beta(:,2:m+1);
figure(1)
subplot(1,2,1); spy(A); title('$A_m$','Interpreter','latex');
subplot(1,2,2); spy(B); title('$B_m$','Interpreter','latex');
t0 = 0;
tf = 1;
h = (tf-t0)/m;

%% Discretization of the space operator:
u0f = @(x) 5*x.*(2-x);
hN = 2/(n-1); x = 0:hN:2;
u0 = u0f(x).';
dplus = @(x,t) gamma(3-alpha).*x.^alpha;
dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
f = @(x,t) -x.*(x-2).*exp(-t);
G = glmatrix(n,alpha);
Gr = G; Grt = G.';
Dplus = diag(dplus(x,0));
Dminus = diag(dminus(x,0));
I = eye(n,n);
e = ones(n,1);
% Matrix:
M = speye(n,n);
L = (Dplus*Gr + Dminus*Grt)/hN^alpha;
%% Building the right-hand side
nk=n*(m+1);
b=zeros(nk,1);
for j=1:m
    b(1+j*n:(j+1)*n)=f(x,t0+j*h);
end %for j
b(n+1:n*(m+1))=h*kron(Beta,speye(n))*b;
b(1:n)=u0;
Am = kron(Alpha(:,1),speye(n))-h*kron(Beta(:,1),L);
b(n+1:nk)=b(n+1:nk)-Am*u0;

% Computing the eigenvalues of the Omega-circulant approximation
omega = -1;
Omega = spdiags(omega.^([0,-(1:m-1)].'./m),0,m,m);
[ro,si] = rosi_bdf(k,1);

omegacol = zeros(m,1);
omegacol(1:fix(k/2)+1) = ro(fix(k/2)+1:-1:1);
omegacol(end-fix(k/2)+1:end) = omega*ro(fix(k/2)+2:end);
phi = fft(Omega'*omegacol);
psi = ones(m,1);

% Build whole preconditioner just for checking
F = dftmtx(m)/sqrt(m);
wA = real(Omega'*F'*diag(phi)*F*Omega);
P = kron(wA,I) -h*kron(B,L);
Mat = kron(A,I) -h*kron(B,L);

% Without preconditioner
[x,~,~,iter,resvec] = gmres(@(x) Mprod(A,B,L,h,x),b(n+1:nk),[],1e-9,100);

% With preconditioner
[xprec,~,~,iterprec,resvecprec] = ...
    gmres(@(x) Mprod(A,B,L,h,x),b(n+1:nk),[],1e-9,100,P);%...
% @(v) Psolve(omega,phi,psi,h,L,v));

fprintf('%d & %d & %d & %d & %1.2e & %d & %1.2e \\\\\n',...
    k,m,n,prod(iter),resvec(end),prod(iterprec),resvecprec(end));

% Did we survive?
%         figure(2)
%         semilogy(0:prod(iter),resvec,'k-',...
%             0:prod(iterprec),resvecprec,'o--','LineWidth',2);

