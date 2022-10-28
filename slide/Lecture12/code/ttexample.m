%% TT-EXAMPLE

clear; clc; close all;


%% Discretization of the time-dependent part:
k = 3;          % Number of GBDF steps
kval = 11;       % Grid power
m = 2^kval;     % Number of time steps
[Alpha,Beta] = mab(k,m);
A = Alpha(:,2:m+1);
B = Beta(:,2:m+1);
t0 = 0;
tf = 1;
h = (tf-t0)/m;

tA = tt_matrix(full(A),1e-14);
tA = tt_reshape(tA,2*ones(kval,2));
tB = tt_eye(2,kval);


%% Discretization of the space operator:
n = 2^(kval+1); % Size of the space problem
alpha = 1.5; % Order of the fractional derivative
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

tL = tt_matrix(L,1e-14);
tL = tt_reshape(tL,2*ones(kval+1,2));
tM = tt_eye(2,kval+1);

tMat = tkron(tA,tM)-h*tkron(tB,tL);

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

sz = size(tMat);
ttb = tt_tensor(b(n+1:nk),eps);
ttb = tt_reshape(ttb,sz(:,1));

fprintf('\t %d & %d & %d & \\\\\n',k,m,n);
tx = amen_solve2(tMat,ttb,1e-6);