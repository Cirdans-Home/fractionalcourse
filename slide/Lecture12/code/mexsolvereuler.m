%% MATRIX EQUATION SOLVER

clear; clc; close all;

for kval = [5:12]
    n = 2^kval;
    alpha = 1.5;
    %% Problem
    w0 = @(x) 5*x.*(2-x);
    hN = 2/(n-1); x = 0:hN:2;
    dt = hN; t = 0:dt:1; m = length(t);
    dplus = @(x,t) gamma(3-alpha).*x.^alpha;
    dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
    f = @(x,t) -x.*(x-2).*exp(-t);
    %% Discretize Space Operator
    G = glmatrix(n,alpha);
    Gr = G; Grt = G.';
    Dplus = diag(dplus(x,0));
    Dminus = diag(dminus(x,0));
    I = eye(n,n);
    e = ones(n,1);
    % Left-hand side
    nu = hN^alpha/dt;
    M = speye(n,n);
    L = -dt*(Dplus*Gr + Dminus*Grt)/hN^alpha;
    %% Discretize Time operator
    A = spdiags([-e,e],-1:0,m,m);
    B = speye(m,m);

    [X,T] = meshgrid(x,t);
    C = dt*f(X,T);
    C(1,:) = w0(x) + C(1,:);
    C = -C';
    [U,S,V] = svd(C);
    C1 = U(:,1:2)*sqrt(S(1:2,1:2));
    C2 = (sqrt(S(1:2,1:2))*V(:,1:2).').';

    %% Matrix equation solver
    maxit = 100;
    tol = 1e-9;
    [LL,UL] = lu(L);
    [LA,UA] = lu(A);
    [X1,X2,res]=kpik_sylv(L,LL,UL,A,LA,UA,C1,C2,maxit,tol);

    fprintf('%d & %d & %d & %e \n',m,n,length(res),res(end));
end
%% Linear system solver
Mat = kron(A,M) + kron(B,L);
bvec = -C(:);
xvec = Mat\bvec;

%% Plotting
figure(1)
subplot(1,3,1)
mesh(t,x,X1*X2')
xlabel('t')
ylabel('x')
subplot(1,3,2)
mesh(t,x,reshape(xvec,n,m))
subplot(1,3,3)
mesh(t,x,log10(abs(X1*X2'-reshape(xvec,n,m))));