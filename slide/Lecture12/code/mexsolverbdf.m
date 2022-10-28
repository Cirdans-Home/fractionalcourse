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
    r = zeros(m,1);
    c = zeros(m,1);
    r(1:7) = [147/60,-6,15/2,-20/3,15/4,-6/5,1/6];
    c(1) = 147/60;
    A = toeplitz(r,c);
    % Fix boundary conditions
    A(1,1) = 1;
    A(2,1) = -2; A(2,2) = 3/2;
    A(3,1) = 3/2; A(3,2) = -3; A(3,3) = 11/6;
    A(4,1) = -4/3; A(4,2) = 3; A(4,3) = -4; A(4,4) = 25/12;
    A(5,1) = 5/4; A(5,2) = -10/3; A(5,3) = 5; A(5,4) = -5; A(5,5) = 137/60;
    A = sparse(A);
    B = speye(m,m);

    [X,T] = meshgrid(x,t);
    C = dt*f(X,T);
    C(1:6,:) = [1,-1/2,1/3,-1/4,1/5,-1/6]'*w0(x) + C(1,:);
    C = -C';
    [U,S,V] = svd(C);
    C1 = U(:,1:2)*sqrt(S(1:2,1:2));
    C2 = (sqrt(S(1:2,1:2))*V(:,1:2).').';

    %% Matrix equation solver
    maxit = 600;
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