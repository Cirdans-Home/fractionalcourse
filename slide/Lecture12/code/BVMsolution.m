%% BVM SOLUTION

clear; clc; close all;

for k = 2:8
    for kval = [5:12]
        m = 2^kval;     % Number of time steps
        %k = 2;          % Number of GBDF steps
        n = 2^(kval+1); % Size of the space problem
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
        %% Solution as a linear system
        if nk <= 10000
            Mat = kron(A,M) - h*kron(B,L);
            rhs = b(n+1:nk);
            u = Mat\rhs;

            U = [u0,reshape(u,n,m)];
            t = t0:h:tf;
            [TT,UU] = ode15s(@(t,y) L*y + f(x.',t),t,u0);
            figure(2)
            subplot(1,3,1)
            mesh(t,x,U);
            xlabel('t');
            ylabel('x');
            title('GBDF(6,100) on 100')
            subplot(1,3,2)
            mesh(t,x,reshape(UU,m+1,n).')
            xlabel('t');
            ylabel('x');
            title('ode15s')
            subplot(1,3,3)
            mesh(t,x,log10(abs(U-reshape(UU,m+1,n).')))
            xlabel('t');
            ylabel('x');
            title('Error')
        end

        %% Matrix equation solution
        C = reshape(b(n+1:nk),n,m);
        [U,S,V] = svds(-C,1);
        C1 = U*sqrt(S);
        C2 = (sqrt(S)*V.').';

        maxit = 100;
        tol = 1e-9;
        [LL,UL] = lu(-h*L);
        [LA,UA] = lu(A);
        [X1,X2,res]=kpik_sylv(-h*L,LL,UL,A,LA,UA,C1,C2,maxit,tol);

        fprintf('%d & %d & %d & %d & %1.2e\\\\\n',...
            k,m,n,length(res),res(end));

    end
end

if nk <= 20000
    U = [u0,reshape(X1*X2',n,m)];
    figure(4)
    t = t0:h:tf;
    mesh(t,x,U);
    xlabel('t');
    ylabel('x');
    title('Matrix equation solution')
end