%% TEST WITH MATRIX EQUATION - 2D Solve

clear; clc; close all;
iind = 1;
for alpha = [1.2,1.3,1.4,1.5,1.6,1.7,1.8]
    jind = 1;
    for N  = 2.^(5:11)
        w0 = @(x,y) 5*x.*(2-x).*y.*(2-y);
        hN = 2/(N-1); x = 0:hN:2; y = 0:hN:2;
        dt = hN; t = 0:dt:1; M = length(t);
        dplus = @(x,t) gamma(3-alpha).*x.^alpha;
        dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
        f = @(x,y,t) sin(pi*x).*sin(pi*y).*exp(-t);
        %% Discretize
        G = glmatrix(N,alpha);
        Gr = G; Grt = G.';
        Dplus = diag(dplus(x,0));
        Dminus = diag(dminus(x,0));
        I = eye(N,N);
        e = ones(N,1);
        % Left-hand side
        A = 0.5*I -dt*(Dplus*Gr + Dminus*Grt)/hN^alpha;
        B = (0.5*I -dt*(Dplus*Gr + Dminus*Grt)/hN^alpha).';
        [X,Y] = meshgrid(x,y);
        C = w0(X,Y) + dt*f(X,Y,t(1));
        C = -C';
        [U,S,V] = svd(C);
        C1 = U(:,1:2)*sqrt(S(1:2,1:2));
        C2 = (sqrt(S(1:2,1:2))*V(:,1:2).').';

        %% Solving
        m = 100;
        tol = 1e-9;
        [LA,UA] = lu(A);
        [LB,UB] = lu(B);
        [X1,X2,res]=kpik_sylv(A,LA,UA,B,LB,UB,C1,C2,m,tol);

        ITI(iind,jind) = length(res);

        if N == 2^5
            fprintf(' %1.1f & $2^{%d}$ & %d & %e \\\\\n', ...
                alpha,log2(N),length(res),res(end));
        else
            fprintf('& $2^{%d}$ & %d & %e \\\\\n', ...
                log2(N),length(res),res(end));
        end

        if N == 2^9 && alpha == 1.5
            SOL = X1*X2';

            figure(2)
            mesh(Y',X',SOL)
            xlabel('Y')
            ylabel('X')
        end
        jind = jind + 1;
    end
    iind = iind + 1;
end