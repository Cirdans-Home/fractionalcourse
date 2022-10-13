%% Circulant preconditioner for Fractional

clear; clc; close all;

%% Problem data
theta = 0.5; % Enforce symmetric case
alpha = 1.2;
w0 = @(x) 5*x.*(1-x);
%% Discretization data
for N  = 2.^(5:11)
    hN = 1/(N-1); x = 0:hN:1;
    dt = hN; t = 0:dt:1;
    %% Discretize
    G = glmatrix(N,alpha);
    Gr = G(2:N-1,2:N-1); Grt = Gr.';
    I = eye(N-2,N-2);
    % Left-hand side
    nu = hN^alpha/dt;
    A = nu*I - (theta*Gr + (1-theta)*Grt);
    % Right-hand side
    w = w0(x).';

    W = zeros(N-2,N-2);
    W(:,1) = w(2:N-1);

    figure(2)
    plot(x,w,'-','LineWidth',2);
    axis([0 1 0 2])
    %pause();
    for i=2:N
        w(2:N-1) = A\(nu*w(2:N-1));
        W(:,i) = w(2:N-1);
        plot(x,w,'-','LineWidth',2);
        axis([0 1 0 2]);
        pause(dt);
    end

    %% Circulant preconditioners

    % Tony Chan Preconditioner
    c = gencircprec( A(:,1),'tchan' );
    C = gallery('circul',ifft(c));
    evtchan = eig(A,C);
    % Strang Preconditioner
    c = gencircprec( A(:,1),'strang' );
    C = gallery('circul',ifft(c));
    evstrang = eig(A,C);
    % Raymond Chan Preconditioner
    c = gencircprec( A(:,1),'rchan');
    C = gallery('circul',ifft(c));
    evrchan = eig(A,C);
    % Jackson preconditioner
    c  = jacksonprec( A(:,1),2 );
    C = gallery('circul',ifft(c));
    evjackson = eig(A,C);

    ev = eig(A);
    figure(2)
    plot(ev,1,'k*',evtchan,2,'k*',evstrang,3,'k*',evrchan,4,'k*',...
        evjackson,5,'k*');
    yticks(1:5);
    yticklabels({'Original','T.Chan','Strang','R.Chan','Jackson'});
    ylim([0 6])

    %% Convergence
    w = w0(x).';
    % Tony Chan Preconditioner
    c = gencircprec( A(:,1),'tchan' );
    P = @(x) cprec(c,x);
    [X,FLAGchan,RELRESchan,ITERchan,RESVECchan] = ...
        pcg(A,(nu*w(2:N-1)),1e-9,N,P);
    % Strang Preconditioner
    c = gencircprec( A(:,1),'strang' );
    P = @(x) cprec(c,x);
    [X,FLAGstrang,RELRESstrang,ITERstrang,RESVECstrang] = ...
        pcg(A,(nu*w(2:N-1)),1e-9,N,P);
    % Raymond Chan Preconditioner
    c = gencircprec( A(:,1),'rchan');
    P = @(x) cprec(c,x);
    [X,FLAGrchan,RELRESrchan,ITERrchan,RESVECrchan] = ...
        pcg(A,(nu*w(2:N-1)),1e-9,N,P);
    % Jackson preconditioner
    c  = jacksonprec( A(:,1),2 );
    P = @(x) cprec(c,x);
    [X,FLAGjack,RELRESjack,ITERjack,RESVECjack] = ...
        pcg(A,(nu*w(2:N-1)),1e-9,N,P);
    % PCG
    [X,FLAG,RELRES,ITER,RESVEC] = ...
        pcg(A,(nu*w(2:N-1)),1e-9,N);

    %% PLOT
    figure(3)
    loglog(0:ITER,RESVEC,'o',...
        0:ITERjack,RESVECjack,'--',...
        0:ITERchan,RESVECchan,'--',...
        0:ITERrchan,RESVECrchan,'-.',...
        0:ITERstrang,RESVECstrang,'x--','LineWidth',2)
    legend('PCG','Jackson','T. Chan','R. Chan','Strang')
    xlabel('Iteration')
    ylabel('Residual')
    if N == 2^11
        set(gcf,'Color','none')
        cd ..
        export_fig(sprintf('N=%dalpha=%1.1f.png',N,alpha));
        cd Code/
    end

    fprintf(' %1.1f & $2^{%d}$ & %d & %d & %d & %d & %d \\\\\n', ...
        alpha,log2(N),ITER,ITERjack,ITERchan,ITERrchan,ITERstrang);
end