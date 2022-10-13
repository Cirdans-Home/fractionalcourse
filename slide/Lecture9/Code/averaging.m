%% THE NON SYMMETRIC CASE

clear; clc; close all;
dofov = false;

for alpha = [1.2,1.4,1.6,1.8]
    for N  = 2.^(5:11)
        w0 = @(x) 5*x.*(1-x);
        hN = 1/(N-1); x = 0:hN:1;
        dt = hN; t = 0:dt:1;
        dplus = @(x,t) gamma(3-alpha).*x.^alpha;
        dminus = @(x,t) gamma(3-alpha).*(2-x).^alpha;
        %% Discretize
        G = glmatrix(N,alpha);
        Gr = G; Grt = G.';
        Dplus = diag(dplus(x,0));
        Dminus = diag(dminus(x,0));
        I = eye(N,N);
        % Left-hand side
        nu = hN^alpha/dt;
        A = nu*I - (Dplus*Gr + Dminus*Grt);

        if dofov
            %% Computation
            F = fov(A);
            ev = eig(A);
            t = linspace(0,pi,100);

            %% PLOT
            figure(1)
            hold on
            plot(F,'r-','LineWidth',2);
            plot(real(ev),imag(ev),'b.')
            axis([-1 4 -1 1])
            hold off
            legend({'Field of Value','Eigenvalues'},'FontSize',14)
            hline(0,'k--')
            vline(0,'k--')
            set(gcf,'Color','none')
            export_fig('fovcase4.pdf')
        end

        w = w0(x).';
        % Wei-Sun preconditioner
        [ev,evt] = sunprec(N,alpha);
        c = nu + mean(dplus(x,0))*ev + mean(dminus(x,0))*evt;
        P = @(x) cprec(c,x);
        [X,FLAGsun,RELRESsun,ITERsun,RESVECsun] = ...
            gmres(A,(nu*w),[],1e-9,N,P);
        % GMRES
        [X,FLAG,RELRES,ITER,RESVEC] = ...
            gmres(A,(nu*w),[],1e-9,N);

        %% PLOT
        figure(3)
        loglog(0:prod(ITER),RESVEC,'o',...
            0:prod(ITERsun),RESVECsun,'--','LineWidth',2)
        legend('GMRES','Circulant')
        xlabel('Iteration')
        ylabel('Residual')
        if N == 2^11
            set(gcf,'Color','none')
            %cd ..
            export_fig(sprintf('gmresN=%dalpha=%1.1f.png',N,alpha));
            %cd Code/
        end
        if N == 2^5
            fprintf(' %1.1f & $2^{%d}$ & %d & %d \\\\\n', ...
                alpha,log2(N),prod(ITER),prod(ITERsun));
        else
            fprintf('& $2^{%d}$ & %d & %d \\\\\n', ...
                log2(N),prod(ITER),prod(ITERsun));
        end
    end
end