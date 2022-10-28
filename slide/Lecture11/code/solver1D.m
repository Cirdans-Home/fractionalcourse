%% SOLVING 1D Problem in HODLR

addpath('/home/cirdan/matlabauxiliarycodes/hm-toolbox')
clear; clc; close all;

index = 1;
kvec = 7:14;
for k=kvec
    %% Discretization
    N = 2^k; hN = 1/(N-1);  x = 0:hN:1; dt = hN;
    %% Coefficients
    alpha = 1.5;
    dplus=@(x)gamma(3-alpha).*x.^alpha;
    dminus=@(x)gamma(3-alpha).*(1-x).^alpha;
    w = @(x) 5*x.*(1-x);
    %% HODLR building
    tol = 1e-9;
    tic;
    G = glhodlrmatrix(N,alpha,tol);
    Dplus = hodlr('diagonal',dplus(x));
    Dminus = hodlr('diagonal',dminus(x));
    I = hodlr('eye', N);
    nu = hN^alpha/dt;
    A = nu*I -(Dplus*G + Dminus*G');
    buildtime(index) = toc;
    %% Solving
    [L,U] = lu(A);
    flu = @() lu(A);
    timelu(index) = timeit(flu,2);
    w = w(x).';
    plot(x,w); pause(0.05);
    solvetime(index) = 0;
    for i=1:N
        tic;
        w = U\(L\(nu*w));
        solvetime(index) = solvetime(index) + toc;
        plot(x,w);
        axis([0 1 0 2]);
        pause(dt)
    end
    solvetime(index) = solvetime(index)/N;

    fprintf('$2^{%d}$ & %e & %e & %e \\\\\n',...
        log2(N),buildtime(index),timelu(index),solvetime(index));
    index = index + 1;
end

%% save data
save('1Dscaling.mat','buildtime',"timelu","solvetime","kvec");

%% Plot results

figure(2)
loglog(2.^kvec,buildtime,'s-',2.^kvec,timelu,'o-',...
    2.^kvec,solvetime,'x-','LineWidth',2);
xlabel('$N$','Interpreter','latex');
ylabel('Time (second)');
legend({'Build time (s)','LU time (s)','Avg. solve time (s)'},...
    "Location","best");

figure(3)
semilogy(2.^kvec,buildtime./buildtime(1),'s-',...
    2.^kvec,timelu./timelu(1),'o-',...
    2.^kvec,solvetime./solvetime(1),'x-','LineWidth',2);
xlabel('$N$','Interpreter','latex');
ylabel('Speedup');
legend({'Build time','LU time','Avg. solve time'},...
    "Location","southeast");
axis tight
