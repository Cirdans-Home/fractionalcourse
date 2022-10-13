%% Asymptotic distribution for Fractional Differential Equation

clear; close; clc;

for N = [10:10:80].^2

    alpha = 1.5;
    hN = 1/(N-1);  x = 0:hN:1; dt = hN;
    dplus=@(x)gamma(3-alpha).*x.^alpha;
    dminus=@(x)gamma(3-alpha).*(1-x).^alpha;
    G = glmatrix(N,alpha); % Discretize
    Gr = G; Grt = G.'; I = eye(N,N);
    Dplus = diag(dplus(x));
    Dminus = diag(dminus(x));
    nu = hN^alpha/dt;
    A = nu*I -(Dplus*Gr + Dminus*Grt);
    f = @(theta) -exp(-1i*theta).*...
        (1-exp(1i*theta)).^alpha;
    xsq = linspace(0,1,sqrt(N));
    tsq = linspace(-pi,pi,sqrt(N));
    [X,THETA] = meshgrid(xsq,tsq);
    h = @(x,theta) nu + (dplus(x).*f(theta) ...
        + dminus(x).*f(-theta));
    sv = svd(full(A));
    figure(1)
    index = linspace(0,1,N);
    VAL = abs(h(X,THETA));
    plot(index,sort(VAL(:),'descend'),'r-',....
        index,sort(sv,'descend'),'b:','LineWidth',2)
    axis tight
    set(gcf,'Color','none')
    export_fig(sprintf("distrN%d.pdf",N));
end
