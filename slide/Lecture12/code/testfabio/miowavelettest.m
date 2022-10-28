%% TEST DELLA MIA ROUTINE PER PRECONDIZIONAMENTO WAVELET

clear
clc

%N=input('Introduci N:');
results = zeros(8,10);
j = 1;
for N=[100,500,1000,1500,2000,2500,5000,10000]
    x0 = 0; y0 = 0;
    x_max = 2; yN = 0;
    l = @(x) 0.6*ones(size(x));
    r = @(x) 0.5*ones(size(x));
    b = @(x) zeros(size(x));
    c = @(x) zeros(size(x));
    fracexp = 1.2;
    T = 2;
    t0 = 0;
    [L,F,dx] = FADE1Dgen(x0,x_max,N,y0,yN,l,r,b,c,fracexp);
    A = -L;
    xx = x0:dx:x_max;
    g21 = @(t) (-xx(2:end-1).*(xx(2:end-1) - 1))' + F;
    IVPfun = 0;
    y0 = exp(-(xx(2:end-1) - 1.2).^2./(2*(0.08)^2));
    
    tic;
    [x,~,~,results(j,1)] = bicgstab(A,y0',1e-6,N);
    results(j,2) = toc;
    lev = 6;
    dwt = 2;
    for type = 1:4
        [ L,U,Am,bm ] = genwaveprec( A,y0',lev,dwt,type );
        tic;
        [x2,~,~,results(j,2*(type)+1)] = bicgstab(Am,bm,1e-6,N,L,U);
        results(j,2*(type)+2) = toc;
    end
    j=j+1;
end