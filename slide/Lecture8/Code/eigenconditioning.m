%% Eigenvector matrix conditioning

iind = 1;
jind = 1;
for alpha = [1.2,1.4,1.6,1.8]
    for N  = 2.^(5:11)
        theta = 0.2; % Enforce non-symmetric case
        w0 = @(x) 5*x.*(1-x);
        hN = 1/(N-1); x = 0:hN:1;
        dt = hN; t = 0:dt:1;
        %% Discretize
        G = glmatrix(N,alpha);
        %Gr = G(2:N-1,2:N-1); Grt = Gr.';
        %I = eye(N-2,N-2);
        Gr = G; Grt = G.';
        I = eye(N,N);
        % Left-hand side
        nu = hN^alpha/dt;
        A = nu*I - (theta*Gr + (1-theta)*Grt);
        [V,D,W] = eig(A);
        data(iind,jind) = norm(V,2)*norm(W,2);
        jind = jind+1;
    end 
    iind = iind + 1;
    jind = 1;
end