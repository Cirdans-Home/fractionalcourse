function [y,res,err,time] = MatFunLanczos(A,v,alpha,tol,maxit)
%% MatFunLanczos 
% This Algorithm Uses the Polynomial Lanczos algorithm
% for computing the matrix function times vector product A^(\alpha/2)v
% Input:    A Square Matrix (Symmetric and Positive Definite)
%           v Vector
%           alpha Fractional Power
%           tolerance for the residual of the Krylvo subspace
%           maxit maximum number of allowed iteration
% Output:   y Result of the computation
%           res residual of the Krylov subspace
%           time overall time of soluton

if size(A,1) < 20000
    %yexact = mpower(full(A),alpha/2)*v;
    yexact = sqrtm(full(A))\v;
else
    [V,D] = eig(A,'vector');
    yexact = (V * (spdiags(D.^(alpha/2),0,size(A,1),size(A,1)) * (V\v)));
end

tic;
% Allocation for the Lanczos Algorithm
n = size(A,1);
V = zeros(n,maxit);      % Contains the Basis of the Krylov Subspace
T = sparse(maxit,maxit); % Tridiagonal Projection
% Allocation of auxiliar stuff
res = zeros(maxit,1);
stopping = zeros(maxit,1);

beta = norm(v,2);
V(:,1) = v/beta;
for j=1:maxit
    z = A*V(:,j);
    T(j,j) = V(:,j)'*z;
    if j > 1
        vhat = z - T(j,j)*V(:,j) - T(j-1,j)*V(:,j-1);
    else
        vhat = z - T(j,j)*V(:,j);
    end
    T(j+1,j) = norm(vhat,2);
    T(j,j+1) = T(j+1,j);
    V(:,j+1) = vhat/T(j+1,j);
    T1full = full(T(1:j,1:j));
    ym = mpower(T1full,alpha/2)*eye(j,1);
    nrmy = norm(ym,2);
    % Error Check
    if nrmy > eps
        stopping(j) = beta*abs(ym(end))/abs(ym(1));
        if stopping(j) > .5
            res(j+1) = stopping(j);
        else
            res(j+1) = abs(ym(end));
        end
    else
        res(j+1)=1;
    end
    if res(j+1) < tol
        break
    end
end
y = beta*V(:,1:j)*ym;
time = toc;
res(res == 0) = [];
fprintf('Iter %d Error %1.2e Residual %1.2e Time %1.2e\n',j,norm(yexact-y),res(end),time);
err = norm(yexact-y)/norm(yexact);
fprintf(' %d & %1.2e & %1.2e & %1.2e \n',j,err,res(end),time);
end