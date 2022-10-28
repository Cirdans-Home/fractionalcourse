function [x,flag,res,its] = fgmresrecycling(A,b,restart,maxit,tol)
%%FGMRES+GMRES(m) Flexible-GMRES algorithm with restart preconditioned with
% dbprec(v) algorithm (using the Bertaccini implementation).
% It is based on restarted GMRES with Krylov subspace of dim = im.

% Internal solution with susbspace recycling

machine_tolerance = eps;  % Spacing of double numbers.
verbose = 0;              % Use verbose mode, print all the iteration
n = size(b,1);
its = 0;                  % Number of iterations
x = sparse(n,1);
flag = 1;

% Inizialize some quantities in memory to make routine run faster:
res = -ones(maxit,1);
c = zeros(maxit,1);
s = zeros(maxit,1);
rs = zeros(maxit,1);
%%  MAIN LOOP OF THE ALGORITHM
while (its < maxit)
    vv(1:n,1) = b - A(x)  ;
    ro = norm(vv(1:n,1),2)  ;
    res(its+1) = ro;
    if (its  == 0)
        tol1=tol*ro  ;
    end  ;
    if (ro <= tol1 || its >= maxit)
        return
    end
    t = 1.0/ ro   ;
    vv(1:n,1) = vv(1:n,1) * t  ;
    %------------ initialize 1-st term  of rhs of hessenberg system.. ;
    rs(1) = ro  ;
    if (verbose)    % In verbose mode prints residual informations
        fprintf(1,' its %d  res %e \n',its,ro)
    end
    i = 0  ;
    while (i < restart  &&  (ro  >  tol1)  &&  its < maxit)
        i=i+1  ;
        its = its + 1  ;
        i1 = i + 1 ;
        %z = feval(precfun,PRE,vv(:,i));
        %[z,~] = gmres(A,vv(:,i),RESTART,TOL,MAXIT);
        z = dbprec(vv(:,1));
        Z(:,i) = z;
        %% Arnoldi with modified GS algorithm:
        vv(1:n,i1) = A(z);
        for j=1:i
            t = vv(1:n,j)'*vv(1:n,i1)  ;
            hh(j,i) = t  ;
            vv(1:n,i1) = vv(1:n,i1) - t*vv(1:n,j)  ;
        end  ;
        t = norm(vv(1:n,i1),2)  ;
        hh(i1,i) = t  ;
        if (t  ~= 0.0)
            t = 1.0 / t  ;
            vv(1:n,i1) = vv(1:n,i1)*t  ;
        end
        if (i ~= 1) %  perfrom previous transformations  on i-th column of h
            for k=2:i
                k1 = k-1  ;
                t = hh(k1,i)  ;
                hh(k1,i) = c(k1)*t + s(k1)*hh(k,i);
                hh(k,i) = -s(k1)*t + c(k1)*hh(k,i);
            end
        end
        gam = sqrt(hh(i,i)^2 + hh(i1,i)^2)  ;
        if (gam  == 0.0)
            gam = machine_tolerance  ;
        end 
        c(i) = hh(i,i)/gam; % determine plane rotation 
                            % and update rhs of ls pb;
        s(i) = hh(i1,i)/gam;
        rs(i1) = -s(i)*rs(i);
        rs(i) =  c(i)*rs(i);
        %
        %       determine res. norm. and test for convergence-  ;
        %
        hh(i,i) = c(i)*hh(i,i) + s(i)*hh(i1,i)  ;
        ro = abs(rs(i1));
        res(its+1) = ro;
        if (verbose)
            fprintf(1,' its %d  res %e \n',its,ro)
        end
    end
    %
    %       now compute solution. first solve upper triangular system.  ;
    %
    rs(i) = rs(i)/hh(i,i)  ;
    for  k=i-1:-1:1
        t=rs(k)  ;
        for j=k+1:i
            t = t-hh(k,j)*rs(j)  ;
        end  ;
        rs(k) = t/hh(k,k)  ;
    end  ;
    %
    %       done with back substitution..  ;
    %       now form linear combination to get solution  ;
    %         ;
    for j=1:i
        x = x +rs(j)*Z(1:n,j);
    end  ;
    if ((ro  <=  tol1) || (its >= maxit))
        res(res==-1)=[];
        flag = 0;
        return 
    end
end;  %% end while -- restart

end
