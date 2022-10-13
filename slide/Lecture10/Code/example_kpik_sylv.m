nh=30;
T=diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I=speye(nh);
A=-(kron(T,I)+kron(I,T));
n=nh^2;
B=-spdiags(rand(n,1),0,n,n);
[LA,UA]=lu(A);
[LB,UB]=lu(B);
rhs1=randn(n,2);
rhs2=randn(n,2);
m=100;
tol=1e-9;
[X1,X2,r]=kpik_sylv(A,LA,UA,B,LB,UB,rhs1,rhs2,m,tol);

fprintf('final true absolute residual norm: \n')
disp(norm(A*X1*X2'+X1*X2'*B'+rhs1*rhs2'))    %this matrix should never be formed for n large 
