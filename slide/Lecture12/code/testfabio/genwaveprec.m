function [ L,U,Am,bm ] = genwaveprec( A,rhs,lev,dwt,type )
%GENWAVEPREC this routine builds a Wavelet--based preconditioner for the
%matrix A.
%   lev := Number of levels for the DWT transform
%   dwt := 
%   type := 1 Simple DWT 
%           2 DWT with permutation type 0
%           3 DWT with permutation type 2
%           4 DWT with permutation type 3

%% Compute information about the matrix
n = size(A,1);
nnzA = nnz(A);
%% General setup
mu=min(max(floor(2*nnzA/n),0), (dwt-1)*2^(lev-1));
mu=max(mu,5);
A1=fwts(A,dwt,1,1, lev); 
b1=fwt(rhs,dwt,1,0, lev);

switch type
    case 1
        r=n:-1:1;
        B=banda(A1(r,r),mu,mu);
        Am = A1(r,r);
        bm = b1;
        [L,U]=lu(B);
    case 2
        Ap = perm2(A1,1,1,lev);
        bm = perm2(b1,1,0,lev);
        r=n:-1:1;
        Am=Ap(r,r);
        B=banda(Am,mu,mu,mu,mu);
        [L,U]=lu(B);
    case 3
        Am = perm0(A1,1,1,lev,0);
        bm = perm0(b1,1,0,lev,0);
        r=n:-1:1;
        Am=Am(r,r);
        B=banda(Am,mu,mu,mu,mu);
        [L,U]=lu(B);
    case 4
        Am = perm0(A1,1,1,lev,1);
        bm = perm0(b1,1,0,lev,1);
        r=n:-1:1;
        Am=Am(r,r);
        B=banda(Am,mu,mu,mu,mu);
        [L,U]=lu(B);
    otherwise
        error('Type unrecognised')
end
end

function B=banda(A,mu,mu2,gam1,gam2,tol)
%B = banda(A)  % extract a diagonal banded arrow matrix out of A
%B = banda(A,mu)
%B = banda(A,mu_L,mu_U,gam1,gam2)    % Bordering gam1 + gam2
%B = banda(A,mu_L,mu_U,gam1,gam2,tol)  % tol is for optional thresh
%
% e.g.  Ex_1=rand(64); B=banda(Ex_1,2,5,4,2); spy(B)
%       hold on, spy(diag(diag(Ex_1)),'ro'); hold off
%______________________________________
%ALSO  banda => banded arrow extraction
%      bandg => banded matrix with a gap
if nargin<1,help banda; return; end
if nargin<2,mu=1; end
if nargin<3,mu2=mu; end
if nargin<4,gam1=0; gam2=0; end
if nargin<5,gam2=0; end
[m,n]=size(A);
B=tril(tril(A(1:m-gam1,1:n-gam2),mu2)',mu)';
B=[B; A(m-gam1+1:m,1:n-gam2)];
B=[B  A(:,n-gam2+1:n)];
if nargin>5
    B=sparse(B);
    B=B.*(abs(B)+diag(1:length(B))*tol>tol);
end
end

function M=fwts(A,m_in,d1,d2,t,formy,iout) % formy 0/normal 1/NS
%% fwts.m --- for Daubechies' order m FWT transforms (sparse version) vs fwt.m
%% Usage: [1] y =fwts(x)           FWT for vector x
%%        [2] y =fwts(x,m)         FWT for vector x (m order)
%%        [3] B =fwts(A)           FWT for matrix A
%%        [4] B =fwts(A,m)         FWT for matrix A (m order)
%%        [5] B =fwts(A,m,d1,d2)   FWT for matrix A with 1-sided options
%%        [6] B =fwts(A,m,d1,d2,t) for t-levels of FWT for matrix A
%%                                [dj=1 on,=0 off] with 1-sided options
%%     [7] B =fwts(A,m,d1,d2,t,s,iout) for t-levels and s=0 (normal) 1 (NS)
%%e.g. A=cz(128,2);                    % iout = 0 (default) or =1 (printout)
%      B = fwts(A,4, 1, 1, 3); % for 3-levels of Daub_4 (both sides) A=A*W^T
%      B2= fwts(A,6, 1, 0, 2); % for 2-levels of Daub_6 (left only)  A=W*A
%___________________________________________ Internal Control
if nargin <1, eval(['help ' mfilename]),
    disp('! At least one parameter needed');  return, end
[n1,n2]=size(A); n=max(n1,n2); M=A;
if nargin <2, m_in=4; end %------- Set Daub4 (as default)
if nargin <3 %-------------------- Set d1 / d2
    d1=1; d2=1;
    if n1==1, d1=0; d2=1; end
    if n2==1, d2=0; d1=1; end
end
if nargin <6,formy=0; end, if nargin <7,iout=0; end
t_max=round( log(n)/log(2) );          % n = 2^t_max
if nargin <5 || isempty(t) %-------------------- Set t
    t=t_max+1 - ceil( log(m_in)/log(2) );  % Max active levels up to m
end
if iout==1
    Pat='A=WAW'''; if d2==0; Pat='A=WA'; end, if d1==0; Pat='A=AW'''; end
    fprintf(' fwts for %d levels of Daub%d to A of %d x %d [%s,%d]\n\t\tDone levs:',...
        t,m_in,n1,n2, Pat,formy);
end
N = 2*floor(n/2);

for iter = 1:t
    k=floor( n/2^iter ) * 2;
    T=speye(n,n); p=1:n;
    % NS form ----------------------------------------
    if formy==1
        T(N-k+1:N,N-k+1:N)=w_matrix(m_in,k);
        if (d1==1),M=T*M; end
        if (d2==1),M=M*T';end
        p(N+1-k/2:N)=N-k+1:2:N;     % Odd ...... | A B |
        p(N-k+1:N-k/2)=N-k+2:2:N;   % Even ..... | C T |
        % Normal form ----------------------------------------
    else
        T(1:k,1:k)=w_matrix(m_in,k);
        if (d1==1),M=T*M; end
        if (d2==1),M=M*T';end
        p(1:k/2)=1:2:k;     % Odd ...... | T B |
        p(k/2+1:k)=2:2:k;   % Even ..... | C A |
    end
    
    if (d1==1 &&  d2==1),M=M(p,p);end
    if (d1==1 &&  d2==0),M=M(p,:);end
    if (d1==0 &&  d2==1),M=M(:,p);end
    k1=k/2; if k1*2 ~= k, iter=1; end
    if iout==1,fprintf(' %d',iter),end
end % iter - t levels;
if iout==1,disp(' '),end
return
end

function W = w_matrix(m,n)
%% ________________________________________________________
%% function W=w_matrix(m, n) - matrix size n with m filters
p = find([2 4 6 8 16] == m);
if nargin<1, help w_matrix; return; end
if m<2, W=speye(n); return; end
%  Note these DIAG repeating patterns:   (0,c_{m-1}), (c_{m-2} d_{m-1})
%   (c_{m-4} d_{m-3}) (c_{m-6} d_{m-5}) ... (c_{0} d_1) and (c_{m-1} 0)
if m==2
    c1=sqrt(2)/2; c0=sqrt(2)/2;
    d1=zeros(n,1); d1(1:2:n)=-c0*ones(n/2,1);
    d3=c0*ones(n,1);
    d5=zeros(n,1); d5(2:2:n)= c0*ones(n/2,1);
    W=spdiags([d1 d3 d5],-1:1,n,n);
elseif m==4      %====================================================
    c3=sqrt(2)*(1+sqrt(3))/8; c2=sqrt(2)*(3+sqrt(3))/8;
    c1=sqrt(2)*(3-sqrt(3))/8; c0=sqrt(2)*(1-sqrt(3))/8;
    
    d1=zeros(n,1); d1(1:2:n)=c0*ones(n/2,1);
    d2=zeros(n,1); d2(1:2:n)=c3*ones(n/2,1); d2(2:2:n)=-c1*ones(n/2,1);
    d3=c2*ones(n,1);
    d4=zeros(n,1); d4(1:2:n)=c1*ones(n/2,1); d4(2:2:n)=-c3*ones(n/2,1);
    d5=zeros(n,1); d5(2:2:n)=c0*ones(n/2,1);
    
    W=spdiags([d1 d2 d3 d4 d5],-1:3,n,n);
    W(n-1:n,1:2)=[c1 c0; c2 -c3];
elseif m==6      %====================================================
    c0= 3.326705529500826d-01; c1= 8.068915093110927d-01;
    c2= 4.598775021184915d-01; c3=-1.350110200102546d-01;
    c4=-8.544127388202666d-02; c5= 3.522629188570955d-02;
    
    d1=zeros(n,1); d1(1:2:n)=c5*ones(n/2,1);
    d2=zeros(n,1); d2(1:2:n)=c0*ones(n/2,1); d2(2:2:n)=-c4*ones(n/2,1);
    d3=zeros(n,1); d3(1:2:n)=c3*ones(n/2,1); d3(2:2:n)= c1*ones(n/2,1);
    d4=zeros(n,1); d4(1:2:n)=c2*ones(n/2,1); d4(2:2:n)=-c2*ones(n/2,1);
    d5=zeros(n,1); d5(1:2:n)=c1*ones(n/2,1); d5(2:2:n)= c3*ones(n/2,1);
    d6=zeros(n,1); d6(1:2:n)=c4*ones(n/2,1); d6(2:2:n)=-c0*ones(n/2,1);
    d7=zeros(n,1);                           d7(2:2:n)= c5*ones(n/2,1);
    
    W=spdiags([d1 d2 d3 d4 d5 d6 d7],-1:5,n,n);
    W(n-3:n-2,1:2)=[c4 c5; c1 -c0];
    W(n-1:n,3:4)=[c4 c5; c1 -c0];
    W(n-1:n,1:2)=[c2 c3; c3 -c2];
elseif m==8      %====================================================
    c0=-7.5765714789357d-2; c1=-2.9635527645960d-2;
    c2= 4.9761866763256d-1; c3= 8.0373875180539d-1;
    c4= 2.9785779560561d-1; c5=-9.9219543576956d-2;
    c6=-1.2603967262264d-2; c7= 3.2223100604078d-2;
    
    d1=zeros(n,1); d1(1:2:n)=c7*ones(n/2,1);
    d2=zeros(n,1); d2(1:2:n)=c0*ones(n/2,1); d2(2:2:n)=-c6*ones(n/2,1);
    d3=zeros(n,1); d3(1:2:n)=c5*ones(n/2,1); d3(2:2:n)= c1*ones(n/2,1);
    d4=zeros(n,1); d4(1:2:n)=c2*ones(n/2,1); d4(2:2:n)=-c4*ones(n/2,1);
    d5=c3*ones(n,1);
    d6=zeros(n,1); d6(1:2:n)=c4*ones(n/2,1); d6(2:2:n)=-c2*ones(n/2,1);
    d7=zeros(n,1); d7(1:2:n)=c1*ones(n/2,1); d7(2:2:n)= c5*ones(n/2,1);
    d8=zeros(n,1); d8(1:2:n)=c6*ones(n/2,1); d8(2:2:n)=-c0*ones(n/2,1);
    d9=zeros(n,1);                           d9(2:2:n)= c7*ones(n/2,1);
    
    W=spdiags([d1 d2 d3 d4 d5 d6 d7 d8 d9],-1:7,n,n);
    W(n-5:n-4,1:2)=[c6 c7; c1 -c0];
    W(n-3:n-2,3:4)=[c6 c7; c1 -c0];
    W(n-1:n  ,5:6)=[c6 c7; c1 -c0];
    W(n-3:n-2,1:2)=[c4 c5; c3 -c2];
    W(n-1:n,3:4)=[c4 c5; c3 -c2];
    W(n-1:n,1:2)=[c2 c3; c5 -c4];
elseif m==8      %====================================================
    c0  = 1.8899503329009d-03; c1  =-3.0292051455166d-04;
    c2  =-1.4952258336794d-02; c3  = 3.8087520140604d-03;
    c4  = 4.9137179673481d-02; c5  =-2.7219029916816d-02;
    c6  =-5.1945838107879d-02; c7  = 3.6444189483599d-01;
    c8  = 7.7718575169981d-01; c9  = 4.8135965125924d-01;
    c10 =-6.1273359067914d-02; c11 =-1.4329423835107d-01;
    c12 = 7.6074873252854d-03; c13 = 3.1695087810348d-02;
    c14 =-5.4213233163559d-04; c15 =-3.3824159513597d-03;
    d1=zeros(n,1); d1(1:2:n)=c15*ones(n/2,1);
    d2=zeros(n,1); d2(1:2:n)=c0 *ones(n/2,1); d2(2:2:n)=-c14*ones(n/2,1);
    d3=zeros(n,1); d3(1:2:n)=c13*ones(n/2,1); d3(2:2:n)= c1 *ones(n/2,1);
    d4=zeros(n,1); d4(1:2:n)=c2 *ones(n/2,1); d4(2:2:n)=-c12*ones(n/2,1);
    d5=zeros(n,1); d5(1:2:n)=c11*ones(n/2,1); d5(2:2:n)=c3  *ones(n/2,1);
    d6=zeros(n,1); d6(1:2:n)=c4 *ones(n/2,1); d6(2:2:n)=-c10*ones(n/2,1);
    d7=zeros(n,1); d7(1:2:n)=c9 *ones(n/2,1); d7(2:2:n)= c5 *ones(n/2,1);
    d8=zeros(n,1); d8(1:2:n)=c6 *ones(n/2,1); d8(2:2:n)=-c8 *ones(n/2,1);
    d9=c7*ones(n,1);
    d10=zeros(n,1); d10(1:2:n)=c8 *ones(n/2,1); d10(2:2:n)=-c6*ones(n/2,1);
    d11=zeros(n,1); d11(1:2:n)=c5 *ones(n/2,1); d11(2:2:n)=c9 *ones(n/2,1);
    d12=zeros(n,1); d12(1:2:n)=c10*ones(n/2,1); d12(2:2:n)=-c4*ones(n/2,1);
    d13=zeros(n,1); d13(1:2:n)=c3 *ones(n/2,1); d13(2:2:n)=c11*ones(n/2,1);
    d14=zeros(n,1); d14(1:2:n)=c12*ones(n/2,1); d14(2:2:n)=-c2*ones(n/2,1);
    d15=zeros(n,1); d15(1:2:n)=c1 *ones(n/2,1); d15(2:2:n)=c13*ones(n/2,1);
    d16=zeros(n,1); d16(1:2:n)=c14*ones(n/2,1); d16(2:2:n)=-c0*ones(n/2,1);
    d17=zeros(n,1);                             d17(2:2:n)=c15*ones(n/2,1);
    
    W=spdiags([d1 d2 d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17],-1:15,n,n);
    W=W+spdiags([d3 d4 d5 d6 d7 d8 d9 d10 d11 d12 d13 d14 d15 d16 d17],-n+1:-n+15,n,n);
end
end

function M=fwt(A,m_in,d1,d2,t)
%% fwt.m --- for Daubechies' order m FWT transforms (see also sparse fwts.m)
%% Usage: [1]   y =fwt(x)           FWT for vector x
%%        [2]   y =fwt(x,m)         FWT for vector x (m order)
%%        [3]   B =fwt(A)           FWT for matrix A
%%        [4]   B =fwt(A,m)         FWT for matrix A (m order)
%%        [5]   B =fwt(A,m,d1,d2)   FWT for matrix A with 1-sided options
%%        [6]   B =fwt(A,m,d1,d2,t) for t-levels of FWT for matrix A
%%                                  [dj=1 on,=0 off] with 1-sided options 
%%e.g. A=cz(128,2);
%      B = fwt(A,4, 1, 1, 3); % for 3-levels of Daub_4 (both sides) A=A*W^T 
%      B2= fwt(A,6, 1, 0, 2); % for 2-levels of Daub_6 (left only)  A=W*A
%_______________________________________ Internal Control
if nargin <1, eval(['help ' mfilename]), 
              disp('! At least one parameter needed');  return, end
   [n1 n2]=size(A); n=max(n1,n2); 
if nargin <2, m_in=4; end %------- Set Daub4 (as default)
if nargin <3 %-------------------- Set d1 / d2
   d1=1; d2=1; 
   if n1==1, d1=0; d2=1; end
   if n2==1, d2=0; d1=1; end
end
   t_max=round( log(n)/log(2) );          % n = 2^t_max
if nargin <5 %-------------------- Set t
   t=t_max+1 - ceil( log(m_in)/log(2) );  % Max active levels up to m
end
   Pat='A=WAW'''; if d2==0; Pat='A=WA'; end, if d1==0; Pat='A=AW'''; end
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Internal Control
fprintf('  fwt.m for %d levels of Daub%d to A of %d x %d [%s]\n\t\tDone levs:',...
               t,m_in,n1,n2, Pat);

  [c d]=get_filters(m_in);

for iter = 1:t

  k=floor( n/2^iter ); %k=Half way  % Allow A to have odd dimension=n
  K=2*k;  %Whole block

if d1==1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Left WT
  for col=1:n2
     Sum=zeros(K,1);  % tmp storage for whole col

   for row=1:k        %%%%% Sums for c 

     k_bar=row-1;          % Normal range (2*row-1)   to    (2*row+m-2)
     pos = m_in + 2*k_bar; % i.e. start=2*(row-1) +1  to    2*(row-1)+m
     if pos <= K,
        j1=1; j2=m_in;   j3=1; j4=0; Shift=0;
     else
        j1=1; j2=m_in - (pos-K);
        j3=j2+1;  j4=m_in; Shift=K;  % as (row-Shift) will be wrap around
     end 
     for j=j1:j2  % c_j Sums
       pos = j + 2*k_bar; 
       Sum(row) = Sum(row) + c(j)*A(pos,col); 
     end
     for j=j3:j4  % c_j Sums
       pos = j + 2*k_bar - Shift; 
       Sum(row) = Sum(row) + c(j)*A(pos,col); 
     end

   end %row

   for f_2nd=(1+k):K    %%%%% Differences d
       row = f_2nd - k;

     k_bar=row-1;          % Normal range (2*row-1)   to    (2*row+m-2)
     pos = m_in + 2*k_bar; % i.e. start=2*(row-1) +1  to    2*(row-1)+m
     if pos <= K,
        j1=1; j2=m_in;   j3=1; j4=0; Shift=0;
     else
        j1=1; j2=m_in - (pos-K);
        j3=j2+1;  j4=m_in; Shift=K;  % as (row-Shift) will be wrap around
     end 
     for j=j1:j2  % d_j Difs
       pos = j + 2*k_bar; 
       Sum(f_2nd) = Sum(f_2nd) + d(j)*A(pos,col); 
     end
     for j=j3:j4  % d_j Difs
       pos = j + 2*k_bar - Shift; 
       Sum(f_2nd) = Sum(f_2nd) + d(j)*A(pos,col); 
     end

   end %f_2nd

   A(1:K,col) = Sum;  % Updated - Done Col

   end %col

end % d1=1
if d2==1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Right WT i.e. A=AW^T

  for row=1:n1
     Sum=zeros(1,K);  % tmp storage for whole row

   for col=1:k        %%%%% Sums for c 

     k_bar=col-1;          % Normal range (2*col-1)   to    (2*col+m-2)
     pos = m_in + 2*k_bar; % i.e. start=2*(col-1) +1  to    2*(col-1)+m
     if pos <= K,
        j1=1; j2=m_in;   j3=1; j4=0; Shift=0;
     else
        j1=1; j2=m_in - (pos-K);
        j3=j2+1;  j4=m_in; Shift=K;  % as (col-Shift) will be wrap around
     end 
     for j=j1:j2  % c_j Sums
       pos = j + 2*k_bar; 
       Sum(col) = Sum(col) + c(j)*A(row,pos); 
     end
     for j=j3:j4  % c_j Sums
       pos = j + 2*k_bar - Shift; 
       Sum(col) = Sum(col) + c(j)*A(row,pos); 
     end

   end %col

   for f_2nd=(1+k):K    %%%%% Differences d
       col = f_2nd - k;

     k_bar=col-1;          % Normal range (2*col-1)   to    (2*col+m-2)
     pos = m_in + 2*k_bar; % i.e. start=2*(col-1) +1  to    2*(col-1)+m
     if pos <= K,
        j1=1; j2=m_in;   j3=1; j4=0; Shift=0;
     else
        j1=1; j2=m_in - (pos-K);
        j3=j2+1;  j4=m_in; Shift=K;  % as (col-Shift) will be wrap around
     end 
     for j=j1:j2  % d_j Difs
       pos = j + 2*k_bar; 
       Sum(f_2nd) = Sum(f_2nd) + d(j)*A(row,pos); 
     end
     for j=j3:j4  % d_j Difs
       pos = j + 2*k_bar - Shift; 
       Sum(f_2nd) = Sum(f_2nd) + d(j)*A(row,pos); 
     end

   end %f_2nd

   A(row,1:K) = Sum;  % Updated - Done Col

   end %row

end % id=2

   fprintf(' %d',iter)
end % iter - t levels;
   M = A;    disp(' ')
return
end

%===================================================================
function [c,d] = get_filters(m)
   p = find([2 4 6 8 16] == m);
if isempty(p)
   m,error('Wrong m asked or Requested m not yet implemented');
   return
end,  c=[];d=[]; Sign=-1;
if p==1
    c(1)=sqrt(2)/2;
    c(2)=sqrt(2)/2; 
elseif p==2
    c(1)=sqrt(2)*(1+sqrt(3))/8;
    c(2)=sqrt(2)*(3+sqrt(3))/8;
    c(3)=sqrt(2)*(3-sqrt(3))/8;
    c(4)=sqrt(2)*(1-sqrt(3))/8;
elseif p==3
    c(1)= 3.326705529500826d-01;
    c(2)= 8.068915093110927d-01;
    c(3)= 4.598775021184915d-01;
    c(4)=-1.350110200102546d-01;
    c(5)=-8.544127388202666d-02;
    c(6)= 3.522629188570955d-02;
elseif p==4
    c(1)=-7.5765714789357d-2;
    c(2)=-2.9635527645960d-2; 
    c(3)= 4.9761866763256d-1;
    c(4)= 8.0373875180539d-1;
    c(5)= 2.9785779560561d-1;
    c(6)=-9.9219543576956d-2;
    c(7)=-1.2603967262264d-2;
    c(8)= 3.2223100604078d-2;
elseif p==5
    c(1) = 1.8899503329009d-03;
    c(2) =-3.0292051455166d-04;
    c(3) =-1.4952258336794d-02;
    c(4) = 3.8087520140604d-03;
    c(5) = 4.9137179673481d-02;
    c(6) =-2.7219029916816d-02;
    c(7) =-5.1945838107879d-02;
    c(8) = 3.6444189483599d-01;
    c(9) = 7.7718575169981d-01;
    c(10)= 4.8135965125924d-01;
    c(11)=-6.1273359067914d-02;
    c(12)=-1.4329423835107d-01;
    c(13)= 7.6074873252854d-03;
    c(14)= 3.1695087810348d-02;
    c(15)=-5.4213233163559d-04;
    c(16)=-3.3824159513597d-03;
else
    error('More filters to be added')
end
%========= Set d
    for j=1:m
      d(m+1-j)=Sign*c(j); Sign=-Sign;
    end
return
end

function [Ao,x] = perm0(Ai,d1,d2,lev,levuse,L2)
%% New Wavelet Permutation Ordering   perm0.m
% Usage:  [Ao,x] = perm0(Ai,d1,d2,lev,lev_border,id)
%    or      Ao  = perm0(Ai,d1,d2,lev)      % Exclude the lev_border
%                  for  x = ordering  and lev_border in [1,lev] (id=1 extra)
% If lev-border=lev, perm0 is the same as DWTPer (perm2) !
% If lev-border=0,   perm0 is the same as DWTPerMod (omit J1)
%
%e.g.  n=128,A=cz(n,1);  % n x n CZ matrix %                      from perm2.m 
%        A1=fwt(A,4, 1,1,4); spy(abs(A1)>1/n), figure
%       A0=perm0(A1, 1,1,4); spy(abs(A0)>1/n), figure,
%       A0=perm0(A1, 1,1,4,0); spy(abs(A0)>1/n), figure
%       A0=perm0(A1, 1,1,4,1); spy(abs(A0)>1/n), figure
%       A0=perm0(A1, 1,1,4,2); spy(abs(A0)>1/n) 

if nargin ==1, d1 = 1; d2 = 1; end
if nargin <1, help perm0; return; end
     [m1 m2] = size(Ai); n=max([m1 m2]);
     levmax=floor(log(n)/log(2))-1;
if nargin <4, lev=ceil(log(n)/log(2)); end
if nargin <5, levuse=lev; end
       if m1 ==1, d1=0; end
         if m2 ==1, d2=0; end
n_sav = n; 
exclud = n/2^lev; %n=n-exclud;
    p=speye(n_sav);           x = 1:n; n1=n;
  Pat='A=WAW'''; if d2==0; Pat='A=PA'; end, if d1==0; Pat='A=AP'''; end 
fprintf('  perm0.m for %d levels of A permutation %d x %d [%s]\n\t\tDone levs:',...
               lev,m1,m2, Pat);

for k=lev:-1:1     
    n1 = floor(n_sav/2^k)*2;         % Active length- allow odd n
    n2 = floor(n_sav/2^k);           % Active's Half
         x0 = x(1:n1);   % Part 1 of x
         x2 = x(n1+1:n); % Part 2 of x - fixed
    x1 = [];
    for j=1:n2; x1=[x1 j n2+j]; end
    x1 = x0(x1);   %% Re-order this Part 1
    x = [x1 x2];   %% Re-set x
    fprintf(' %d',k)
end
Ao = Ai;
P=p(:,x);  % for A_new = P'*A_dwt*P
   if d1==1, Ao = P'*Ao; end
   if d2==1, Ao = Ao*P; end
%%% ========================= Possibly A Dense Block at top-left corner :
if nargin<6, L2=0; end
if nargin>4
fprintf('  [used border levels = %d]\n',levuse)
   kb=1:2^(lev-levuse):n; %%%%%%%%%%%%

if L2>0 % Once more recursion (could carry on)
   lb=length(kb); 
   [A1,x]=perm0(Ai(kb,kb),d1,d2,lev,levuse); kb=kb(x);
end
   k1=setdiff(1:n,kb); 
   x=[kb k1]; %=====================================
   P=p(:,x);  % for A_new = P'*A_dwt*P
   if d1==1, Ao = P'*Ao; end,  if d2==1, Ao = Ao*P; end
else,disp(' '),end
end

function [Ao,x] = perm2(Ai,d1,d2,lev)
%% New Wavelet Permutation Ordering  perm2.m
% Usage:  [Ao,x] = Mr_perm(Ai,d1,d2,lev)
%    or      Ao  = Mr_perm(Ai,d1,d2,lev) for x ordering
%
%e.g. A=cz(64,3);  % 64 x 64 CZ matrix %                        from FWT2.m 
%     A1=fwt(A,4, 1,1,3); spyc(A1), figure
%    A2=perm2(A1, 1,1,3); spyc(A2)
 
if nargin ==1,
   d1 = 1; d2 = 1;
end
if nargin <1, help perm2; return; end
if nargin <4, lev=floor(log(n)/log(2)); end
     [m1 m2] = size(Ai); n=max([m1 m2]);
       if m1 ==1, d1=0; end
         if m2 ==1, d2=0; end
    p=speye(n);           x = 1:n; n1=n;
  Pat='A=P''AP'; if d2==0; Pat='A=P''A'; end, if d1==0; Pat='A=AP'; end 
fprintf('  perm2.m for %d levels of A permutation %d x %d [%s]\n\t\tDone levs:',...
               lev,m1,m2, Pat);

for k=lev:-1:1     
    n1 = floor(n/2^k)*2;         % Active length- allow odd n
    n2 = floor(n/2^k);           % Active's Half
         x0 = x(1:n1);   % Part 1 of x
         x2 = x(n1+1:n); % Part 2 of x - fixed
    x1 = [];
    for j=1:n2; x1=[x1 j n2+j]; end
    x1 = x0(x1);   %% Re-order this Part 1
    x = [x1 x2];   %% Re-set x
    fprintf(' %d',k)
end,disp(' ')
   %% P=p(:,x);  % for A_new = P'*A_dwt*P
   %% Ao = Ai;
   %%    if d1==1, Ao = P'*Ao; end
   %%    if d2==1, Ao = Ao*P; end
Ao = Ai;
   if d1==1, Ao = Ao(x,:); end
   if d2==1, Ao = Ao(:,x); end
   %% check = norm(An-Ao,1)
%--------------------------------------------------------------
end
